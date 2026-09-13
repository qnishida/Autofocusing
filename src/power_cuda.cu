#include "power.h"
#include "cuda_support.h"
#include <algorithm>
#include <cstdint>
#include <mutex>

namespace {
using namespace autofocus_cuda;
struct Config {
  unsigned stations,bins,windows,batches,first_bin,shared_weights,shared_geometry,bias;
  float angular_df;
};
__global__ void phases_kernel(const float2 *xy, const float4 *point, float2 *phase, Config c) {
  unsigned id=blockIdx.x*blockDim.x+threadIdx.x;
  unsigned count=(c.shared_geometry?1:c.batches)*c.stations*c.bins;
  if (id>=count) return;
  unsigned k=id%c.bins, station=(id/c.bins)%c.stations, batch=id/(c.stations*c.bins);
  float4 p=point[batch];
  float ex=cosf(p.y), ey=sinf(p.y), cot=cosf(p.z)/sinf(p.z);
  float2 pos=xy[station];
  float eta=(ex*pos.x+ey*pos.y)/6371.f, zeta=(-ey*pos.x+ex*pos.y)/6371.f;
  float l=(-eta+zeta*zeta*cot/2.f+eta*zeta*zeta*(1.f/6.f+cot*cot/2.f))*6371.f;
  float angle=l*(p.x+p.w*l/2.f)*c.angular_df*float(k+c.first_bin);
  phase[id]=make_float2(cosf(angle),sinf(angle));
}
__global__ void powers_kernel(const float2 *spec, const float *weights,
                               const float2 *phase, float2 *partial, Config c) {
  unsigned id=blockIdx.x*blockDim.x+threadIdx.x, wb=c.windows*c.bins;
  if (id>=c.batches*wb) return;
  unsigned b=id/wb, w=(id%wb)/c.bins, k=id%c.bins;
  unsigned base=(c.shared_weights?0:b)*c.windows*c.stations+w*c.stations;
  unsigned ph=(c.shared_geometry?0:b)*c.stations*c.bins;
  float2 sum=make_float2(0,0); float bias=0;
  for (unsigned s=0;s<c.stations;++s) {
    float2 v=spec[(w*c.stations+s)*c.bins+k], r=phase[ph+s*c.bins+k];
    float weight=weights[base+s];
    sum.x+=(v.x*r.x-v.y*r.y)*weight;
    sum.y+=(v.x*r.y+v.y*r.x)*weight;
    if (c.bias) bias+=(v.x*v.x+v.y*v.y)*weight*weight;
  }
  partial[id]=make_float2(sum.x*sum.x+sum.y*sum.y,bias);
}
__global__ void reduce_kernel(const float2 *partial, float *out, Config c) {
  unsigned t=threadIdx.x, b=blockIdx.x, count=c.windows*c.bins;
  __shared__ float2 sums[256];
  float2 total=make_float2(0,0);
  for (unsigned k=t;k<count;k+=256) { total.x+=partial[b*count+k].x; total.y+=partial[b*count+k].y; }
  sums[t]=total; __syncthreads();
  for (unsigned step=128;step>0;step/=2) {
    if (t<step) { sums[t].x+=sums[t+step].x; sums[t].y+=sums[t+step].y; }
    __syncthreads();
  }
  if (t==0) out[b]=sums[0].x-sums[0].y;
}
struct Backend {
  Stream stream;
  Buffer<float2> spectra,xy,phases,partial;
  Buffer<float4> points;
  Buffer<float> weights,output;
  std::mutex mutex;
};
} // namespace

std::vector<double> cuda_power_batch(const PowerData &d, const std::vector<PowerPoint> &p,
    const std::vector<double> &weights, bool shared_weights, bool shared_geometry, bool bias) {
  const size_t ws=product(d.windows,d.stations), wb=product(d.windows,d.bins);
  const size_t sb=product(d.stations,d.bins), samples=product(ws,d.bins);
  if (!d.windows || !d.stations || !d.bins || p.empty() || p.size()>UINT32_MAX ||
      samples>UINT32_MAX || sb>UINT32_MAX/32 || wb>UINT32_MAX/32 || ws>UINT32_MAX/32 ||
      uint64_t(d.first_bin)+d.bins-1>UINT32_MAX ||
      d.spectra.size()!=samples || d.xy.size()!=product(d.stations,2) ||
      weights.size()!=product(ws,shared_weights?1:p.size()) || !std::isfinite(d.df) || d.df<=0)
    throw std::invalid_argument("CUDA power: invalid dimensions");
  const float angular_df=narrow(2*M_PI*d.df);
  static Backend b;
  std::lock_guard<std::mutex> lock(b.mutex);
  b.stream.activate();
  Drain drain(b.stream.value);
  auto stream=b.stream.value;
  b.spectra.reserve(samples); b.xy.reserve(d.stations);
  for (size_t i=0;i<samples;++i)
    b.spectra.host[i]=make_float2(narrow(d.spectra[i].real()),narrow(d.spectra[i].imag()));
  for (size_t i=0;i<d.stations;++i) b.xy.host[i]=make_float2(narrow(d.xy[2*i]),narrow(d.xy[2*i+1]));
  b.spectra.upload(samples,stream); b.xy.upload(d.stations,stream);
  std::vector<double> result(p.size());
  KernelTimer timer;
  for (size_t first=0;first<p.size();first+=32) {
    const size_t n=std::min(size_t(32),p.size()-first), np=shared_geometry?1:n, nw=shared_weights?1:n;
    b.points.reserve(np); b.weights.reserve(nw*ws);
    b.phases.reserve(np*sb,false); b.partial.reserve(n*wb,false); b.output.reserve(n);
    for (size_t i=0;i<np;++i) {
      const auto &v=p[shared_geometry?0:first+i];
      b.points.host[i]=make_float4(narrow(v.p),narrow(v.theta),narrow(v.delta),narrow(v.curvature));
    }
    for (size_t i=0;i<nw*ws;++i) b.weights.host[i]=narrow(weights[(shared_weights?0:first*ws)+i]);
    b.points.upload(np,stream); b.weights.upload(nw*ws,stream);
    Config c{d.stations,d.bins,d.windows,unsigned(n),d.first_bin,unsigned(shared_weights),
             unsigned(shared_geometry),unsigned(bias),angular_df};
    timer.begin(stream);
    phases_kernel<<<unsigned((np*sb+255)/256),256,0,stream>>>(b.xy.device,b.points.device,b.phases.device,c);
    check(cudaGetLastError(),"launch power phases");
    powers_kernel<<<unsigned((n*wb+255)/256),256,0,stream>>>(b.spectra.device,b.weights.device,b.phases.device,b.partial.device,c);
    check(cudaGetLastError(),"launch power sums");
    reduce_kernel<<<unsigned(n),256,0,stream>>>(b.partial.device,b.output.device,c);
    check(cudaGetLastError(),"launch power reduction");
    timer.stop(stream);
    b.output.download(n,stream);
    check(cudaStreamSynchronize(stream),"complete power batch");
    timer.collect();
    for (size_t i=0;i<n;++i) {
      if (!std::isfinite(b.output.host[i])) throw std::runtime_error("CUDA power: nonfinite objective");
      result[first+i]=double(b.output.host[i]);
    }
  }
  timer.report("power",b.spectra.bytes()+b.xy.bytes()+b.phases.bytes()+b.partial.bytes()+
                         b.points.bytes()+b.weights.bytes()+b.output.bytes());
  return result;
}
