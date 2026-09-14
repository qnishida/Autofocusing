#include "slant_stack.h"
#include "cuda_support.h"
#include <cstdint>
#include <mutex>

namespace {
using namespace autofocus_cuda;
struct Params {
  unsigned stations, bins, width, components;
  int half_width, first_bin;
  float step, px0, py0, angular_df, factor;
};
__device__ float norm2(float2 v) { return v.x*v.x + v.y*v.y; }
__device__ float2 multiply(float2 a, float2 b) {
  return make_float2(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x);
}
__global__ void bias_kernel(const float2 *s, float4 *power, Params p) {
  unsigned k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= p.bins) return;
  float4 v = make_float4(0,0,0,0);
  unsigned stride = p.stations*p.bins;
  for (unsigned j=0; j<p.stations; ++j) {
    float2 e=s[j*p.bins+k], n=s[stride+j*p.bins+k];
    v.x+=norm2(e); v.y+=norm2(n); v.w+=e.x*n.x+e.y*n.y;
    if (p.components==3) v.z+=norm2(s[2*stride+j*p.bins+k]);
  }
  power[k]=v;
}
__global__ void stack_kernel(const float2 *s, const float2 *xy,
                              const float4 *power, float *out, Params p) {
  unsigned cell=blockIdx.x, lane=threadIdx.x;
  int ix=int(cell/p.width)-p.half_width, iy=int(cell%p.width)-p.half_width;
  float px=p.px0+float(ix)*p.step, py=p.py0+float(iy)*p.step;
  float ex=0.7071067811865475f, ey=ex;
  if (ix!=0 || iy!=0) { float length=sqrtf(px*px+py*py); ex=px/length; ey=py/length; }
  float3 result=make_float3(0,0,0);
  unsigned stride=p.stations*p.bins;
  for (unsigned k=lane; k<p.bins; k+=256) {
    float2 e=make_float2(0,0), n=e, u=e;
    float omega=p.angular_df*float(p.first_bin+int(k));
    for (unsigned j=0; j<p.stations; ++j) {
      float angle=-(px*xy[j].x+py*xy[j].y)*omega;
      float sine, cosine;
      sincosf(angle,&sine,&cosine);
      float2 rotation=make_float2(cosine,sine);
      float2 v=multiply(s[j*p.bins+k],rotation); e.x+=v.x; e.y+=v.y;
      v=multiply(s[stride+j*p.bins+k],rotation); n.x+=v.x; n.y+=v.y;
      if (p.components==3) { v=multiply(s[2*stride+j*p.bins+k],rotation); u.x+=v.x; u.y+=v.y; }
    }
    float4 b=power[k];
    float2 r=make_float2(ex*e.x+ey*n.x,ex*e.y+ey*n.y);
    float2 t=make_float2(ey*e.x-ex*n.x,ey*e.y-ex*n.y);
    result.x+=norm2(r)-ex*ex*b.x-2*ex*ey*b.w-ey*ey*b.y;
    result.y+=norm2(t)-ey*ey*b.x+2*ex*ey*b.w-ex*ex*b.y;
    if (p.components==3) result.z+=norm2(u)-b.z;
  }
  __shared__ float3 partial[256];
  partial[lane]=result;
  __syncthreads();
  for (unsigned offset=128; offset>0; offset/=2) {
    if (lane<offset) {
      partial[lane].x+=partial[lane+offset].x;
      partial[lane].y+=partial[lane+offset].y;
      partial[lane].z+=partial[lane+offset].z;
    }
    __syncthreads();
  }
  if (lane==0) {
    unsigned cells=p.width*p.width;
    out[cell]=partial[0].x*p.factor;
    out[cells+cell]=partial[0].y*p.factor;
    if (p.components==3) out[2*cells+cell]=partial[0].z*p.factor;
  }
}
struct Backend {
  Stream stream;
  Buffer<float2> spectra, coords;
  Buffer<float4> power;
  Buffer<float> output;
  std::mutex mutex;
};
Backend &backend() { static Backend b; return b; }
} // namespace

std::string slant_stack_cuda_device() { return backend().stream.device_name; }

void slant_stack_cuda(const std::complex<double> *input, size_t capacity,
                      size_t count, const double *dx, const double *dy,
                      const SlantStackGrid &g, double *rtu) {
  const int64_t nf64=int64_t(g.last_bin)-g.first_bin+1;
  if (count<2 || nf64<=0) return;
  if (count>capacity || g.half_width<0 || g.half_width>16383 ||
      g.first_bin<0 || nf64>INT32_MAX ||
      count>UINT32_MAX/size_t(nf64)/3 ||
      !std::isfinite(g.df) || g.df<=0 || !std::isfinite(g.step) || g.step<=0 ||
      !input || !dx || !dy || !rtu)
    throw std::invalid_argument("CUDA: unsupported slant-stack dimensions or grid");
  const unsigned nf=unsigned(nf64), components=g.horizontal_only?2:3;
  const size_t width=2*g.half_width+1, cells=width*width;
  product(product(capacity,nf),3); // Check host indexing before dereferencing.
  Params p{unsigned(count),nf,unsigned(width),components,g.half_width,g.first_bin,
           narrow(g.step),narrow(g.px0),narrow(g.py0),narrow(2*M_PI*g.df),
           narrow(1./count/(count-1)/nf)};
  auto &b=backend();
  std::lock_guard<std::mutex> lock(b.mutex);
  b.stream.activate();
  Drain drain(b.stream.value);
  b.spectra.reserve(components*count*nf); b.coords.reserve(count);
  b.power.reserve(nf,false); b.output.reserve(components*cells);
  for (unsigned c=0;c<components;++c)
    for (size_t j=0;j<count;++j)
      for (unsigned k=0;k<nf;++k) {
        auto v=input[(c*capacity+j)*nf+k];
        b.spectra.host[(c*count+j)*nf+k]=make_float2(narrow(v.real()),narrow(v.imag()));
      }
  for (size_t j=0;j<count;++j) b.coords.host[j]=make_float2(narrow(dx[j]),narrow(dy[j]));
  auto stream=b.stream.value;
  b.spectra.upload(components*count*nf,stream); b.coords.upload(count,stream);
  KernelTimer timer;
  timer.begin(stream);
  bias_kernel<<<(nf+255)/256,256,0,stream>>>(b.spectra.device,b.power.device,p);
  check(cudaGetLastError(),"launch slant-stack bias");
  stack_kernel<<<unsigned(cells),256,0,stream>>>(b.spectra.device,b.coords.device,b.power.device,b.output.device,p);
  check(cudaGetLastError(),"launch slant stack");
  timer.stop(stream);
  b.output.download(components*cells,stream);
  check(cudaStreamSynchronize(stream),"complete slant stack");
  for (size_t i=0;i<components*cells;++i)
    if (!std::isfinite(b.output.host[i])) throw std::runtime_error("CUDA: nonfinite slant-stack output");
  for (size_t i=0;i<components*cells;++i) rtu[i]+=double(b.output.host[i]);
  timer.collect();
  timer.report("stack",b.spectra.bytes()+b.coords.bytes()+b.power.bytes()+b.output.bytes());
}
