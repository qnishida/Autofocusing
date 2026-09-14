#pragma once
static const char *metal_power_source = R"METAL(
#include <metal_stdlib>
using namespace metal;
struct Config { uint stations, bins, windows, batches, first_bin, shared_weights, shared_geometry, bias; float angular_df; };
kernel void phases(device const float2 *xy [[buffer(0)]], device const float4 *point [[buffer(1)]],
                   device float2 *phase [[buffer(2)]], constant Config &c [[buffer(3)]], uint id [[thread_position_in_grid]]) {
  uint count=(c.shared_geometry?1:c.batches)*c.stations*c.bins;
  if(id>=count) return;
  uint k=id%c.bins, station=(id/c.bins)%c.stations, batch=id/(c.stations*c.bins);
  float4 p=point[batch]; float ex=cos(p.y), ey=sin(p.y), cot=cos(p.z)/sin(p.z);
  float2 pos=xy[station]; float eta=(ex*pos.x+ey*pos.y)/6371.f, zeta=(-ey*pos.x+ex*pos.y)/6371.f;
  float l=(-eta+zeta*zeta*cot/2.f+eta*zeta*zeta*(1.f/6.f+cot*cot/2.f))*6371.f;
  float angle=l*(p.x+p.w*l/2.f)*c.angular_df*float(k+c.first_bin);
  phase[id]=float2(cos(angle),sin(angle));
}
kernel void powers(device const float2 *spec [[buffer(0)]], device const float *weights [[buffer(1)]],
                   device const float2 *phase [[buffer(2)]], device float2 *partial [[buffer(3)]],
                   constant Config &c [[buffer(4)]], uint id [[thread_position_in_grid]]) {
  uint wb=c.windows*c.bins;
  if(id>=c.batches*wb) return;
  uint b=id/wb, w=(id%wb)/c.bins, k=id%c.bins;
  uint base=(c.shared_weights?0:b)*c.windows*c.stations+w*c.stations;
  uint ph=(c.shared_geometry?0:b)*c.stations*c.bins;
  float2 sum=0; float bias=0;
  for(uint s=0;s<c.stations;++s) {
    float2 v=spec[(w*c.stations+s)*c.bins+k], r=phase[ph+s*c.bins+k];
    float weight=weights[base+s];
    sum+=float2(v.x*r.x-v.y*r.y,v.x*r.y+v.y*r.x)*weight;
    if(c.bias) bias+=dot(v,v)*weight*weight;
  }
  partial[id]=float2(dot(sum,sum),bias);
}
kernel void reduce_power(device const float2 *partial [[buffer(0)]], device float *out [[buffer(1)]],
                         constant Config &c [[buffer(2)]], uint t [[thread_index_in_threadgroup]],
                         uint b [[threadgroup_position_in_grid]]) {
  threadgroup float2 sums[256];
  uint count=c.windows*c.bins; float2 total=0;
  for(uint k=t;k<count;k+=256) total+=partial[b*count+k];
  sums[t]=total; threadgroup_barrier(mem_flags::mem_threadgroup);
  for(uint step=128;step>0;step/=2) {
    if(t<step) sums[t]+=sums[t+step];
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }
  if(t==0) out[b]=sums[0].x-sums[0].y;
}
)METAL";
