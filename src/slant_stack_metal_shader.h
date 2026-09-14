// Embedded source: runtime compilation also works with Command Line Tools only.
static const char *slant_stack_metal_source = R"METAL(
#include <metal_stdlib>
using namespace metal;
struct Params {
 uint stations, bins, width, components;
 int half_width, first_bin;
 float step, px0, py0, angular_df, factor;
};
float2 multiply(float2 a, float2 b) {
 return float2(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x);
}
kernel void bias(device const float2 *s [[buffer(0)]],
 device float4 *power [[buffer(1)]], constant Params &p [[buffer(2)]],
 uint k [[thread_position_in_grid]]) {
 if(k>=p.bins) return;
 float4 v=0;
 uint stride=p.stations*p.bins;
 for(uint j=0;j<p.stations;++j) {
  float2 e=s[j*p.bins+k], n=s[stride+j*p.bins+k];
  v.x+=dot(e,e); v.y+=dot(n,n); v.w+=dot(e,n);
  if(p.components==3) {float2 u=s[2*stride+j*p.bins+k]; v.z+=dot(u,u);}
 }
 power[k]=v;
}
kernel void stack(device const float2 *s [[buffer(0)]],
 device const float2 *xy [[buffer(1)]], device const float4 *power [[buffer(2)]],
 device float *out [[buffer(3)]], constant Params &p [[buffer(4)]],
 uint cell [[threadgroup_position_in_grid]], uint lane [[thread_index_in_threadgroup]]) {
 int ix=int(cell/p.width)-p.half_width, iy=int(cell%p.width)-p.half_width;
 float2 slow=float2(p.px0+float(ix)*p.step,p.py0+float(iy)*p.step);
 float2 direction=(ix==0 && iy==0)?float2(0.7071067811865475f):normalize(slow);
 float ex=direction.x,ey=direction.y;
 float3 result=0;
 uint stride=p.stations*p.bins;
 for(uint k=lane;k<p.bins;k+=256) {
  float2 e=0,n=0,u=0;
  float omega=p.angular_df*float(p.first_bin+int(k));
  for(uint j=0;j<p.stations;++j) {
   float angle=-dot(slow,xy[j])*omega;
   float cosine;
   float sine=sincos(angle,cosine);
   float2 rotation=float2(cosine,sine);
   e+=multiply(s[j*p.bins+k],rotation);
   n+=multiply(s[stride+j*p.bins+k],rotation);
   if(p.components==3) u+=multiply(s[2*stride+j*p.bins+k],rotation);
  }
  float4 b=power[k];
  float2 r=ex*e+ey*n,t=ey*e-ex*n;
  result.x+=dot(r,r)-ex*ex*b.x-2*ex*ey*b.w-ey*ey*b.y;
  result.y+=dot(t,t)-ey*ey*b.x+2*ex*ey*b.w-ex*ex*b.y;
  if(p.components==3) result.z+=dot(u,u)-b.z;
 }
 threadgroup float3 partial[256];
 partial[lane]=result;
 threadgroup_barrier(mem_flags::mem_threadgroup);
 for(uint offset=128;offset>0;offset/=2) {
  if(lane<offset) partial[lane]+=partial[lane+offset];
  threadgroup_barrier(mem_flags::mem_threadgroup);
 }
 if(lane==0) {
  uint cells=p.width*p.width;
  out[cell]=partial[0].x*p.factor;
  out[cells+cell]=partial[0].y*p.factor;
  if(p.components==3) out[2*cells+cell]=partial[0].z*p.factor;
 }
}
)METAL";
