#include "metal_power.h"
#include "metal_power_shader.h"
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <stdexcept>
namespace {
struct F2 { float x,y; };
struct F4 { float x,y,z,w; };
struct Config { uint32_t stations,bins,windows,batches,first_bin,shared_weights,shared_geometry,bias; float angular_df; };
static_assert(sizeof(Config)==36,"Metal power config layout");
struct Backend {
  id<MTLDevice> device;
  id<MTLCommandQueue> queue;
  id<MTLComputePipelineState> phase, power, reduce;
  id<MTLBuffer> spectra, xy, points, weights, phases, partial, output;
  std::mutex mutex;
  Backend() {
    device=MTLCreateSystemDefaultDevice();
    if(!device) throw std::runtime_error("Metal power: no GPU");
    queue=[device newCommandQueue];
    MTLCompileOptions *options=[MTLCompileOptions new];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    options.fastMathEnabled=NO;
#pragma clang diagnostic pop
    NSError *error=nil;
    auto lib=[device newLibraryWithSource:[NSString stringWithUTF8String:metal_power_source] options:options error:&error];
    if(!lib) throw std::runtime_error(std::string("Metal power shader: ")+error.description.UTF8String);
    phase=[device newComputePipelineStateWithFunction:[lib newFunctionWithName:@"phases"] error:&error];
    power=[device newComputePipelineStateWithFunction:[lib newFunctionWithName:@"powers"] error:&error];
    reduce=[device newComputePipelineStateWithFunction:[lib newFunctionWithName:@"reduce_power"] error:&error];
    if(!queue || !phase || !power || !reduce || reduce.maxTotalThreadsPerThreadgroup<256)
      throw std::runtime_error("Metal power: pipeline creation failed");
  }
  void reserve(id<MTLBuffer> __strong &buffer,size_t bytes) {
    if(!buffer || buffer.length<bytes) buffer=[device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
    if(!buffer) throw std::runtime_error("Metal power: allocation failed");
  }
};
}
std::vector<double> metal_power_batch(const PowerData &d,const std::vector<PowerPoint> &p,
    const std::vector<double> &weights,bool shared_weights,bool shared_geometry,bool bias) {
  const size_t ws=size_t(d.windows)*d.stations, wb=size_t(d.windows)*d.bins;
  if(!d.windows || !d.stations || !d.bins || p.empty() || p.size()>UINT32_MAX ||
     ws*d.bins>UINT32_MAX || size_t(d.stations)*d.bins>UINT32_MAX/32 ||
     wb>UINT32_MAX/32 || d.spectra.size()!=ws*d.bins || d.xy.size()!=size_t(d.stations)*2 ||
     weights.size()!=ws*(shared_weights?1:p.size()) || !std::isfinite(d.df) || d.df<=0)
    throw std::invalid_argument("Metal power: invalid dimensions");
  @autoreleasepool {
    static Backend b;
    std::lock_guard<std::mutex> lock(b.mutex);
    b.reserve(b.spectra,d.spectra.size()*sizeof(F2)); b.reserve(b.xy,d.stations*sizeof(F2));
    auto *s=static_cast<F2 *>(b.spectra.contents);
    for(size_t i=0;i<d.spectra.size();++i) s[i]={float(d.spectra[i].real()),float(d.spectra[i].imag())};
    auto *xy=static_cast<F2 *>(b.xy.contents);
    for(size_t i=0;i<d.stations;++i) xy[i]={float(d.xy[2*i]),float(d.xy[2*i+1])};
    std::vector<double> result(p.size());
    // Bound temporary phase/weight memory even for expanded search grids.
    for(size_t first=0;first<p.size();first+=32) {
      size_t n=std::min(size_t(32),p.size()-first), np=shared_geometry?1:n, nw=shared_weights?1:n;
      b.reserve(b.points,np*sizeof(F4)); b.reserve(b.weights,nw*ws*sizeof(float));
      b.reserve(b.phases,np*d.stations*d.bins*sizeof(F2));
      b.reserve(b.partial,n*wb*sizeof(F2)); b.reserve(b.output,n*sizeof(float));
      auto *point=static_cast<F4 *>(b.points.contents);
      for(size_t i=0;i<np;++i) { const auto &v=p[shared_geometry?0:first+i];point[i]={float(v.p),float(v.theta),float(v.delta),float(v.curvature)}; }
      auto *w=static_cast<float *>(b.weights.contents);
      for(size_t i=0;i<nw*ws;++i) w[i]=float(weights[(shared_weights?0:first*ws)+i]);
      Config c{d.stations,d.bins,d.windows,uint32_t(n),d.first_bin,uint32_t(shared_weights),uint32_t(shared_geometry),uint32_t(bias),float(2*M_PI*d.df)};
      auto command=[b.queue commandBuffer];
      auto encoder=[command computeCommandEncoder];
      if(!command || !encoder) throw std::runtime_error("Metal power: command creation failed");
      [encoder setComputePipelineState:b.phase];
      [encoder setBuffer:b.xy offset:0 atIndex:0]; [encoder setBuffer:b.points offset:0 atIndex:1];
      [encoder setBuffer:b.phases offset:0 atIndex:2]; [encoder setBytes:&c length:sizeof(c) atIndex:3];
      [encoder dispatchThreadgroups:MTLSizeMake((np*d.stations*d.bins+255)/256,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
      [encoder endEncoding];
      encoder=[command computeCommandEncoder];
      [encoder setComputePipelineState:b.power];
      [encoder setBuffer:b.spectra offset:0 atIndex:0]; [encoder setBuffer:b.weights offset:0 atIndex:1];
      [encoder setBuffer:b.phases offset:0 atIndex:2]; [encoder setBuffer:b.partial offset:0 atIndex:3];
      [encoder setBytes:&c length:sizeof(c) atIndex:4];
      [encoder dispatchThreadgroups:MTLSizeMake((n*wb+255)/256,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
      [encoder endEncoding];
      encoder=[command computeCommandEncoder];
      [encoder setComputePipelineState:b.reduce]; [encoder setBuffer:b.partial offset:0 atIndex:0];
      [encoder setBuffer:b.output offset:0 atIndex:1]; [encoder setBytes:&c length:sizeof(c) atIndex:2];
      [encoder dispatchThreadgroups:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
      [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
      if(command.status!=MTLCommandBufferStatusCompleted) throw std::runtime_error("Metal power: execution failed");
      auto *out=static_cast<float *>(b.output.contents);
      for(size_t i=0;i<n;++i) {
        if(!std::isfinite(out[i])) throw std::runtime_error("Metal power: nonfinite objective");
        result[first+i]=out[i];
      }
    }
    return result;
  }
}
