#include "slant_stack.h"
#include "slant_stack_metal_shader.h"
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace {
struct Float2 {
  float x, y;
};
struct Params {
  uint32_t stations, bins, width, components;
  int32_t half_width, first_bin;
  float step, px0, py0, angular_df, factor;
};
static_assert(sizeof(Params) == 44, "Metal parameter layout");
std::runtime_error failure(const char *context, NSError *error = nil) {
  return std::runtime_error(
      std::string("Metal: ") + context +
      (error ? ": " + std::string(error.description.UTF8String) : ""));
}
struct Backend {
  id<MTLDevice> device;
  id<MTLCommandQueue> queue;
  id<MTLComputePipelineState> bias, stack;
  id<MTLBuffer> spectra, coords, power, output;
  std::mutex mutex;
  Backend() {
    @autoreleasepool {
      device = MTLCreateSystemDefaultDevice();
      if (!device)
        throw failure("no accessible GPU device");
      queue = [device newCommandQueue];
      if (!queue)
        throw failure("cannot create command queue");
      NSError *error = nil;
      MTLCompileOptions *options = [MTLCompileOptions new];
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
      if (@available(macOS 15.0, *)) {
        options.mathMode = MTLMathModeSafe;
      } else {
#endif
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
        options.fastMathEnabled = NO;
#pragma clang diagnostic pop
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
      }
#endif
      id<MTLLibrary> library =
          [device newLibraryWithSource:
                      [NSString stringWithUTF8String:slant_stack_metal_source]
                               options:options
                                 error:&error];
      if (!library)
        throw failure("shader compilation failed", error);
      bias = [device
          newComputePipelineStateWithFunction:[library
                                                  newFunctionWithName:@"bias"]
                                        error:&error];
      if (!bias)
        throw failure("bias pipeline failed", error);
      stack = [device
          newComputePipelineStateWithFunction:[library
                                                  newFunctionWithName:@"stack"]
                                        error:&error];
      if (!stack || stack.maxTotalThreadsPerThreadgroup < 256)
        throw failure("stack pipeline requires 256 threads", error);
    }
  }
  void reserve(id<MTLBuffer> __strong &buffer, std::size_t bytes) {
    if (!buffer || buffer.length < bytes)
      buffer = [device newBufferWithLength:bytes
                                   options:MTLResourceStorageModeShared];
    if (!buffer)
      throw failure("buffer allocation failed");
  }
};
Backend &backend() {
  static Backend state;
  return state;
}
} // namespace
std::string slant_stack_metal_device() {
  @autoreleasepool {
    return backend().device.name.UTF8String;
  }
}
void slant_stack_metal(const std::complex<double> *input, std::size_t capacity,
                       std::size_t count, const double *dx, const double *dy,
                       const SlantStackGrid &g, double *rtu) {
  const int nf = g.last_bin - g.first_bin + 1;
  if (count < 2 || nf <= 0)
    return;
  if (count > capacity || g.half_width < 0 || g.half_width > 16383 ||
      count > std::numeric_limits<uint32_t>::max() / std::size_t(nf) / 3)
    throw std::invalid_argument("Metal: unsupported slant-stack dimensions");
  @autoreleasepool {
    auto &b = backend();
    std::lock_guard<std::mutex> lock(b.mutex);
    const std::size_t width = 2 * g.half_width + 1, cells = width * width;
    const unsigned components = g.horizontal_only ? 2 : 3;
    b.reserve(b.spectra, components * count * nf * sizeof(Float2));
    b.reserve(b.coords, count * sizeof(Float2));
    b.reserve(b.power, nf * 4 * sizeof(float));
    b.reserve(b.output, components * cells * sizeof(float));
    auto *s = static_cast<Float2 *>(b.spectra.contents);
    for (unsigned c = 0; c < components; ++c)
      for (std::size_t j = 0; j < count; ++j)
        for (int k = 0; k < nf; ++k) {
          auto v = input[(c * capacity + j) * nf + k];
          s[(c * count + j) * nf + k] = {float(v.real()), float(v.imag())};
        }
    auto *xy = static_cast<Float2 *>(b.coords.contents);
    for (std::size_t j = 0; j < count; ++j)
      xy[j] = {float(dx[j]), float(dy[j])};
    Params p{uint32_t(count),
             uint32_t(nf),
             uint32_t(width),
             components,
             g.half_width,
             g.first_bin,
             float(g.step),
             float(g.px0),
             float(g.py0),
             float(2 * M_PI * g.df),
             float(1. / count / (count - 1) / nf)};
    id<MTLCommandBuffer> command = [b.queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    if (!command || !encoder)
      throw failure("cannot create command encoder");
    [encoder setComputePipelineState:b.bias];
    [encoder setBuffer:b.spectra offset:0 atIndex:0];
    [encoder setBuffer:b.power offset:0 atIndex:1];
    [encoder setBytes:&p length:sizeof(p) atIndex:2];
    [encoder dispatchThreadgroups:MTLSizeMake((nf + 255) / 256, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
    [encoder endEncoding];
    encoder = [command computeCommandEncoder];
    if (!encoder)
      throw failure("cannot create stack encoder");
    [encoder setComputePipelineState:b.stack];
    [encoder setBuffer:b.spectra offset:0 atIndex:0];
    [encoder setBuffer:b.coords offset:0 atIndex:1];
    [encoder setBuffer:b.power offset:0 atIndex:2];
    [encoder setBuffer:b.output offset:0 atIndex:3];
    [encoder setBytes:&p length:sizeof(p) atIndex:4];
    [encoder dispatchThreadgroups:MTLSizeMake(cells, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted)
      throw failure("GPU execution failed", command.error);
    const auto *out = static_cast<const float *>(b.output.contents);
    for (std::size_t i = 0; i < components * cells; ++i)
      if (!std::isfinite(out[i]))
        throw failure("nonfinite GPU output");
    for (std::size_t i = 0; i < components * cells; ++i)
      rtu[i] += double(out[i]);
  }
}
