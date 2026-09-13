#pragma once
#include <cuda_runtime.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

// Internal CUDA resources; no CUDA headers escape into the portable interface.
namespace autofocus_cuda {
inline void check(cudaError_t result, const char *operation) {
  if (result != cudaSuccess)
    throw std::runtime_error(std::string("CUDA: ") + operation + ": " +
                             cudaGetErrorString(result));
}
inline size_t product(size_t a, size_t b) {
  if (b && a > std::numeric_limits<size_t>::max() / b)
    throw std::invalid_argument("CUDA: dimension overflow");
  return a * b;
}
inline float narrow(double value) {
  float result = static_cast<float>(value);
  if (!std::isfinite(result))
    throw std::invalid_argument("CUDA: input is nonfinite or outside FP32 range");
  return result;
}
struct Stream {
  cudaStream_t value = nullptr;
  std::string device_name;
  int device = 0;
  Stream() {
    int count = 0;
    check(cudaGetDeviceCount(&count), "enumerate devices");
    if (!count) throw std::runtime_error("CUDA: no accessible GPU device");
    check(cudaGetDevice(&device), "get device");
    cudaDeviceProp properties{};
    check(cudaGetDeviceProperties(&properties, device), "get device properties");
    if (properties.maxThreadsPerBlock < 256)
      throw std::runtime_error("CUDA: backend requires 256 threads per block");
    device_name = properties.name;
    check(cudaStreamCreateWithFlags(&value, cudaStreamNonBlocking), "create stream");
  }
  void activate() const { check(cudaSetDevice(device), "select backend device"); }
  ~Stream() { if (value) { cudaSetDevice(device); cudaStreamDestroy(value); } }
  Stream(const Stream &) = delete;
  Stream &operator=(const Stream &) = delete;
};
// Also drain queued work on exceptions before any host staging storage is reused.
struct Drain {
  cudaStream_t stream;
  explicit Drain(cudaStream_t s) : stream(s) {}
  ~Drain() { cudaStreamSynchronize(stream); }
};
template <class T> struct Buffer {
  T *device = nullptr, *host = nullptr;
  size_t capacity = 0;
  Buffer() = default;
  Buffer(const Buffer &) = delete;
  Buffer &operator=(const Buffer &) = delete;
  ~Buffer() { if (device) cudaFree(device); if (host) cudaFreeHost(host); }
  void reserve(size_t count, bool staging = true) {
    if (count <= capacity) return;
    const size_t bytes = product(count, sizeof(T));
    T *new_device = nullptr, *new_host = nullptr;
    check(cudaMalloc(reinterpret_cast<void **>(&new_device), bytes), "allocate device buffer");
    if (staging) {
      cudaError_t status = cudaMallocHost(reinterpret_cast<void **>(&new_host), bytes);
      if (status != cudaSuccess) { cudaFree(new_device); check(status, "allocate staging buffer"); }
    }
    if (device) cudaFree(device);
    if (host) cudaFreeHost(host);
    device = new_device; host = new_host; capacity = count;
  }
  void upload(size_t count, cudaStream_t stream) {
    check(cudaMemcpyAsync(device, host, product(count, sizeof(T)),
                          cudaMemcpyHostToDevice, stream), "upload buffer");
  }
  void download(size_t count, cudaStream_t stream) {
    check(cudaMemcpyAsync(host, device, product(count, sizeof(T)),
                          cudaMemcpyDeviceToHost, stream), "download buffer");
  }
  size_t bytes() const { return capacity * sizeof(T); }
};
struct KernelTimer {
  cudaEvent_t start = nullptr, end = nullptr;
  double seconds = 0;
  KernelTimer() {
    const char *setting = std::getenv("AUTOFOCUSING_CUDA_PROFILE");
    if (setting && std::string(setting) == "1") {
      check(cudaEventCreate(&start), "create profiling event");
      auto status = cudaEventCreate(&end);
      if (status != cudaSuccess) { cudaEventDestroy(start); check(status, "create profiling event"); }
    }
  }
  ~KernelTimer() { if (start) cudaEventDestroy(start); if (end) cudaEventDestroy(end); }
  void begin(cudaStream_t stream) { if (start) check(cudaEventRecord(start, stream), "start profile"); }
  void stop(cudaStream_t stream) { if (end) check(cudaEventRecord(end, stream), "stop profile"); }
  void collect() {
    if (!end) return;
    float milliseconds = 0;
    check(cudaEventElapsedTime(&milliseconds, start, end), "read profile");
    seconds += milliseconds / 1000.;
  }
  void report(const char *stage, size_t bytes) const {
    if (end) std::cerr << "#CUDA_PROFILE stage=" << stage << " kernel_s=" << seconds
                       << " reserved_device_bytes=" << bytes << '\n';
  }
};
} // namespace autofocus_cuda
