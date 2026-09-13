#pragma once
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <omp.h>

// Inclusive call times, not additive pipeline wall time. Calls made by outer
// OpenMP workers are reported separately from calls made by the main thread.
namespace cpu_profile {
enum Stage { objective, hessian, rotation, fitting, count };
struct Counters {
  bool enabled = std::getenv("AUTOFOCUSING_PROFILE") != nullptr;
  std::atomic<unsigned long long> calls[count][2]{};
  std::atomic<unsigned long long> ns[count][2]{};
  ~Counters() {
    if (!enabled) return;
    const char *names[] = {"objective", "hessian", "rotation", "fitting"};
    for (int s = 0; s < count; ++s) for (int p = 0; p < 2; ++p)
      std::cerr << "#CPU_PROFILE stage=" << names[s]
                << " context=" << (p ? "outer_parallel" : "serial_caller")
                << " calls=" << calls[s][p].load()
                << " inclusive_s=" << ns[s][p].load() * 1e-9 << '\n';
  }
};
inline Counters &counters() { static Counters result; return result; }
struct Timer {
  Stage stage;
  bool enabled;
  int parallel = 0;
  std::chrono::steady_clock::time_point start;
  explicit Timer(Stage s) : stage(s), enabled(counters().enabled) {
    if (enabled) {
      parallel = omp_in_parallel() ? 1 : 0;
      start = std::chrono::steady_clock::now();
    }
  }
  ~Timer() {
    if (!enabled) return;
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start).count();
    counters().ns[stage][parallel].fetch_add(elapsed, std::memory_order_relaxed);
    counters().calls[stage][parallel].fetch_add(1, std::memory_order_relaxed);
  }
};
}
