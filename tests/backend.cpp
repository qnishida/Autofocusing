#include "../src/slant_stack.h"
#include "../src/power.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

static void require(bool condition) {
  if (!condition) throw std::runtime_error("Backend configuration regression");
}
template <class F> static void fails(F action, const char *message) {
  try { action(); }
  catch (const std::exception &e) { require(std::string(e.what()).find(message)!=std::string::npos); return; }
  throw std::runtime_error("Expected explicit backend failure");
}
int main() try {
  unsetenv("AUTOFOCUSING_BACKEND"); unsetenv("AUTOFOCUSING_GPU_POWER"); unsetenv("AUTOFOCUSING_METAL_POWER");
  require(slant_stack_backend()==SlantStackBackend::Cpu && gpu_power_mode()=="off");
  setenv("AUTOFOCUSING_GPU_POWER","all",1);
  require(!gpu_power_enabled("bootstrap",true));
  require(!gpu_power_enabled("bootstrap",false));
  for (const char *name : {"metal","cuda"}) {
    const bool cuda=std::string(name)=="cuda";
    setenv("AUTOFOCUSING_BACKEND",name,1);
    require(gpu_power_enabled("bootstrap",true) && gpu_power_enabled("grid",true));
    require(gpu_power_enabled("bootstrap",false)==cuda);
    require(!gpu_power_enabled("grid",false));
    setenv("AUTOFOCUSING_GPU_POWER","bootstrap",1);
    require(gpu_power_enabled("bootstrap",true) && !gpu_power_enabled("grid",true));
    require(gpu_power_enabled("bootstrap",false)==cuda && !gpu_power_enabled("grid",false));
    setenv("AUTOFOCUSING_GPU_POWER","grid",1);
    require(!gpu_power_enabled("bootstrap",true) && gpu_power_enabled("grid",true));
    require(!gpu_power_enabled("bootstrap",false) && !gpu_power_enabled("grid",false));
    setenv("AUTOFOCUSING_GPU_POWER","off",1);
    require(!gpu_power_enabled("grid",true));
    require(!gpu_power_enabled("bootstrap",false));
    setenv("AUTOFOCUSING_GPU_POWER","all",1);
  }
  unsetenv("AUTOFOCUSING_GPU_POWER");
  setenv("AUTOFOCUSING_METAL_POWER","all",1);
  require(gpu_power_mode()=="off"); // CUDA ignores the legacy Metal setting.
  setenv("AUTOFOCUSING_BACKEND","metal",1);
  require(gpu_power_mode()=="all");
  PowerData fp64{}; fp64.double_precision=true;
  fails([&]{gpu_power_batch(fp64,{}, {},true,true,true);},"FP64 power evaluation requires CUDA");
  setenv("AUTOFOCUSING_GPU_POWER","off",1);
  require(!gpu_power_enabled("grid",true)); // New setting wins, including off.
  setenv("AUTOFOCUSING_METAL_POWER","invalid",1);
  require(gpu_power_mode()=="off");
  unsetenv("AUTOFOCUSING_GPU_POWER");
  fails([]{gpu_power_mode();},"AUTOFOCUSING_METAL_POWER");
  setenv("AUTOFOCUSING_GPU_POWER","invalid",1);
  fails([]{gpu_power_mode();},"AUTOFOCUSING_GPU_POWER");
  setenv("AUTOFOCUSING_BACKEND","invalid",1);
  fails([]{slant_stack_backend();},"AUTOFOCUSING_BACKEND");
#if !AUTOFOCUSING_TEST_HAVE_CUDA
  setenv("AUTOFOCUSING_BACKEND","cuda",1);
  fails([]{slant_stack_backend_name();},"DELTAP_ENABLE_CUDA=ON");
  fails([]{cuda_power_batch(PowerData{}, {}, {}, false, false, false);},"not compiled");
#endif
#if !AUTOFOCUSING_TEST_HAVE_METAL
  setenv("AUTOFOCUSING_BACKEND","metal",1);
  fails([]{slant_stack_backend_name();},"DELTAP_ENABLE_METAL=ON");
#endif
  std::cout << "PASS: backend selection, power modes, compatibility and unavailable builds\n";
} catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
