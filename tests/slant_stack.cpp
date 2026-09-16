#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main
#include "slant_stack_reference.h"
#include <iomanip>
#include <random>
#include <stdexcept>

static void run_case(int stations, int half_width, int threads, bool horizontal,
                     int repeats, bool benchmark, bool offset = false, bool primary = false) {
  omp_set_dynamic(0); omp_set_num_threads(threads);
  STATION::if1=primary?51:102; STATION::if2=primary?102:256; STATION::df=1./1024;
  ipmax=half_width; dp=.005; px0=offset?.007:0; py0=offset?-.012:0;
  const int capacity=stations+3, nf=STATION::if2-STATION::if1+1;
  array3c spectra(boost::extents[3][capacity][range3c(STATION::if1,STATION::if2+1)]);
  dvector dx(capacity),dy(capacity);
  std::mt19937 rng(1837); std::uniform_real_distribution<double> noise(-1,1);
  for(int i=0;i<capacity;++i) {
    dx[i]=noise(rng)*1500; dy[i]=noise(rng)*1500;
    for(int k=STATION::if1;k<=STATION::if2;++k)
      for(int c=0;c<3;++c)
        spectra[c][i][k]=(horizontal && c==2)?std::complex<double>(0.,0.):
            std::complex<double>(noise(rng)*1e-8,noise(rng)*1e-8);
  }
  array3d old(boost::extents[3][range3d(-ipmax,ipmax+1)][range3d(-ipmax,ipmax+1)]);
  array3d current(boost::extents[3][range3d(-ipmax,ipmax+1)][range3d(-ipmax,ipmax+1)]);
  std::fill_n(old.data(),old.num_elements(),1e-20);
  std::fill_n(current.data(),current.num_elements(),1e-20);
  const SlantStackGrid grid{ipmax,dp,px0,py0,STATION::if1,STATION::if2,STATION::df,horizontal};
  double old_seconds=0, new_seconds=0;
  // Repeated accumulation tests preservation of pre-existing segment values.
  for(int rep=0;rep<repeats;++rep) {
    auto start=std::chrono::steady_clock::now();
    slant_stack_reference(spectra,dx,dy,stations,old);
    old_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    start=std::chrono::steady_clock::now();
    slant_stack_cpu(spectra.data(),capacity,stations,&dx[0],&dy[0],grid,current.data());
    new_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  }
  double scale=0,error=0;
  for(std::size_t i=0;i<old.num_elements();++i) {
    scale=std::max(scale,std::abs(old.data()[i]));
    error=std::max(error,std::abs(current.data()[i]-old.data()[i]));
    if(!std::isfinite(current.data()[i])) throw std::runtime_error("nonfinite slant stack");
  }
  if(error>scale*1e-12+1e-30) throw std::runtime_error("slant stack differs from legacy reference");
  if(!benchmark) {
    auto before=std::vector<double>(current.data(),current.data()+current.num_elements());
    slant_stack_cpu(spectra.data(),capacity,1,&dx[0],&dy[0],grid,current.data());
    if(!std::equal(before.begin(),before.end(),current.data())) throw std::runtime_error("one-station guard");
  }
  std::cout << std::setprecision(9) << "mode=" << (horizontal?"horizontal":"3c")
            << " stations=" << stations << " grid=" << 2*ipmax+1 << " bins=" << nf
            << " threads=" << threads << " repeats=" << repeats
            << " reference_s=" << old_seconds/repeats << " cpu_s=" << new_seconds/repeats
            << " speedup=" << old_seconds/new_seconds << " max_abs_error=" << error << '\n';
}
int main(int argc,char **argv) {
  try {
    if(argc==6 && std::string(argv[1])=="--benchmark") {
      for(bool horizontal:{true,false})
        run_case(std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]),horizontal,std::stoi(argv[5]),true);
    } else {
      for(bool horizontal:{true,false})
        for(int threads:{1,4})
          for(bool offset:{false,true})
            for(bool primary:{false,true}) run_case(17,6,threads,horizontal,2,false,offset,primary);
    }
  } catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
