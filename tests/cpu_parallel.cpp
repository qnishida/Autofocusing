// Exercise production internals against an independently frozen CPU baseline.
#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main
#include "cpu_parallel_reference.h"
#include <iomanip>
#include <random>

namespace {
double max_abs = 0, max_scaled = 0;
size_t checked = 0, unequal = 0;
void check(double actual, double expected) {
  const double error = std::abs(actual - expected);
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      error > 1e-12 * std::abs(expected) + 1e-30)
    throw std::runtime_error("CPU parallel regression differs from frozen reference");
  max_abs = std::max(max_abs, error);
  max_scaled = std::max(max_scaled, error / (std::abs(expected) + 1e-30));
  ++checked; unequal += actual != expected;
}
void check_complex(std::complex<double> actual, std::complex<double> expected) {
  check(actual.real(), expected.real()); check(actual.imag(), expected.imag());
}
struct Fixture {
  int windows, stations;
  array3c spec;
  array2d weights;
  dvector x, y;
  PARAM point{};
  Fixture(int w, int n, int bins, bool coherent, bool zero = false)
      : windows(w), stations(n),
        spec(boost::extents[w][n][range3c(102, 102 + bins)]),
        weights(boost::extents[w][n]), x(n), y(n) {
    STATION::df = 1. / 1024; STATION::if1 = 102;
    STATION::if2 = 101 + bins; STATION::nfreq = 102 + bins;
    point.p = .06; point.θ = .73; point.Δ = .85; point.dp_Δ = -2e-6;
    std::mt19937 gen(1837); std::normal_distribution<double> normal;
    std::uniform_real_distribution<double> coord(-1500, 1500);
    for (int i = 0; i < n; ++i) {
      x[i] = coord(gen); y[i] = coord(gen);
      const double ex = cos(point.θ), ey = sin(point.θ), cot = cos(point.Δ) / sin(point.Δ);
      const double eta = (ex*x[i]+ey*y[i])/6371, zeta = (-ey*x[i]+ex*y[i])/6371;
      const double l = (-eta+zeta*zeta*cot/2+eta*zeta*zeta*(1./6+cot*cot/2))*6371;
      const double tau = l*(point.p+point.dp_Δ*l/2);
      for (int b = 0; b < w; ++b) {
        weights[b][i] = zero ? 0 : ((i+b)%7 ? 1.+(i%3)*.2 : 0);
        for (int k = 102; k < 102+bins; ++k) {
          auto noise = std::complex<double>(normal(gen), normal(gen));
          spec[b][i][k] = 1e-8 * (coherent ? std::polar(1., -tau*2*M_PI*STATION::df*k)+.01*noise : noise);
        }
      }
    }
  }
};
double objective(Fixture &f, bool reference, int bias = 0) {
  return (reference ? cpu_reference::cal_S : cal_S)(f.point, f.spec, f.weights,
                                                   f.x, f.y, f.windows, bias);
}
double hessian(Fixture &f, bool reference, dvector &d, dmatrix &dd) {
  return (reference ? cpu_reference::cal_HessianS : cal_HessianS)(
      f.point, f.spec, f.weights, f.x, f.y, d, dd, f.windows);
}
void compare(Fixture &f) {
  for (int bias : {0, 1}) check(objective(f, false, bias), objective(f, true, bias));
  dvector d(4), dr(4); dmatrix dd(4,4), ddr(4,4);
  check(hessian(f, false, d, dd), hessian(f, true, dr, ddr));
  for (int m = 0; m < 4; ++m) {
    check(d[m], dr[m]);
    for (int n = 0; n < 4; ++n) check(dd(m,n), ddr(m,n));
  }
}
struct Rotation {
  std::vector<STATION> sta;
  array4c spec;
  array3d weights;
  PARAM point;
  Rotation(const Fixture &f) : sta(f.stations),
      spec(boost::extents[5][f.windows][f.stations][range4c(102, STATION::nfreq)]),
      weights(boost::extents[5][f.windows][f.stations]), point(f.point) {
    STATION::lat_ary = 36; STATION::lon_ary = 138;
    for (int i=0;i<f.stations;++i) sta[i].set_station("Hi-net", "TEST", 34.+i*.001, 137.+i*.001, 0.);
    for (int c=0;c<5;++c) for(int w=0;w<f.windows;++w) for(int i=0;i<f.stations;++i) {
      weights[c][w][i] = f.weights[w][i]*(c+1);
      for(int k=102;k<STATION::nfreq;++k) spec[c][w][i][k] = f.spec[w][i][k]*double(c+1);
    }
  }
  void run(bool reference) {
    (reference ? cpu_reference::rotate_EN_RT : rotate_EN_RT)(sta, spec, weights, point);
  }
};
void rotation_check(Fixture &f) {
  Rotation a(f), b(f); a.run(false); b.run(true);
  for(size_t i=0;i<a.spec.num_elements();++i) check_complex(a.spec.data()[i],b.spec.data()[i]);
  for(size_t i=0;i<a.weights.num_elements();++i) check(a.weights.data()[i],b.weights.data()[i]);
}
volatile double sink = 0;
int benchmark(const std::string &stage, int threads, int repeats, int windows) {
  omp_set_num_threads(threads);
  Fixture f(windows,650,155,false); Rotation rotation(f);
  dvector d(4); dmatrix dd(4,4);
  auto run = [&](bool reference) {
    if(stage=="hessian") sink=hessian(f,reference,d,dd);
    else if(stage=="objective") sink=objective(f,reference);
    else if(stage=="rotation") rotation.run(reference);
    else if(stage=="grid") {
      double values[40];
#pragma omp parallel for
      for(int i=0;i<40;++i) {
        PARAM p=f.point; p.Δ=.1+i*.06;
        values[i]=(reference ? cpu_reference::cal_S : cal_S)(p,f.spec,f.weights,f.x,f.y,f.windows,0);
      }
      sink=values[0];
    } else throw std::runtime_error("unknown benchmark stage");
  };
  compare(f); rotation_check(f); run(true); run(false);
  for(int r=0;r<repeats;++r) {
    double times[2];
    for(int order=0;order<2;++order) {
      const bool reference=(order+r)%2;
      const auto start=std::chrono::steady_clock::now();run(reference);
      times[reference]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    std::cout << std::setprecision(17) << "{\"stage\":\"" << stage << "\",\"threads\":" << threads
              << ",\"windows\":" << windows << ",\"repeat\":" << r
              << ",\"reference_s\":" << times[1] << ",\"current_s\":" << times[0] << "}\n";
  }
  return 0;
}
}
int main(int argc,char **argv) try {
  omp_set_dynamic(0); STATION::horizontal_only=true;
  if(argc==6 && std::string(argv[1])=="--benchmark") {
    int threads=std::stoi(argv[3]),repeats=std::stoi(argv[4]),windows=std::stoi(argv[5]);
    if(threads<1||repeats<1||windows<1) throw std::runtime_error("positive benchmark dimensions required");
    return benchmark(argv[2],threads,repeats,windows);
  }
  if(argc!=1) throw std::runtime_error("usage: test_cpu_parallel [--benchmark STAGE THREADS REPEATS WINDOWS]");
  for(int w : {0,1,3,16,17,48}) for(int bins : {1,155}) for(bool coherent : {false,true}) {
    Fixture f(w,31,bins,coherent);
    for(int threads : {1,4,8,16}) { omp_set_num_threads(threads);compare(f); }
    for(bool horizontal : {false,true}) { STATION::horizontal_only=horizontal;rotation_check(f); }
  }
  for(int n : {1,2,650}) for(bool zero : {false,true}) {
    Fixture f(17,n,155,false,zero);
    for(int threads : {1,4,8,16}) { omp_set_num_threads(threads);compare(f); }
  }
  // Concurrent outer callers must use local scratch without nested teams.
  Fixture nested(17,31,31,true);
  double actual[8], expected[8];
#pragma omp parallel for num_threads(4)
  for(int i=0;i<8;++i) {
    dvector d(4); dmatrix dd(4,4);
    actual[i]=hessian(nested,false,d,dd)+objective(nested,false);
    expected[i]=hessian(nested,true,d,dd)+objective(nested,true);
  }
  for(int i=0;i<8;++i) check(actual[i],expected[i]);
  std::cout << std::setprecision(17) << "PASS CPU parallel: checked=" << checked << " unequal=" << unequal
            << " max_abs=" << max_abs << " max_relative=" << max_scaled << '\n';
  return 0;
} catch(const std::exception &e) { std::cerr << e.what() << '\n';return 1; }
