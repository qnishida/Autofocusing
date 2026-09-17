// Exercise the actual double objective and the experimental FP32 equivalent.
#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main
#include <iomanip>
#include <random>

static void bootstrap_components_check() {
  STATION::horizontal_only=false;
  STATION::df=1./1024; STATION::if1=102; STATION::if2=132;
  const int n=31, samples=101;
  std::mt19937 gen(71837); std::normal_distribution<double> noise;
  double max_scaled=0;
  for(int windows:{1,3,17,48}) {
    array3c spec(boost::extents[windows][n][range3c(STATION::if1,STATION::if2+1)]);
    array2d weight(boost::extents[windows][n]); dvector x(n),y(n);
    for(int i=0;i<n;++i) { x[i]=500*noise(gen); y[i]=500*noise(gen); }
    // The production call packs one selected R, T or U spectrum per event.
    for(int component=0;component<3;++component) {
      PARAM prm{}; prm.p=.03+.02*component; prm.θ=.7; prm.Δ=.9; prm.dp_Δ=-1e-6;
      for(int w=0;w<windows;++w) for(int i=0;i<n;++i)
        for(int k=STATION::if1;k<=STATION::if2;++k)
          spec[w][i][k]=std::complex<double>(noise(gen),noise(gen))*double(component+1)*1e-8;
      std::vector<double> packed,reference,scales;
      for(int b=0;b<samples;++b) {
        for(int w=0;w<windows;++w) for(int i=0;i<n;++i)
          weight[w][i]=(i%7==0)?0:((b==100)?1:((i+3*w+b)%5));
        append_weights(packed,weight,windows,n);
        reference.push_back(cal_S(prm,spec,weight,x,y,windows,1));
        scales.push_back(cal_S(prm,spec,weight,x,y,windows,0));
      }
      auto values=gpu_power_batch(power_data(spec,x,y,windows),
          std::vector<PowerPoint>(samples,power_point(prm)),packed,false,true,true);
      if(values.size()!=samples) throw std::runtime_error("Bootstrap sample count");
      for(int b=0;b<samples;++b) {
        double error=std::abs(values[b]-reference[b]);
        if(!std::isfinite(values[b]) || error>1e-12*scales[b]+1e-30)
          throw std::runtime_error("3c Bootstrap 101-sample batch precision/indexing");
        max_scaled=std::max(max_scaled,error/(scales[b]+1e-30));
      }
    }
  }
  std::cout<<"3c_bootstrap_values="<<4*3*samples<<" max_scaled="<<max_scaled<<'\n';
  STATION::horizontal_only=true;
}

static int grid_benchmark() {
  STATION::horizontal_only=true; // GPU initial-grid evaluation is horizontal-only.
  STATION::df=1./1024;STATION::if1=102;STATION::if2=256;
  const int n=650,windows=3;
  array3c spec(boost::extents[windows][n][range3c(STATION::if1,STATION::if2+1)]);
  array2d weights(boost::extents[windows][n]); dvector x(n),y(n);
  std::mt19937 gen(98231);std::uniform_real_distribution<double> coord(-1500,1500);
  for(int i=0;i<n;++i) {x[i]=coord(gen);y[i]=coord(gen);}
  PARAM truth{};truth.p=.06;truth.θ=.73;truth.Δ=M_PI/4;truth.dp_Δ=0;
  for(int w=0;w<windows;++w) for(int i=0;i<n;++i) {
    weights[w][i]=(i%7)?1:0;
    double eta=(cos(truth.θ)*x[i]+sin(truth.θ)*y[i])/6371;
    double zeta=(-sin(truth.θ)*x[i]+cos(truth.θ)*y[i])/6371;
    double l=(-eta+zeta*zeta/2+eta*zeta*zeta*(1./6+.5))*6371;
    for(int k=STATION::if1;k<=STATION::if2;++k) spec[w][i][k]=std::polar(1e-8,-l*truth.p*2*M_PI*STATION::df*k);
  }
  // Warm shader compilation before timing; all packing/sync remains measured.
  std::vector<PARAM> warm(1,truth);dvector warm_values(1);
  grid_powers(warm_values,warm,1,spec,weights,x,y,windows);
  for(bool wider:{false,true}) for(int factor:{1,10,100}) for(int repeat=0;repeat<5;++repeat) {
    PARAM results[2];double times[2];
    for(int order=0;order<2;++order) {
      bool gpu=(order+(repeat%2))%2;PARAM selected=truth;selected.Δ=-1;selected.dp_Δ=0;
      auto start=std::chrono::steady_clock::now();
      for(int stage=0;stage<2;++stage) {
        int count=(stage?40:35)*factor;std::vector<PARAM> points(count,selected);dvector values(count);
        for(int i=0;i<count;++i) {
          if(stage==0) points[i].Δ=((wider?1.:5.)+(wider?178.:170.)*i/(count-1))*M_PI/180;
          else points[i].dp_Δ=-2e-5+(wider?6e-5:3e-5)*i/count;
        }
        if(gpu) grid_powers(values,points,count,spec,weights,x,y,windows);
        else {
#pragma omp parallel for
          for(int i=0;i<count;++i) values[i]=cal_S(points[i],spec,weights,x,y,windows,0);
        }
        double max=0;for(int i=0;i<count;++i) if(max<values[i]) {max=values[i];selected=points[i];}
      }
      times[gpu]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();results[gpu]=selected;
    }
    if(results[0].Δ!=results[1].Δ || results[0].dp_Δ!=results[1].dp_Δ) throw std::runtime_error("Expanded grid winner mismatch");
    std::cout<<std::setprecision(17)<<"{\"factor\":"<<factor<<",\"wider\":"<<(wider?"true":"false")<<",\"repeat\":"<<repeat
      <<",\"cpu_s\":"<<times[0]<<",\""<<AUTOFOCUSING_TEST_BACKEND<<"_s\":"<<times[1]<<",\"same_winner\":true}"<<std::endl;
  }
  return 0;
}

int main(int argc,char **argv) try {
  omp_set_dynamic(0);
  if(!std::getenv("OMP_NUM_THREADS")) omp_set_num_threads(16);
  setenv("AUTOFOCUSING_BACKEND",AUTOFOCUSING_TEST_BACKEND,1);
  if(argc==2 && std::string(argv[1])=="--benchmark-grid") return grid_benchmark();
  setenv("AUTOFOCUSING_GPU_POWER","all",1);
  setenv("AUTOFOCUSING_BACKEND","cpu",1); STATION::horizontal_only=true;
  if(power_enabled("bootstrap")) throw std::runtime_error("CPU must retain CPU objectives");
  setenv("AUTOFOCUSING_BACKEND",AUTOFOCUSING_TEST_BACKEND,1); STATION::horizontal_only=false;
  const bool cuda=std::string(AUTOFOCUSING_TEST_BACKEND)=="cuda";
  if(power_enabled("grid") || power_enabled("bootstrap")!=cuda)
    throw std::runtime_error("3c requires CUDA FP64 Bootstrap and CPU initial grids");
  STATION::horizontal_only=true;
  if(!power_enabled("bootstrap") || !power_enabled("grid")) throw std::runtime_error("Opt-in dispatch");
  setenv("AUTOFOCUSING_GPU_POWER","off",1);
  if(power_enabled("grid")) throw std::runtime_error("Off dispatch");
  if(cuda) bootstrap_components_check();
  STATION::df=1./1024;
  STATION::if1=102; STATION::if2=256;
  std::mt19937 generator(1837);
  std::normal_distribution<double> normal;
  std::uniform_real_distribution<double> coordinate(-1500,1500);
  unsigned cases=0;
  double max_float=0, max_double=0;
  for(int stations:{2,31,650}) for(int windows:{1,3}) {
    array3c spectra(boost::extents[windows][stations][range3c(STATION::if1,STATION::if2+1)]);
    array2d weights(boost::extents[windows][stations]);
    dvector dx(stations),dy(stations);
    for(int i=0;i<stations;++i) {dx[i]=coordinate(generator);dy[i]=coordinate(generator);}
    for(double delta:{5.*M_PI/180,M_PI/2,175.*M_PI/180}) {
      PARAM prm{}; prm.p=.165;prm.θ=.73;prm.Δ=delta;prm.dp_Δ=-2e-5;
      for(bool coherent:{false,true}) for(double scale:{1e-12,1.,1e6}) {
        for(int w=0;w<windows;++w) for(int i=0;i<stations;++i) {
          weights[w][i]=(i%7==0 && stations>2)?0:1+(i%3);
          double ex=cos(prm.θ),ey=sin(prm.θ),cot=cos(delta)/sin(delta);
          double eta=(ex*dx[i]+ey*dy[i])/6371, zeta=(-ey*dx[i]+ex*dy[i])/6371;
          double l=(-eta+zeta*zeta*cot/2+eta*zeta*zeta*(1./6+cot*cot/2))*6371;
          double tau=l*(prm.p+prm.dp_Δ*l/2);
          for(int k=STATION::if1;k<=STATION::if2;++k) {
            auto noise=std::complex<double>(normal(generator),normal(generator));
            spectra[w][i][k]=scale*(coherent?std::polar(1.,-tau*2*M_PI*STATION::df*k)+.01*noise:noise);
          }
        }
        double reference_scale=cal_S(prm,spectra,weights,dx,dy,windows,0);
        for(int bias:{0,1}) for(bool fp64:{false,true}) {
          if(fp64 && !cuda) continue;
          double reference=cal_S(prm,spectra,weights,dx,dy,windows,bias);
          auto data=power_data(spectra,dx,dy,windows);
          data.double_precision=fp64;
          std::vector<double> packed; append_weights(packed,weights,windows,stations);
          double actual=gpu_power_batch(data,{power_point(prm)},packed,true,true,bias)[0];
          double a=std::abs(actual-reference)/reference_scale;
          if(!std::isfinite(a)||a>(fp64?1e-12:5e-4))
            throw std::runtime_error("Objective smoke-test tolerance exceeded");
          if(fp64) max_double=std::max(max_double,a);
          else max_float=std::max(max_float,a);
          ++cases;
        }
      }
    }
  }
  {
    const int n=31, windows=3;
    array3c spec(boost::extents[windows][n][range3c(STATION::if1,STATION::if2+1)]);
    array2d weight(boost::extents[windows][n]); dvector x(n),y(n);
    for(int i=0;i<n;++i) {x[i]=coordinate(generator); y[i]=coordinate(generator);}
    for(int w=0;w<windows;++w) for(int i=0;i<n;++i) for(int k=STATION::if1;k<=STATION::if2;++k)
      spec[w][i][k]={normal(generator),normal(generator)};
    auto data=power_data(spec,x,y,windows);
    std::vector<PowerPoint> points; std::vector<double> packed,refs,scales;
    for(int b=0;b<65;++b) {
      PARAM p{}; p.p=.06; p.θ=.1*b; p.Δ=(5+2.5*b)*M_PI/180; p.dp_Δ=1e-6;
      points.push_back(power_point(p));
      for(int w=0;w<windows;++w) for(int i=0;i<n;++i) weight[w][i]=(i+b)%5;
      append_weights(packed,weight,windows,n);
      refs.push_back(cal_S(p,spec,weight,x,y,windows,1));scales.push_back(cal_S(p,spec,weight,x,y,windows,0));
    }
    auto values=gpu_power_batch(data,points,packed,false,false,true);
    for(int b=0;b<65;++b) if(std::abs(values[b]-refs[b])>5e-4*scales[b]+1e-30)
      throw std::runtime_error("Batched objective indexing/precision");
    // Every sharing combination spans multiple 32-candidate chunks, including
    // Bootstrap's shared geometry with different weights for each replicate.
    for (bool shared_weights:{false,true}) for (bool shared_geometry:{false,true}) {
      auto shared_points=points;
      if(shared_geometry) std::fill(shared_points.begin(),shared_points.end(),points[0]);
      std::vector<double> input_weights=packed;
      if(shared_weights) input_weights.resize(windows*n);
      auto batch=gpu_power_batch(data,shared_points,input_weights,shared_weights,shared_geometry,true);
      for(int b=0;b<65;++b) {
        const auto &v=shared_points[b]; PARAM p{};p.p=v.p;p.θ=v.theta;p.Δ=v.delta;p.dp_Δ=v.curvature;
        for(int w=0;w<windows;++w) for(int i=0;i<n;++i)
          weight[w][i]=input_weights[(shared_weights?0:b*windows*n)+w*n+i];
        double reference=cal_S(p,spec,weight,x,y,windows,1);
        double scale=cal_S(p,spec,weight,x,y,windows,0);
        if(!std::isfinite(batch[b]) || std::abs(batch[b]-reference)>5e-4*scale+1e-30)
          throw std::runtime_error("Shared geometry/weight chunk boundary");
      }
    }
    // Exact ties and zero signal must preserve the original strict-max policy.
    std::fill_n(spec.data(),spec.num_elements(),std::complex<double>(0,0));
    std::vector<PARAM> tied(40);for(auto &p:tied) {p.p=.06;p.θ=.2;p.Δ=1.;p.dp_Δ=0;}
    dvector values_tied(40);grid_powers(values_tied,tied,40,spec,weight,x,y,windows);
    for(double v:values_tied) if(v!=0) throw std::runtime_error("Zero/tie objective");
  }
  if(std::string(AUTOFOCUSING_TEST_BACKEND)=="cuda") {
    PowerData invalid{};
    bool rejected=false;
    try { cuda_power_batch(invalid,{}, {},true,true,false); }
    catch(const std::invalid_argument &) { rejected=true; }
    if(!rejected) throw std::runtime_error("Invalid CUDA power dimensions accepted");
    PowerData overflow{UINT32_MAX,UINT32_MAX,UINT32_MAX,0,1.,{}, {}};
    rejected=false;
    try { cuda_power_batch(overflow,{{.1,0,1,0}}, {},true,true,false); }
    catch(const std::invalid_argument &) { rejected=true; }
    if(!rejected) throw std::runtime_error("CUDA power dimension overflow accepted");
  }
  std::cout<<std::setprecision(17)<<"cases="<<cases<<" max_float_scaled="<<max_float
           <<" max_double_scaled="<<max_double
           <<" (scaled by uncorrected double power)\n";
  return 0;
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n';return 1; }
