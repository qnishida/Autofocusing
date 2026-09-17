#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main

static void expect(const char *name, bool actual, bool wanted) {
  if(actual!=wanted) throw std::runtime_error(name);
  std::cout<<"PASS "<<name<<'\n';
}

int main() try {
  const double eps=std::numeric_limits<double>::epsilon();
  const Eigen::Vector4d eig=-Eigen::Vector4d::Ones();
  const Eigen::Vector4d tiny=Eigen::Vector4d::Constant(1e-9);
  auto test=[&](const char *name,double f,double trial,const Eigen::Vector4d &e,
                const Eigen::Vector4d &step,double gain,bool wanted) {
    expect(name,newton_roundoff_converged(f,trial,e,step,gain),wanted);
  };
  // Concave quadratic f(x)=1-|x|^2/2 at x=(1e-9,...): the exact Newton
  // improvement is positive but both objective values round to one.
  const double gain=.5*tiny.squaredNorm();
  test("rounded concave maximum",1-gain,1,eig,tiny,gain,true);
  test("exact stationary maximum",1,1,eig,Eigen::Vector4d::Zero(),0,true);
  test("tiny negative rounding error",1,1-eps,eig,tiny,gain,true);
  test("observed gain and update bounds",1.1304596832815571e-6,1.1304596832815571e-6,
       eig,Eigen::Vector4d::Constant(1.0541300957242586e-8/2),
       5.5164869323820266e-16*1.1304596832815571e-6,true);
  for(double scale:{1e-20,1e20})
    test("power-unit invariance",scale,scale,eig*scale,tiny,gain*scale,true);
  // A flat quadratic can round to the same power while its Newton step is
  // still physically meaningful. A tie alone must not imply convergence.
  test("flat region with large full step",1,1,eig*1e-20,Eigen::Vector4d::Ones(),2e-20,false);
  test("resolvable predicted improvement",1,1,eig,tiny,1e-9,false);
  test("genuine objective decrease",1,1-1e-9,eig,tiny,gain,false);
  Eigen::Vector4d saddle=eig;saddle(0)=1;
  test("saddle point",1,1,saddle,tiny,gain,false);
  saddle(0)=0;
  test("singular Hessian",1,1,saddle,tiny,gain,false);
  test("negative predicted gain",1,1,eig,tiny,-gain,false);
  test("zero objective",0,0,eig,tiny,gain,false);
  test("negative objective",-1,-1,eig,tiny,gain,false);
  for(double invalid:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    test("nonfinite objective",invalid,1,eig,tiny,gain,false);
    test("nonfinite trial",1,invalid,eig,tiny,gain,false);
    test("nonfinite prediction",1,1,eig,tiny,invalid,false);
    Eigen::Vector4d vector=tiny;vector(0)=invalid;
    test("nonfinite update",1,1,eig,vector,gain,false);
    vector=eig;vector(0)=invalid;
    test("nonfinite curvature",1,1,vector,tiny,gain,false);
  }
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
