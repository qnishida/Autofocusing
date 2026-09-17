// Exercise the production conversion, including cross terms and output scaling.
#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main

static void check(const char *name, const Eigen::Matrix4d &h, double sigma) {
  dmatrix hessian(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) hessian(i,j)=h(i,j);
  PARAM p{};
  set_parameter_covariance(p,hessian,sigma);

  // Independent, equilibrated LU solve, rather than the production cofactor
  // inverse. long double may equal double on this platform. Undo equilibration
  // before output scaling.
  using Matrix = Eigen::Matrix<long double,4,4>;
  Matrix d=Matrix::Zero();
  for(int i=0;i<4;++i) d(i,i)=1/std::sqrt(std::abs(static_cast<long double>(h(i,i))));
  const Matrix normalized=d*h.cast<long double>()*d;
  const Matrix inverse=normalized.fullPivLu().solve(Matrix::Identity());
  const Matrix physical=-static_cast<long double>(sigma)*d*inverse*d;
  const double w[]={static_cast<float>(.06),static_cast<float>(M_PI/2),
                    static_cast<float>(M_PI/2),static_cast<float>(.04/(30.*111))};
  double error=0;
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
    const double expected=static_cast<double>(physical(i,j))*w[i]*w[j];
    if(!std::isfinite(p.cov[i][j])) throw std::runtime_error("Nonfinite covariance");
    error=std::max(error,std::abs(p.cov[i][j]-expected)/std::max(std::abs(expected),1e-30));
  }
  if(error>1e-8) throw std::runtime_error("Covariance or normalization mismatch");
  std::cout<<"PASS "<<name<<" maximum relative error "<<error<<'\n';
}

int main() try {
  Eigen::Matrix4d h;
  h << -4,1,0,0, 1,-3,1,0, 0,1,-3,1, 0,0,1,-2;
  check("ordinary Hessian",h,.5);
  check("bootstrap sigma scaling",h,2.);

  // Observed H/sigma for the T seed (0.063,-0.054) s/km, window beginning
  // 2005-03-16T11:51:28. All entries fit in float; its determinant does not.
  h << -201770755.95122093,118081233.05462556,-17113381529.027145,-2543079696637.4634,
       118081233.05462556,-64266904.13880946,12049242140.781115,1744873280461.91,
       -17113381529.027145,12049242140.781115,-1279808196644.9146,-194735632162049.44,
       -2543079696637.4634,1744873280461.91,-194735632162049.44,-29489748715964160.;
  const Eigen::Matrix4f old_input=h.cast<float>();
  const Eigen::Matrix4f old_inverse=old_input.inverse();
  if(!old_input.allFinite() || old_inverse.allFinite())
    throw std::runtime_error("Fixture did not reproduce the float inverse failure");
  check("observed determinant overflow",h,1.);
  const double sigma=1.1892318963033488e-8;
  check("observed raw Hessian and bootstrap sigma",(h*sigma).eval(),sigma);
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
