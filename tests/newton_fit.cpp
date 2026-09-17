#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main

// Recover a known curved wavefront from a perturbed starting point. This tests
// the complete production objective/Hessian/fit, not a copied stopping rule.
int main() try {
  omp_set_dynamic(0);
  omp_set_num_threads(2);
  STATION::horizontal_only = false;
  STATION::df = 1. / 1024;
  STATION::if1 = 102;
  STATION::if2 = 132;
  STATION::nfreq = 133;
  const int stations = 49, windows = 3;
  array3c spec(boost::extents[windows][stations][range3c(102,133)]);
  array2d weights(boost::extents[windows][stations]);
  dvector x(stations), y(stations);
  PARAM truth{};
  truth.p = .06; truth.θ = .73; truth.Δ = .85; truth.dp_Δ = -2e-6;
  const double ex = cos(truth.θ), ey = sin(truth.θ);
  const double cot = cos(truth.Δ) / sin(truth.Δ);
  for (int i = 0; i < stations; ++i) {
    x[i] = (i % 7 - 3) * 370. + (i / 7) * 13.;
    y[i] = (i / 7 - 3) * 320. + (i % 7) * 17.;
    const double eta = (ex*x[i] + ey*y[i])/6371;
    const double zeta = (-ey*x[i] + ex*y[i])/6371;
    const double distance = 6371 * (-eta + zeta*zeta*cot/2
                                      + eta*zeta*zeta*(1./6 + cot*cot/2));
    const double tau = distance * (truth.p + truth.dp_Δ*distance/2);
    for (int w = 0; w < windows; ++w) {
      weights[w][i] = (i+w)%7 ? 1. + .2*(i%3) : 0;
      for (int k = 102; k < 133; ++k)
        spec[w][i][k] = std::polar(1e-8 * (1 + .1*w), -2*M_PI*STATION::df*k*tau);
    }
  }
  for (double direction : {-1., 1.}) {
    PARAM fit = truth;
    fit.p += direction * 2e-4;
    fit.θ += direction * .002;
    fit.Δ += direction * .01;
    fit.dp_Δ += direction * 2e-8;
    const double before = cal_S(fit,spec,weights,x,y,windows,0);
    const int result = est_dist_grad(fit,spec,weights,x,y,windows);
    const double after = cal_S(fit,spec,weights,x,y,windows,0);
    const double error = std::max({std::abs(fit.p-truth.p)/.06,
        std::abs(fit.θ-truth.θ)/(M_PI/2), std::abs(fit.Δ-truth.Δ)/(M_PI/2),
        std::abs(fit.dp_Δ-truth.dp_Δ)/(.04/(30.*111))});
    std::cout << std::setprecision(17) << "direction=" << direction
              << " iterations=" << result << " error=" << error
              << " improvement=" << after/before-1 << '\n';
    dvector gradient(4); dmatrix hessian(4,4);
    cal_HessianS(fit,spec,weights,x,y,gradient,hessian,windows);
    Eigen::Vector4d scales(.06,M_PI/2,M_PI/2,.04/(30.*111));
    Eigen::Matrix4d scaled_hessian;
    for (int i=0;i<4;++i) for (int j=0;j<4;++j)
      scaled_hessian(i,j)=hessian(i,j)*scales(i)*scales(j);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> curvature(scaled_hessian);
    if (result<=1 || !std::isfinite(error) || error>1e-5 || !(after>before) ||
        curvature.info()!=Eigen::Success || !(curvature.eigenvalues().array()<0).all())
      throw std::runtime_error("Known wavefront was not recovered");
  }
  PARAM invalid=truth;
  invalid.Δ=0;
  if (est_dist_grad(invalid,spec,weights,x,y,windows)!=-1)
    throw std::runtime_error("Invalid initial distance accepted");
  std::fill(spec.data(),spec.data()+spec.num_elements(),std::complex<double>(0,0));
  PARAM flat=truth;
  if (est_dist_grad(flat,spec,weights,x,y,windows)!=-1)
    throw std::runtime_error("Zero power accepted");
  std::cout << "PASS Newton wavefront recovery and rejection checks\n";
} catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
