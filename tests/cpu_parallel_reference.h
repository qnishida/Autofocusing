// Frozen numerical kernels from 977736f. Keep independent of optimized helpers.
// Includes the original rotation parallelism; objective and Hessian are serial.
namespace cpu_reference {
static double cal_S(const PARAM prm, const array3c &buf_spec,
                    const array2d &w_spec, const dvector &dx, const dvector &dy,
                    const int num_ss, const int flag_red) {
  const int num_sta = buf_spec.shape()[1];
  cvector phi(STATION::if2 + 1);

  double p = prm.p;
  double ex = cos(prm.θ);
  double ey = sin(prm.θ);

  dvector tau(num_sta, 0);
  for (int i = 0; i < num_sta; i++) {
    double η = (ex * dx[i] + ey * dy[i]) / 6371.;
    double ζ = (-ey * dx[i] + ex * dy[i]) / 6371.;
    double cotΔ = cos(prm.Δ) / sin(prm.Δ);
    double dp_Δ = prm.dp_Δ;
    double l =
        (-η + ζ * ζ * cotΔ / 2. + η * ζ * ζ * (1. / 6. + cotΔ * cotΔ / 2.)) *
        6371.;

    tau[i] = l * (p + dp_Δ * l / 2);
  }

  double S = 0;
  for (int ibuf = 0; ibuf < num_ss; ibuf++) {
    for (int k = STATION::if1; k <= STATION::if2; ++k)
      phi[k] = 0;
    for (int i = 0; i < num_sta; i++) {
      double phase = tau[i] * 2. * M_PI * STATION::df;
      double dc = cos(phase);
      double ds = sin(phase);
      double cp0 = cos(phase * STATION::if1);
      double sp0 = sin(phase * STATION::if1);

      for (int k = STATION::if1; k <= STATION::if2; ++k) {
        double cp1 = cp0 * dc - sp0 * ds;
        double sp1 = sp0 * dc + cp0 * ds;
        phi[k] += (buf_spec[ibuf][i][k] * std::complex<double>(cp0, sp0) *
                   w_spec[ibuf][i]); //*(double)flag_ss[ist]);
        cp0 = cp1;
        sp0 = sp1;
      }
    }
    for (int k = STATION::if1; k <= STATION::if2; ++k)
      S += real(conj(phi[k]) * phi[k]);
  }

  double bias = 0;
  if (flag_red == 1) {
    for (int ibuf = 0; ibuf < num_ss; ibuf++) {
      for (int i = 0; i < num_sta; i++) {
        for (int k = STATION::if1; k <= STATION::if2; ++k) {
          bias += real(conj(buf_spec[ibuf][i][k]) * buf_spec[ibuf][i][k] *
                       w_spec[ibuf][i] * w_spec[ibuf][i]);
        }
      }
    }
    // std::cout<<"###"<< S <<" "<<bias<<std::endl;
    S -= bias;
  }
  return (S);
}

static double cal_HessianS(const PARAM prm, const array3c &buf_spec,
                           const array2d &w_spec, const dvector &dx,
                           const dvector &dy, dvector &dS, dmatrix &ddS,
                           const int num_ss) {
  const int num_sta = buf_spec.shape()[1];

  dvector tau(num_sta);
  array2d dtau(boost::extents[num_sta][4]);
  array3d ddtau(boost::extents[num_sta][4][4]);

  cvector phi(STATION::if2 + 1);
  cmatrix dphi(4, STATION::if2 + 1);
  array3c ddphi(boost::extents[4][4][STATION::if2 + 1]);

  const double p = prm.p;
  const double ex = cos(prm.θ);
  const double ey = sin(prm.θ);
  const double dp_Δ = prm.dp_Δ;
  const double cotΔ = cos(prm.Δ) / sin(prm.Δ);
  const double sinΔ = sin(prm.Δ);

  for (int i = 0; i < num_sta; i++) {
    double η = (ex * dx[i] + ey * dy[i]) / 6371.;
    double ζ = (-ey * dx[i] + ex * dy[i]) / 6371.;
    double l =
        (-η + ζ * ζ * cotΔ / 2. + η * ζ * ζ * (1. / 6. + cotΔ * cotΔ / 2.)) *
        6371.;

    // double dl_phi = -ζ*(1+η*cotΔ)*6371.;
    // double dl_Δ = -pow(ζ/sinΔ,2)/2*6371.;
    // double ddl_phi2 = (η+(pow(η,2)-pow(ζ,2))*cotΔ)*6371.;
    // double ddl_phiΔ = ζ*η/(pow(sinΔ,2))*6371.;
    // double ddl_Δ2 = (pow(ζ/sinΔ,2)*cotΔ)*6371.;
    double dl_phi = (-ζ * (1 + η * cotΔ) +
                     (1 / 6. + cotΔ * cotΔ / 2.) * ζ * (ζ * ζ - 2 * η * η)) *
                    6371.;
    double dl_Δ =
        (-pow(ζ / sinΔ, 2) / 2 - η * ζ * ζ * cotΔ / (sinΔ * sinΔ)) * 6371.;
    double ddl_phi2 =
        ((η + (pow(η, 2) - pow(ζ, 2)) * cotΔ) +
         1 / 6. * (1 + 3 * cotΔ * cotΔ) * (2 * η * η * η - 7 * η * ζ * ζ)) *
        6371.;
    double ddl_phiΔ = (ζ * η / (pow(sinΔ, 2)) +
                       ζ / (sinΔ * sinΔ) * cotΔ * (2 * η * η - ζ * ζ)) *
                      6371.;
    double ddl_Δ2 =
        (pow(ζ / sinΔ, 2) * (cotΔ + η * (3 / (sinΔ * sinΔ) - 2))) * 6371.;

    tau[i] = l * (p + dp_Δ * l / 2);

    // 0: p, 1: θ, 2: Δ
    dtau[i][0] = l;
    dtau[i][1] = p * dl_phi + dp_Δ * l * dl_phi;
    dtau[i][2] = p * dl_Δ + dp_Δ * l * dl_Δ;
    dtau[i][3] = pow(l, 2) / 2.;

    ddtau[i][0][0] = 0;
    ddtau[i][0][1] = dl_phi;
    ddtau[i][0][2] = dl_Δ;
    ddtau[i][0][3] = 0.;
    ddtau[i][1][0] = ddtau[i][0][1];
    ddtau[i][1][1] = p * ddl_phi2 + dp_Δ * (pow(dl_phi, 2) + l * ddl_phi2);
    ddtau[i][1][2] = p * ddl_phiΔ + dp_Δ * (dl_phi * dl_Δ + l * ddl_phiΔ);
    ddtau[i][1][3] = l * dl_phi;
    ddtau[i][2][0] = ddtau[i][0][2];
    ddtau[i][2][1] = ddtau[i][1][2];
    ddtau[i][2][2] = p * ddl_Δ2 + dp_Δ * (pow(dl_Δ, 2) + l * ddl_Δ2);
    ddtau[i][2][3] = l * dl_Δ;
    ddtau[i][3][0] = 0.;
    ddtau[i][3][1] = ddtau[i][1][3];
    ddtau[i][3][2] = ddtau[i][2][3];
    ddtau[i][3][3] = 0.;
  }
  double S = 0;
  dS.clear();
  ddS.clear();

  for (int ibuf = 0; ibuf < num_ss; ibuf++) {
    phi.clear();
    dphi.clear();
    fill_n(ddphi.data(), ddphi.num_elements(), 0.);

    for (int i = 0; i < num_sta; i++) {
      double phase = tau[i] * 2. * M_PI * STATION::df;
      double dc = cos(phase);
      double ds = sin(phase);
      double cp0 = cos(phase * STATION::if1);
      double sp0 = sin(phase * STATION::if1);
      /// from python
      for (int k = STATION::if1; k <= STATION::if2; ++k) {
        double omega = k * STATION::df * 2. * M_PI;
        double cp1 = cp0 * dc - sp0 * ds;
        double sp1 = sp0 * dc + cp0 * ds;
        std::complex<double> amp = buf_spec[ibuf][i][k] *
                                   std::complex<double>(cp0, sp0) *
                                   w_spec[ibuf][i];
        phi[k] += amp;
        for (int m = 0; m < 4; m++) {
          dphi(m, k) += (amp * dtau[i][m] * std::complex<double>(0, 1) * omega);
          for (int n = m; n < 4; n++) {
            ddphi[m][n][k] +=
                (-pow(omega, 2) * dtau[i][m] * dtau[i][n] +
                 std::complex<double>(0, 1) * omega * ddtau[i][m][n]) *
                amp;
          }
        }
        cp0 = cp1;
        sp0 = sp1;
      }
    }
    for (int k = STATION::if1; k <= STATION::if2; ++k)
      S += real(conj(phi[k]) * phi[k]);

    for (int k = STATION::if1; k <= STATION::if2; ++k) {
      for (int m = 0; m < 4; m++) {
        dS[m] += 2 * real(conj(phi[k]) * dphi(m, k));
        for (int n = m; n < 4; n++) {
          ddS(m, n) += 2 * real(conj(dphi(m, k)) * dphi(n, k) +
                                conj(phi[k]) * ddphi[m][n][k]);
        }
      }
    }
  }
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < m; n++)
      ddS(m, n) = ddS(n, m);

  return (S);
}

static int rotate_EN_RT(const std::vector<STATION> &sta0,
                        array4c &buf_specENURT, array3d &w_specENURT,
                        const PARAM prm) {
  std::vector<double> cosbaz2(sta0.size()), sinbaz2(sta0.size());
  double evlat, evlon;
  geod.ArcDirect(STATION::lat_ary, STATION::lon_ary, 90 - prm.θ / M_PI * 180.,
                 prm.Δ / M_PI * 180., evlat, evlon);

// Parallelize the loop to compute sinbaz2 and cosbaz2 for each station
#pragma omp parallel for
  for (int ist = 0; ist < (int)sta0.size(); ist++) {
    double s12, az12, az21;
    geod.Inverse(sta0[ist].print_lat(), sta0[ist].print_lon(), evlat, evlon,
                 s12, az12, az21);
    sinbaz2[ist] = sin(az12 / 180. * M_PI);
    cosbaz2[ist] = cos(az12 / 180. * M_PI);
  }
  const int num_buffers = static_cast<int>((buf_specENURT.shape())[1]);
  const int num_stations = static_cast<int>((buf_specENURT.shape())[2]);
  for (int ibuf = 0; ibuf < num_buffers; ++ibuf) {
    for (int ist = 0; ist < num_stations; ++ist) {
      for (int k = STATION::if1; k <= STATION::if2; ++k) { // clang-format off
        buf_specENURT[3][ibuf][ist][k] =
            sinbaz2[ist] * buf_specENURT[0][ibuf][ist][k] + cosbaz2[ist] * buf_specENURT[1][ibuf][ist][k];
        buf_specENURT[4][ibuf][ist][k] =
            cosbaz2[ist] * buf_specENURT[0][ibuf][ist][k] - sinbaz2[ist] * buf_specENURT[1][ibuf][ist][k];
      }
      w_specENURT[3][ibuf][ist] = (w_specENURT[0][ibuf][ist] + w_specENURT[1][ibuf][ist])/2.;
      w_specENURT[4][ibuf][ist] = (w_specENURT[0][ibuf][ist] + w_specENURT[1][ibuf][ist])/2.;
    } // clang-format on
  }
  return 0;
}
}
