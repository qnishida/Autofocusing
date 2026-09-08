// Legacy loops copied from 393a93e; keep independent of the optimized kernel.
static void slant_stack_reference(const array3c &buf_aryENU, const dvector &dx_ary,
                                 const dvector &dy_ary, int num_ary, array3d &ssRTU) {
  array4c sENU(boost::extents[3][range4c(-ipmax, ipmax + 1)][range4c(
      -ipmax, ipmax + 1)][range4c(STATION::if1, STATION::if2 + 1)]);
  array2d er_x(
      boost::extents[range2d(-ipmax, ipmax + 1)][range2d(-ipmax, ipmax + 1)]);
  array2d er_y(
      boost::extents[range2d(-ipmax, ipmax + 1)][range2d(-ipmax, ipmax + 1)]);
  for (int i = -ipmax; i <= ipmax; i++)
    for (int j = -ipmax; j <= ipmax; j++)
      if (i != 0 || j != 0) {
        double px = px0 + i * dp;
        double py = py0 + j * dp;
        er_x[i][j] = px / sqrt(px * px + py * py);
        er_y[i][j] = py / sqrt(px * px + py * py);
      } else { // To avoid singularity
        er_x[i][j] = sqrt(2.) / 2.;
        er_y[i][j] = sqrt(2.) / 2.;
      }


  std::fill_n(sENU.data(), sENU.num_elements(), std::complex<double>(0., 0.));
  dvector msE(STATION::if2+1,0.), msN(STATION::if2+1,0.),
      msZ(STATION::if2+1,0.), msEN(STATION::if2+1,0.);
  for (int ist=0;ist<num_ary;++ist)
    for (int k=STATION::if1;k<=STATION::if2;++k) {
      msE[k] += norm(buf_aryENU[0][ist][k]);
      msN[k] += norm(buf_aryENU[1][ist][k]);
      msZ[k] += norm(buf_aryENU[2][ist][k]);
      msEN[k] += real(conj(buf_aryENU[0][ist][k])*buf_aryENU[1][ist][k]);
    }
  const int count_ss=num_ary;
// Using dynamic scheduling to balance the load among threads as the iterations
// may have varying execution times. The collapse(2) clause is used to collapse
// the nested loops into a single loop for better load balancing. The
// schedule(dynamic) clause is used to dynamically distribute iterations to
// threads to handle varying execution times.
#pragma omp parallel for collapse(2)                                           \
    schedule(dynamic) ////private(specE,specN,specZ)
        for (int ipy = -ipmax; ipy <= ipmax; ipy++) {
          for (int ipx = -ipmax; ipx <= ipmax; ipx++) {
            for (int ist = 0; ist < num_ary; ++ist) {
              double px = px0 + ipx * dp;
              double py = py0 + ipy * dp;

              for (int icmp = 0; icmp < 3; ++icmp) {
                double tau = -(px * dx_ary[ist] + py * dy_ary[ist]);
                double phase = tau * 2. * M_PI * STATION::df;
                double dc = cos(phase);
                double ds = sin(phase);
                double cp0 = cos(phase * STATION::if1);
                double sp0 = sin(phase * STATION::if1);
                for (int k = STATION::if1; k <= STATION::if2; ++k) {
                  double cp1 = cp0 * dc - sp0 * ds;
                  double sp1 = sp0 * dc + cp0 * ds;
                  sENU[icmp][ipx][ipy][k] += (buf_aryENU[icmp][ist][k] *
                                              std::complex<double>(cp0, sp0));
                  cp0 = cp1;
                  sp0 = sp1;
                }
              }
            }
          }
        }

        double fctr = 1. / (count_ss * (count_ss - 1)) /
                      (STATION::if2 - STATION::if1 + 1);
        for (int ipx = -ipmax; ipx <= ipmax; ipx++) {
          for (int ipy = -ipmax; ipy <= ipmax; ipy++) {
            for (int k = STATION::if1; k <= STATION::if2; ++k) {
              ssRTU[0][ipx][ipy] +=
                  ((norm((er_x[ipx][ipy] * sENU[0][ipx][ipy][k] +
                          er_y[ipx][ipy] * sENU[1][ipx][ipy][k])) -
                    er_x[ipx][ipy] * er_x[ipx][ipy] * msE[k] -
                    2 * er_x[ipx][ipy] * er_y[ipx][ipy] * msEN[k] -
                    er_y[ipx][ipy] * er_y[ipx][ipy] * msN[k]) *
                   fctr);
              ssRTU[1][ipx][ipy] +=
                  ((norm((er_y[ipx][ipy] * sENU[0][ipx][ipy][k] -
                          er_x[ipx][ipy] * sENU[1][ipx][ipy][k])) -
                    er_y[ipx][ipy] * er_y[ipx][ipy] * msE[k] +
                    2 * er_x[ipx][ipy] * er_y[ipx][ipy] * msEN[k] -
                    er_x[ipx][ipy] * er_x[ipx][ipy] * msN[k]) *
                   fctr);
              ssRTU[2][ipx][ipy] +=
                  ((norm(sENU[2][ipx][ipy][k]) - msZ[k]) * fctr);
            }
          }
        }
}
