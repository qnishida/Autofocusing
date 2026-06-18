#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <chrono>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <GeographicLib/Geodesic.hpp>
#include <boost/format.hpp>
#include <boost/random.hpp>
// #define DEBUG
#define BOOST_DISABLE_ASSERTS
using namespace GeographicLib;
using namespace boost::posix_time;
using namespace boost::gregorian;

namespace fs = boost::filesystem;
// Now the definition of vector and multiarray is in util.h

// static const int yr0 = 2004;//2015;
static double dp = 5.0E-3;              // 2.5E-3;
static int ipmax = (int)(1.65E-1 / dp); // 1.5E-1
static double px0 = 0;                  //-0.007;
static double py0 = 0;                  // 0.047;
static double rad0 = 1800.;             // The radius of the initial beamforming
static double rad1 = 1800.;             // The radius of the array is 500 km.;
static int event_number = 30;
static const Geodesic &geod = Geodesic::WGS84();

#include "calTT.h"
#include "station_info.h"
#include "util.h"
#include <math.h>
#include <omp.h>

// Prefetch
#include <fcntl.h>
void prefetchFile(const std::string &dir_Hinet, const date &d) {
  date d_next = d + date_duration(1);
  std::string h5file_next;
  std::string h5dir =
      dir_Hinet + (boost::format("%d/%02d%02d") % d_next.year() %
                   (int)(d_next.month()) % d_next.day())
                      .str();
  const fs::path path_1(h5dir);
  if (!fs::exists(path_1))
    return;
  BOOST_FOREACH (const fs::path &p,
                 std::make_pair(fs::directory_iterator(path_1),
                                fs::directory_iterator())) {
    if (!fs::is_directory(p)) {
      h5file_next = p.string();
    }
  }
  int fd = open(h5file_next.c_str(), O_RDONLY);
  if (fd == -1) {
    perror("open");
    return;
  }
  // Tell the OS to “read the entire file sequentially
  if (posix_fadvise(fd, 0, 0, POSIX_FADV_WILLNEED) != 0) {
    perror("posix_fadvise");
  }
  close(fd);
  return;
}

struct PARAM {
  int num_eig;
  double p, θ, Δ, dp_Δ, fmax; // p [s/km], θ [rad], dp_Δ [s/km/km]
  double dp, dD, dS, ε, S_boot, sigma, MS_max;
  double cov[4][4];
};

static int search_events(std::vector<STATION> &sta0, array3d &ssRTU,
                         const int iseg, const int num_segments,
                         std::ofstream &ofs);
static double cal_S(const PARAM prm, const array3c &buf_spec,
                    const array2d &w_spec, const dvector &dx, const dvector &dy,
                    const int num_ss, const int flag_red);
static double est_fmax(const PARAM prm, const array3c &buf_spec,
                       const array2d &w_spec, const dvector &dx,
                       const dvector &dy, const int num_ss);
static double cal_HessianS(const PARAM prm, const array3c &buf_spec,
                           const array2d &w_spec, const dvector &dx,
                           const dvector &dy, dvector &dS, dmatrix &ddS,
                           const int num_ss);
static cmatrix cal_S_matrix(const PARAM prm, const array4c &buf_specENURT,
                            const array3d &w_spec, const dvector &dx,
                            const dvector &dy, const int num_ss);

static int rotate_EN_RT(const std::vector<STATION> &sta0,
                        array4c &buf_specENURT, array3d &w_specENURT,
                        const PARAM prm);
static int est_dist_grid(PARAM &prm, const array3c &buf_spec,
                         const array2d &w_spec, const dvector &dx,
                         const dvector &dy, const int num_ss);
static int est_dist_grad(PARAM &prm, const array3c &buf_spec,
                         const array2d &w_spec, const dvector &dx,
                         const dvector &dy, const int num_ss);
static int est_dist_boot(PARAM &prm, const array3c &buf_spec,
                         const array2d &w_spec, const dvector &dx,
                         const dvector &dy, const int num_ss,
                         const double bias);

static int search_max(const array3d &ssRTU, dvector &px, dvector &py,
                      const int event_number, double &mad,
                      std::vector<double> &quartile, const int icmp);
static void output_result(std::ostream &ofs, const PARAM &prm, int icmp,
                          double max, double mad, int num_ss, const ptime &t_s,
                          const ptime &t_e, cmatrix &S_matrix, int flag2);
// For definition of transients
int eval_deri(std::vector<STATION> &sta0, double th, double &median0,
              double &median1);

static double tmp_integ0[2000], tmp_integ1[2000], integ_pre = 0, integ_old = 0;
static int count_gap = 0;
std::vector<ptime> CMT_ptime;
std::vector<double> CMT_slat, CMT_slon, CMT_sdep, CMT_moment;

int main(int argc, char *argv[]) {
  STATION::dt_msec = 500;  // Sampling interval in millisecond
  STATION::len = 1024 * 2; // 2^x
  STATION::stride = 928;   // 667*2; //667*2*128 = (86400*2-1024*2) //824;
  STATION::df = 1. / (STATION::len * STATION::dt_msec * 1E-3);
  STATION::nfreq = (int)(2.6E-1 / STATION::df); // STATION::len/2;
  STATION::npts = 86400 * 1000 / (STATION::dt_msec);
  STATION::init_Freq();
  STATION::nl_h = 50;
  STATION::nl_v = 50;
  // STATION::flag_ary=0;
  STATION::if1 = (int)(0.1 / STATION::df);
  STATION::if2 = (int)(0.25 / STATION::df);

  // init_ttTable(); // Initialize the travel time table based on IASP91 model

  // Cache setting
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  size_t rdcc_nslots = 1049;       // Hash entry
  size_t rdcc_nbytes = 2147483648; // 2GB cache
  double rdcc_w0 = 0.75;           // Balance of writing
  H5Pset_cache(fapl, 0, rdcc_nslots, rdcc_nbytes, rdcc_w0);

  int yr0;
  std::ofstream ofs; // スコープの外で宣言
  std::string dir_Hinet;
  std::string cmt_catalog;

  if (argc == 6) {
    yr0 = std::stoi(argv[1]);

    std::string param_id = argv[2];
    std::string git_version = argv[3];
    dir_Hinet = argv[4];
    if (!dir_Hinet.empty() && dir_Hinet[dir_Hinet.size() - 1] != '/')
      dir_Hinet += "/";
    cmt_catalog = argv[5];
    ofs.open("output/" + param_id + "/" + git_version + "/" +
             std::to_string(yr0) + "_" + std::to_string(STATION::len) + "_" +
             std::to_string(STATION::if1 * STATION::df) + "-" +
             std::to_string(STATION::if2 * STATION::df) + ".dat");
    if (!ofs.is_open()) {
      std::cerr << "Cannot open output file for param-id=" << param_id
                << " git-version=" << git_version << " year=" << yr0
                << std::endl;
      H5Pclose(fapl);
      return 1;
    }
  } else {
    std::cerr << "Usage:   ./bin/cal_ccf YYYY <param-id> <GIT_version> "
                 "<hinet-root> <cmt-catalog>"
              << std::endl;
    std::cerr << "Output:  ./output/<param-id>/<GIT_version>/YYYY.*.dat"
              << std::endl;
    H5Pclose(fapl);
    return 1;
  }

  if (init_CMT(2004, cmt_catalog) != 0) {
    H5Pclose(fapl);
    return 1;
  }

  get_CMT(CMT_ptime, CMT_slat, CMT_slon, CMT_sdep, CMT_moment);
  std::cout << "#" << STATION::if1 * STATION::df << " "
            << STATION::if2 * STATION::df << std::endl;

  // Definitions of variables for calculations of slant stack.
  array3d ssRTU(boost::extents[3][range2d(-ipmax, ipmax + 1)]
                              [range2d(-ipmax, ipmax + 1)]);

  int sta_num = 0;
  int count0 = 0;

  date d_end(2024, 12, 31);
  auto start_tm = std::chrono::system_clock::now();
  for (int days = 0; days < 366 * 20.75;
       days += 1) { //*14.67 //  for(int days=0;days<365*16.75;days+=1){//*14.67
                    // for(int days=0;days<366;days+=1){//*14.67 //  for(int
    // days=0;days<365*16.75;days+=1){//*14.67
    date d(date(yr0, 1, 1) +
           date_duration(
               days)); // 4/ /date d(date(yr0,4,1)+date_duration(days));//4
    // for(int days=0;days<5;days+=1){//*14.67 //  for(int
    // days=0;days<365*16.75;days+=1){//*14.67 date
    // d(date(2014,12,8)+date_duration(days));//4/ /date
    // d(date(yr0,4,1)+date_duration(days));//4 date
    // d(date(2004,9,9)+date_duration(days));//4//Test std::string h5file =
    // dir_Hinet+ (boost::format("%d/%02d%02d/%d%03d0000.h5")%d.year()%
    // (int)(d.month())%d.day()%d.year() %((d-d0).days()+1)).str(); if(d.year()
    // != yr0)continue;//One year stack
    if (d > d_end)
      break;

    std::string h5file;
    int count_files = 0;
    std::string h5dir = dir_Hinet + (boost::format("%d/%02d%02d") % d.year() %
                                     (int)(d.month()) % d.day())
                                        .str();
    const fs::path path_1(h5dir);
    if (!fs::exists(path_1))
      continue;
    BOOST_FOREACH (const fs::path &p,
                   std::make_pair(fs::directory_iterator(path_1),
                                  fs::directory_iterator())) {
      if (!fs::is_directory(p)) {
        h5file = p.string();
        count_files++;
      }
    }
    prefetchFile(dir_Hinet, d);

    auto now_tm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = now_tm - start_tm;
    std::cerr << "#Start" << h5file << " " << sta_num << " "
              << elapsed_seconds.count() << "s\n";

    const fs::path path(h5file);
    if (count_files == 1 && fs::exists(path)) {
      std::vector<STATION> sta0;
      init_station(h5file, sta0, rad0, rad1, fapl);

      // Clear sac information
      for (int i = 0; i < (int)sta0.size(); i++)
        sta0[i].clear_sac();
      std::fill_n(ssRTU.data(), ssRTU.num_elements(), 0.);

      now_tm = std::chrono::system_clock::now();
      elapsed_seconds = now_tm - start_tm;
      std::cerr << "#Read data: init_station " << h5file << " " << sta_num
                << " " << elapsed_seconds.count() << "s\n";

      sta_num = read_h5(sta0, h5file, fapl);

      now_tm = std::chrono::system_clock::now();
      elapsed_seconds = now_tm - start_tm;
      std::cerr << "#Read data: read_h5 " << h5file << " " << sta_num << " "
                << elapsed_seconds.count() << "s\n";

      const int num_segments = 4;
      for (int iloc = 0; iloc < num_segments; iloc++) {
        if (sta_num > 0)
          count0 += search_events(sta0, ssRTU, iloc, num_segments, ofs);
      }
    }
  }
  H5Pclose(fapl);
  return 0;
}

/**
 * @brief Searches for events in seismic data and performs slant stack analysis.
 *
 * This function processes seismic data from multiple stations to identify
 * events and perform slant stack analysis and auto-focusing method. It
 * calculates various parameters and writes the results to an output file
 * stream. Power of the initial beamforming is normalized by each segment (1024
 * s typically), whereas power in cal_S() is normalized by the num_segments (1/4
 * day typically).
 *
 * @param sta0 A reference to a vector of STATION objects containing seismic
 * data.
 * @param ssRTU A reference to a 3D array for storing slant stack results.
 * @param iseg The segment index to process.
 * @param num_segments Number of segment divisions. 4 means dividing one day
 * into 4 segments.
 * @param ofs A reference to an output file stream for writing results.
 * @return The number of slant stack segments processed.
 */
static int search_events(std::vector<STATION> &sta0, array3d &ssRTU,
                         const int iseg, const int num_segments,
                         std::ofstream &ofs) {
  ptime t1, t2;
  SPCTRM specE, specN, specZ;
  double median0, median1;
  // Determine reference time t0
  ptime t0;
  bool flag = 0;
  for (int i = 0; i < (int)sta0.size(); i++) {
    if (sta0[i].print_sta_num() == 3) {
      if (flag == 0) {
        t0 = sta0[i].print_sacE().ts;
        flag = 1;
      }
      if (t0 > sta0[i].print_sacE().ts) {
        t0 = sta0[i].print_sacE().ts;
        flag = 1;
      }
    }
  }
  // Definition for buffer for slant stack
  array3c buf_aryENU(boost::extents[3][(int)sta0.size()]
                                   [range3c(STATION::if1, STATION::if2 + 1)]);

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

  dvector dx_ary((int)sta0.size()), dy_ary((int)sta0.size());
  std::vector<double> integH_tmp(2 * (int)sta0.size()),
      integV_tmp((int)sta0.size());

  long len = STATION::len * STATION::dt_msec;
  long stride_msec = STATION::stride * STATION::dt_msec;

  dvector msE(STATION::if2 + 1, 0.), msEN(STATION::if2 + 1, 0.),
      msN(STATION::if2 + 1, 0.), msZ(STATION::if2 + 1, 0.);

  int num_ss = 0;
  std::vector<int> flag_nl((int)sta0.size(), 0);
  std::vector<int> flag_ss((int)sta0.size(), 0);

  int count_all = 0;
  for (int ist = 0; ist < (int)sta0.size(); ist++)
    if (sta0[ist].print_rad() < rad1)
      count_all++;
  long pos_start, pos_end;
  if (STATION::stride == 928) {
    if (num_segments == 4) { // 185 = (86400*2-2048)/928 + 1
      pos_start = iseg * 46 * 928 * STATION::dt_msec;
      pos_end = pos_start + ((46 - 1) * 928 + 2048) * STATION::dt_msec;
      if (iseg == 3)
        pos_end += (928 * STATION::dt_msec);
    } else
      return (0);
  } else
    return (0);

  const int len_buf = (pos_end - pos_start) / stride_msec + 1;
  array4c buf_specENURT(boost::extents[5][len_buf][(int)sta0.size()]
                                      [range4c(STATION::if1, STATION::nfreq)]);
  array3d w_specENURT(boost::extents[5][len_buf][(int)sta0.size()]);

  for (long int_t = pos_start; int_t < pos_end; int_t += stride_msec) {
    t1 = t0 + milliseconds(int_t);
    t2 = t1 + milliseconds(len);
    msZ.clear();
    msE.clear();
    msEN.clear();
    msN.clear();
    fill_n(sENU.data(), sENU.num_elements(), 0.);

    if (CMT(t1, t2) == 0) {
#pragma omp parallel for private(specE, specN, specZ)
      for (int ist = 0; ist < (int)sta0.size(); ist++) {
        flag_nl[ist] = 0;
        if (sta0[ist].print_sta_num() && sta0[ist].cal_spec(t1)) {
          sta0[ist].print_spec(specE, specN, specZ);
          if (1.0 < specZ.integ[0] && specZ.integ[0] < 4E6    // 4E6
              && 1.0 < specZ.integ[1] && specZ.integ[1] < 4E5 // 4E5
              && 1.0 < specZ.integ[2] && specZ.integ[2] < 4E5 // 4E5
              && 1.0 < specE.integ[0] && specE.integ[0] < 4E6 &&
              1.0 < specE.integ[1] && specE.integ[1] < 4E5 &&
              1.0 < specE.integ[2] && specE.integ[2] < 4E5 &&
              1.0 < specN.integ[0] && specN.integ[0] < 4E6 &&
              1.0 < specN.integ[1] && specN.integ[1] < 4E5 &&
              1.0 < specN.integ[2] && specN.integ[2] < 4E5) {
            flag_nl[ist] = 1;
          }
        }
      }
      // Calculation of median of integ to reject the outliers
      for (int ist = 0; ist < (int)sta0.size(); ist++) {
        sta0[ist].print_spec(specE, specN, specZ);
        integH_tmp[ist * 2] = specE.integ[2];
        integH_tmp[ist * 2 + 1] = specN.integ[2];
        integV_tmp[ist] = specZ.integ[2];
      }
      double medianH =
          boost::math::statistics::median(integH_tmp); // Change the order
      double medianV = boost::math::statistics::median(integV_tmp);
      for (int ist = 0; ist < (int)sta0.size(); ist++) {
        sta0[ist].print_spec(specE, specN, specZ);
        if (specZ.integ[2] > medianV * 20 || specE.integ[2] > medianH * 20 ||
            specN.integ[2] > medianH * 20) {
          flag_nl[ist] = 0;
        }
      }

      int count_ss = 0;
      int count_est = 0;
      for (int ist = 0; ist < (int)sta0.size(); ist++) {
        if (flag_nl[ist] == 1 && sta0[ist].print_rad() < rad1) {
          flag_ss[ist] = 1;
          count_ss++;
        } else
          flag_ss[ist] = 0;
        if (flag_nl[ist] == 1)
          count_est++;
      }
      int flag_deri = eval_deri(sta0, 1.12, median0, median1);
#ifdef DEBUG
      std::cerr << to_iso_string(t1) << " " << count_ss / (double)count_all
                << " " << count_est / (double)sta0.size() << " " << flag_deri
                << " " << median0 << " " << median1 << std::endl;
#endif
      if (flag_deri && count_ss > count_all * .8 &&
          count_est > 0.8 * (int)(sta0.size())) {
        // Copy of specZ for further loop of the parameter search
        int num_ary = 0;
        for (int ist = 0; ist < (int)sta0.size(); ist++) {
          for (int icmp = 0; icmp < 3; ++icmp) {
            if (flag_nl[ist] == 1)
              w_specENURT[icmp][num_ss][ist] = 1.;
            else
              w_specENURT[icmp][num_ss][ist] = 0.;
          }

          sta0[ist].print_spec(specE, specN, specZ);
          for (int k = STATION::if1; k < STATION::nfreq; ++k) {
            buf_specENURT[0][num_ss][ist][k] = specE.spec[k];
            buf_specENURT[1][num_ss][ist][k] = specN.spec[k];
            buf_specENURT[2][num_ss][ist][k] = specZ.spec[k];
          }
          if (flag_ss[ist] == 1) {
            for (int k = STATION::if1; k <= STATION::if2; ++k)
              buf_aryENU[0][num_ary][k] = specE.spec[k];
            for (int k = STATION::if1; k <= STATION::if2; ++k)
              buf_aryENU[1][num_ary][k] = specN.spec[k];
            for (int k = STATION::if1; k <= STATION::if2; ++k)
              buf_aryENU[2][num_ary][k] = specZ.spec[k];

            dx_ary[num_ary] = sta0[ist].print_dx();
            dy_ary[num_ary] = sta0[ist].print_dy();

            num_ary++;
          }
        }
        num_ss++;
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

        for (int ist = 0; ist < (int)sta0.size(); ist++) {
          if (flag_ss[ist] == 1) {
            sta0[ist].print_spec(specE, specN, specZ);
            for (int k = STATION::if1; k <= STATION::if2; ++k) {
              msE[k] += norm(specE.spec[k]);
              msN[k] += norm(specN.spec[k]);
              msZ[k] += norm(specZ.spec[k]);
              msEN[k] += real(conj(specE.spec[k]) * specN.spec[k]);
            }
          }
        }

        std::cerr << t1 << " " << median0 << " " << median1 << " " << flag_deri
                  << " " << count_ss << std::endl;
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
    }
  }
  if (num_ss == 0)
    return (0);
  for (int icmp = 0; icmp < 3; ++icmp) {
    for (int ipx = -ipmax; ipx <= ipmax; ipx++) {
      for (int ipy = -ipmax; ipy <= ipmax; ipy++) {
        ssRTU[icmp][ipx][ipy] /= num_ss;
      }
    }
  }
  double mad;
  dvector pxV(event_number), pyV(event_number);
  dvector px(event_number * 2), py(event_number * 2);
  dvector dx((int)sta0.size()), dy((int)sta0.size());
  for (int ist = 0; ist < (int)sta0.size(); ist++) {
    dx[ist] = sta0[ist].print_dx();
    dy[ist] = sta0[ist].print_dy();
  }

  ptime t_s = t0 + milliseconds(pos_start);
  ptime t_e = t0 + milliseconds(pos_end);
  std::vector<double> quartile(3, 0.0);
  PARAM prm;

  // Search for maximum events in the slant stack results in vertical
  // components. icmp = 2
  int l2V = search_max(ssRTU, pxV, pyV, event_number, mad, quartile, 2);
  for (int i = 0; i < l2V; i++) {
    double max =
        ssRTU[2][round(pxV[i] / dp + 1E-10)][round(pyV[i] / dp + 1E-10)];

    prm.p = sqrt(pxV[i] * pxV[i] + pyV[i] * pyV[i]);
    prm.θ = atan2(pyV[i], pxV[i]);
    prm.Δ = M_PI / 2.;
    prm.dp_Δ = 0.;
    rotate_EN_RT(sta0, buf_specENURT, w_specENURT, prm);
    est_dist_grid(prm, buf_specENURT[2], w_specENURT[2], dx, dy, num_ss);
    rotate_EN_RT(sta0, buf_specENURT, w_specENURT, prm);
    int flag2 =
        est_dist_grad(prm, buf_specENURT[2], w_specENURT[2], dx, dy, num_ss);
    rotate_EN_RT(sta0, buf_specENURT, w_specENURT, prm);
    if (flag2 > 1) {
      est_dist_boot(prm, buf_specENURT[2], w_specENURT[2], dx, dy, num_ss,
                    quartile[0]);
      cmatrix S_matrix =
          cal_S_matrix(prm, buf_specENURT, w_specENURT, dx, dy, num_ss);
      output_result(ofs, prm, 2, max, mad, num_ss, t_s, t_e, S_matrix, flag2);
    }
  }

  for (int icmp = 0; icmp < 2; ++icmp) { // Horizontal components
    int l2 = search_max(ssRTU, px, py, event_number, mad, quartile, icmp);

    // Add new horizontal event based on the vertical component results
    // Simply scale the pxV and pyV
    for (int i1 = 0; i1 < l2V; i1++) {
      double p0v = sqrt(pxV[i1] * pxV[i1] + pyV[i1] * pyV[i1]);
      if (p0v > 0.04 && p0v * 1.8 < ipmax * dp) {
        double px_new = pxV[i1] * 1.8;
        double py_new = pyV[i1] * 1.8;

        bool isDuplicate = false;
        for (int i2 = 0; i2 < l2; i2++) {
          double dist2 = (px_new - px[i2]) * (px_new - px[i2]) +
                         (py_new - py[i2]) * (py_new - py[i2]);
          if (dist2 < 1E-2) {
            isDuplicate = true;
            break;
          }
        }

        if (!isDuplicate) {
          px[l2] = px_new;
          py[l2] = py_new;
          l2++;
        }
      }
    }
    // Search for maximum events based on the slant stack results in horizontal
    // components.
    for (int i = 0; i < l2; i++) {
      double max =
          ssRTU[icmp][round(px[i] / dp + 1E-10)][round(py[i] / dp + 1E-10)];

      PARAM prm;
      prm.p = sqrt(px[i] * px[i] + py[i] * py[i]);
      prm.θ = atan2(py[i], px[i]);
      prm.Δ = M_PI / 2.;
      prm.dp_Δ = 0.;
      rotate_EN_RT(sta0, buf_specENURT, w_specENURT, prm);
      est_dist_grid(prm, buf_specENURT[icmp + 3], w_specENURT[icmp + 3], dx, dy,
                    num_ss);
      rotate_EN_RT(sta0, buf_specENURT, w_specENURT, prm);
      int flag2 = est_dist_grad(prm, buf_specENURT[icmp + 3],
                                w_specENURT[icmp + 3], dx, dy, num_ss);
      rotate_EN_RT(sta0, buf_specENURT, w_specENURT, prm);
      if (flag2 > 1) {
        est_dist_boot(prm, buf_specENURT[icmp + 3], w_specENURT[icmp + 3], dx,
                      dy, num_ss, quartile[0]);
        cmatrix S_matrix =
            cal_S_matrix(prm, buf_specENURT, w_specENURT, dx, dy, num_ss);
        output_result(ofs, prm, icmp, max, mad, num_ss, t_s, t_e, S_matrix,
                      flag2);
      }
    }
  }
  return (num_ss);
}

void output_result(std::ostream &ofs, const PARAM &prm, int icmp, double max,
                   double mad, int num_ss, const ptime &t_s, const ptime &t_e,
                   cmatrix &S_matrix, int flag2) {

  double lat_s, lon_s;
  geod.ArcDirect(STATION::lat_ary, STATION::lon_ary, 90 - prm.θ / M_PI * 180.,
                 prm.Δ / M_PI * 180., lat_s, lon_s);

  ofs << to_iso_string(t_s) << " " << icmp << " " << prm.Δ / M_PI * 180 << " "
      << lat_s << " " << lon_s << " " << std::cos(prm.θ) * prm.p << " "
      << std::sin(prm.θ) * prm.p << " " << prm.p / prm.dp_Δ << " " << prm.dp
      << " " << -prm.dD / M_PI * 180 << " " << prm.dS << " " << max << " "
      << mad << " " << num_ss << " " << prm.ε << " " << prm.fmax << " "
      << STATION::lat_ary << " " << STATION::lon_ary << " "
      << to_iso_string(t_e) << " " << flag2 << " ";
  ofs << prm.cov[0][0] << " " << prm.cov[0][1] << " " << prm.cov[0][2] << " "
      << prm.cov[0][3] << " " << prm.cov[1][1] << " " << prm.cov[1][2] << " "
      << prm.cov[1][3] << " " << prm.cov[2][2] << " " << prm.cov[2][3] << " "
      << prm.cov[3][3] << " " << prm.S_boot << " " << prm.sigma << " "
      << prm.MS_max;
  for (int icmp2 = 0; icmp2 < 3; ++icmp2)
    ofs << " " << real(S_matrix(icmp2, icmp2));
  ofs << " " << real(S_matrix(0, 2));
  ofs << " " << imag(S_matrix(0, 2));
  ofs << std::endl;
}

/**
 * @brief Rotates the spectral data in the E-W and R-T planes.
 *
 * This function rotates the spectral data in the E-W and R-T planes based on
 * the azimuth and inclination angles.
 *
 * @param sta0 A reference to a vector of STATION objects containing seismic
 * data.
 * @param buf_specENURT A reference to a 4D array of complex numbers containing
 * the spectral data.
 * @param w_specENURT A reference to a 3D array of weights for the spectral
 * data.
 * @param prm The parameters to be updated.
 * @return int The status of the operation.
 */
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

/**
 * @brief Estimates the distance and updates the parameters by grid search.
 *
 * @param prm The parameters to be updated.
 * @param buf_spec 3D array of complex numbers representing the spectral data.
 * @param w_spec 2D array of weights for the spectral data.
 * @param dx Vector of x-coordinates of the stations.
 * @param dy Vector of y-coordinates of the stations.
 * @param num_ss Number of slant stacks.
 * @return int Number of iterations performed during the estimation process.
 */
static int est_dist_grid(PARAM &prm, const array3c &buf_spec,
                         const array2d &w_spec, const dvector &dx,
                         const dvector &dy, const int num_ss) {
  double maxS = 0;

  prm.Δ = -1;
  prm.dp_Δ = 0;
  {
    const int ddeg = 5;
    dvector S(180 / ddeg);
    std::vector<PARAM> prm1(180 / ddeg);
#pragma omp parallel for
    for (int ideg = 0; ideg < 180 / ddeg - 1; ideg++) {
      prm1[ideg] = prm;
      prm1[ideg].Δ = (ideg + 1) * ddeg / 180. * M_PI;
      S[ideg] = cal_S(prm1[ideg], buf_spec, w_spec, dx, dy, num_ss, 0);
    }

    for (int ideg = 0; ideg < 180 / ddeg - 1; ideg++) {
      if (maxS < S[ideg]) {
        maxS = S[ideg];
        prm = prm1[ideg];
      }
    }
    std::cerr << "#deg max=" << prm.Δ / M_PI * 180 << " " << prm.Δ << std::endl;
  }
  {
    maxS = 0.;
    dvector S(40);
    std::vector<PARAM> prm1(40);
    double dp_Δ0, dp_Δ1;
    if (prm.p > 0.04) {
      dp_Δ0 = -2E-5;
      dp_Δ1 = 1E-5;
    } else {
      dp_Δ0 = -2E-5;
      dp_Δ1 = 2E-5;
    }

#pragma omp parallel for
    for (int i = 0; i < 40; i++) {
      prm1[i] = prm;
      prm1[i].dp_Δ = (dp_Δ1 - dp_Δ0) * i / 40. +
                     dp_Δ0; // prm1[i].dp_Δ = (i-20)* .04/ (30.*111)/30.;
      S[i] = cal_S(prm1[i], buf_spec, w_spec, dx, dy, num_ss, 0);
    }

    for (int i = 0; i < 40; i++) {
      if (maxS < S[i]) {
        maxS = S[i];
        prm = prm1[i];
      }
    }
    std::cerr << "#dp_Δ max=" << prm.dp_Δ << std::endl;
  }
  return 0;
}

/**
 * @brief Estimates the distance and updates the parameters by Newton's method.
 *
 * @param prm The parameters to be updated.
 * @param buf_spec 3D array of complex numbers representing the spectral data.
 * @param w_spec 2D array of weights for the spectral data.
 * @param dx Vector of x-coordinates of the stations.
 * @param dy Vector of y-coordinates of the stations.
 * @param num_ss Number of slant stacks.
 * @return int Number of iterations performed during the estimation process.
 */
static int est_dist_grad(PARAM &prm, const array3c &buf_spec,
                         const array2d &w_spec, const dvector &dx,
                         const dvector &dy, const int num_ss) {
  double S0 = cal_S(prm, buf_spec, w_spec, dx, dy, num_ss, 0);
  PARAM prm0 = prm, prm_init, prm_tmp;
  dvector dS(4, 0.);
  dmatrix ddS(4, 4, 0.);

  int num_loop = 0;
  int num_eig = 0;
  if (prm.Δ > 0) {
    double Sinit = S0;
    double ε = 0;
    Eigen::Matrix4f S_tmp;
    Eigen::Vector4f d, iΛ, dprm, W;
    W << 0.06, M_PI / 2, M_PI / 2, .04 / (30. * 111);

    prm0 = prm;
    S0 = cal_S(prm0, buf_spec, w_spec, dx, dy, num_ss, 0);
    double S1 = S0 / 2;

    for (int i = 0; i < 20; ++i) {
      cal_HessianS(prm0, buf_spec, w_spec, dx, dy, dS, ddS, num_ss);
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          ddS(k, l) *= (W(l) * W(k));
          S_tmp(k, l) = ddS(k, l);
        }
        dS(k) *= W(k);
        d(k) = dS(k);
      }
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix4f> eigensolver(S_tmp);
      Eigen::Vector4f Λ = eigensolver.eigenvalues();
      Eigen::Matrix4f Q = eigensolver.eigenvectors();
      Eigen::Vector4f dd = Q.transpose() * d;

      num_eig = 0;
      for (int k = 0; k < 4; k++) {
        if (Λ(k) > 0) {
          iΛ(k) = -dd(k) / Λ(k) / 2.; // Switch to the steepest descent method
        } else {
          iΛ(k) = dd(k) / Λ(k); // Newton's method
          num_eig++;
        }
      }
      prm.num_eig = num_eig; // Return the number of negative eigenvalues
      dprm = W.asDiagonal() * Q * iΛ;

      bool flag_loop = 0;
      for (double r = 1.0; r > 1E-2; r *= 0.8) {
        prm_tmp = prm0;
        prm_tmp.p -= r * dprm(0);    // dS(0);
        prm_tmp.θ -= r * dprm(1);    // dS(1);
        prm_tmp.Δ -= r * dprm(2);    // dS(2);
        prm_tmp.dp_Δ -= r * dprm(3); // dS(3);

        S1 = cal_S(prm_tmp, buf_spec, w_spec, dx, dy, num_ss, 0);
        ε = ((S1 - S0) / S0);
        if (S1 > S0) {
          flag_loop = 1;
          break;
        }
      }
      if (flag_loop == 1) {
        S0 = S1;
        prm0 = prm_tmp;
      } else
        return (-1);

      if (abs(ε) < 1E-9 && dprm(2) * 180. / M_PI < 0.1) {
        num_loop = i + 1;
        break;
      }
    }
    if (num_loop == 0 || num_eig != 4)
      return (-1);
    double fmax = est_fmax(prm0, buf_spec, w_spec, dx, dy, num_ss);
    prm = prm0;
    prm.dS = S0 / Sinit;
    prm.ε = ε;
    prm.fmax = fmax;

    prm.dp = sqrt(pow(prm.p * cos(prm.θ) - prm_init.p * cos(prm_init.θ), 2) +
                  pow(prm.p * sin(prm.θ) - prm_init.p * sin(prm_init.θ), 2));
    // Recalculate the arc distance as to be consistent with the distance
    // between 0 and 180 degrees
    if (tan(prm.Δ) > 0)
      prm.Δ = atan(tan(prm.Δ));
    else
      prm.Δ = atan(tan(prm.Δ)) + M_PI;
    prm.dD = prm.Δ - prm_init.Δ;

    return (num_loop);
  } else
    return (-1);
}

/**
 * @brief Estimates the cov matrix of parameters using bootstrap method.
 *
 * @param prm Reference to a PARAM structure that holds the parameters to be
 * updated.
 * @param buf_spec 3D array of complex numbers representing the spectral data.
 * @param w_spec 2D array of weights for the spectral data.
 * @param dx Vector of x-coordinates of the stations.
 * @param dy Vector of y-coordinates of the stations.
 * @param num_ss Number of slant stacks.
 * @return int Number of iterations performed during the estimation process.
 */
static int est_dist_boot(PARAM &prm, const array3c &buf_spec,
                         const array2d &w_spec, const dvector &dx,
                         const dvector &dy, const int num_ss,
                         const double bias) {
  // Estimate the cov. matrix
  const int loop_num = 100;
  dvector S(loop_num);
  const int num_sta = buf_spec.shape()[1];
  (void)bias;

  Eigen::Matrix4f S_tmp;
  dvector dS(4, 0.);
  dmatrix ddS(4, 4, 0.);
  Eigen::Vector4f W;
  W << 0.06, M_PI / 2, M_PI / 2, .04 / (30. * 111);

  cal_HessianS(prm, buf_spec, w_spec, dx, dy, dS, ddS, num_ss);

  // #pragma omp parallel for
  for (int i = 0; i < loop_num; i++) {
    boost::mt19937 gen(static_cast<unsigned long>(time(0) * i));
    boost::uniform_smallint<> dst(0, num_sta - 1);
    boost::variate_generator<boost::mt19937 &, boost::uniform_smallint<>> rand(
        gen, dst);
    dvector bootstrap(num_sta);
    array2d w_bootstrap(boost::extents[num_ss][num_sta]);

    for (int isg = 0; isg < num_ss; isg++) {
      // Bootstrap resampling
      for (int ist = 0; ist < num_sta; ist++)
        bootstrap[ist] = 0;
      int num_tmp = 0;
      for (int ist = 0; ist < num_sta; ist++)
        if (w_spec[isg][ist] != 0)
          num_tmp++;
      int ist = 0;
      while (1) {
        int bsmp = rand();
        if (w_spec[isg][bsmp] > 0) {
          bootstrap[bsmp]++;
          ist++;
        }
        if (ist == num_tmp)
          break;
      }
      for (int ist = 0; ist < num_sta; ist++)
        w_bootstrap[isg][ist] = w_spec[isg][ist] * bootstrap[ist];
    }

    S[i] = cal_S(prm, buf_spec, w_bootstrap, dx, dy, num_ss, 1);
  }
  double mean = 0, sigma = 0;
  for (int i = 0; i < loop_num; i++)
    mean += S[i];
  mean /= loop_num;
  for (int i = 0; i < loop_num; i++)
    sigma += pow(S[i] - mean, 2);
  sigma = sqrt(sigma / loop_num);

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      S_tmp(i, j) = ddS(i, j) / sigma;
  Eigen::Matrix4f Cov = S_tmp.inverse();
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      prm.cov[i][j] = -Cov(i, j) * W(i) * W(j);

  /// Estimation of max beam power without the bias
  double S_est = cal_S(prm, buf_spec, w_spec, dx, dy, num_ss, 1);
  double weight = 0.;
  for (int isg = 0; isg < num_ss; isg++) {
    double sum1 = 0;
    double sum2 = 0;
    for (int ist = 0; ist < num_sta; ist++) {
      sum1 += (w_spec[isg][ist]);
      sum2 += (w_spec[isg][ist] * w_spec[isg][ist]);
    }
    weight += (sum1 * sum1 - sum2);
  }
  prm.MS_max = S_est / weight * (STATION::df);
  prm.S_boot = mean / weight * (STATION::df);
  prm.sigma = sigma / weight * (STATION::df);

  return (0);
}

/**
 * @brief Calculates the beamforming power S based on station delays and
 * spectral data.
 *
 * This function computes the total power (S) by stacking of spectral data from
 * multiple stations according to the slowness parameters (p, θ, D, dp_Δ)
 * held in \p prm. The arrays \p buf_spec and \p w_spec store spectral data and
 * weights, while \p dx and \p dy hold station coordinates relative to a
 * reference point. If \p flag_red is set, a bias term is subtracted from the
 * final power value.
 *
 * @param prm      The parameter set containing slowness, azimuth, and related
 * fields.
 * @param buf_spec A 3D complex array holding spectral data for each station and
 * frequency.
 * @param w_spec   A 2D array storing the corresponding weights of the spectral
 * data.
 * @param dx       A vector of station x-coordinates relative to the array
 * center.
 * @param dy       A vector of station y-coordinates relative to the array
 * center.
 * @param num_ss   The number of slant-stack segments to process.
 * @param flag_red A flag indicating whether to subtract the bias term from the
 * result.
 * @return The calculated beamforming power (S).
 */
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

/**
 * @brief Calculates the matrix of the beamforming power S based on station
 * delays and spectral data.
 *
 * This function computes the total power (S) by stacking of spectral data from
 * multiple stations according to the slowness parameters (p, θ, D, dp_Δ)
 * held in \p prm. The arrays \p buf_spec and \p w_spec store spectral data and
 * weights, while \p dx and \p dy hold station coordinates relative to a
 * reference point.
 *
 * @param prm      The parameter set containing slowness, azimuth, and related
 * fields.
 * @param buf_specENURT A 4D complex array holding spectral data for each
 * station and frequency.
 * @param w_spec   A 2D array storing the corresponding weights of the spectral
 * data.
 * @param dx       A vector of station x-coordinates relative to the array
 * center.
 * @param dy       A vector of station y-coordinates relative to the array
 * center.
 * @param num_ss   The number of slant-stack segments to process.
 * @return The calculated a matrix of beamforming power (S_RTU) in RTU cmps.
 */
static cmatrix cal_S_matrix(const PARAM prm, const array4c &buf_specENURT,
                            const array3d &w_spec, const dvector &dx,
                            const dvector &dy, const int num_ss) {
  const int num_sta = buf_specENURT.shape()[2];
  cmatrix phi(3, STATION::if2 + 1, 0);

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

  cmatrix S_RTU(3, 3, 0.);
  dmatrix weight(3, 3, 0.);
  const int idx_RTU[3] = {3, 4, 2}; // ENURT => RTU

  for (int ibuf = 0; ibuf < num_ss; ibuf++) {
    for (int icmp = 0; icmp < 3; icmp++) {
      for (int k = STATION::if1; k <= STATION::if2; ++k)
        phi(icmp, k) = 0;
      for (int i = 0; i < num_sta; i++) {
        double phase = tau[i] * 2. * M_PI * STATION::df;
        double dc = cos(phase);
        double ds = sin(phase);
        double cp0 = cos(phase * STATION::if1);
        double sp0 = sin(phase * STATION::if1);

        for (int k = STATION::if1; k <= STATION::if2; ++k) {
          double cp1 = cp0 * dc - sp0 * ds;
          double sp1 = sp0 * dc + cp0 * ds;
          phi(icmp, k) +=
              (buf_specENURT[idx_RTU[icmp]][ibuf][i][k] *
               std::complex<double>(cp0, sp0) *
               w_spec[idx_RTU[icmp]][ibuf][i]); //*(double)flag_ss[ist]);
          cp0 = cp1;
          sp0 = sp1;
        }
      }
    }
    // Compute the bias of the matrix S_RTU
    // S_RTU = phi^* phi = phi^* phi - bias
    for (int icmp1 = 0; icmp1 < 3; icmp1++) {
      for (int icmp2 = 0; icmp2 < 3; icmp2++) {
        double bias = 0;
        for (int ist = 0; ist < num_sta; ist++) {
          for (int k = STATION::if1; k <= STATION::if2; ++k) {
            bias += real(conj(buf_specENURT[idx_RTU[icmp1]][ibuf][ist][k]) *
                         buf_specENURT[idx_RTU[icmp2]][ibuf][ist][k] *
                         w_spec[idx_RTU[icmp1]][ibuf][ist] *
                         w_spec[idx_RTU[icmp2]][ibuf][ist]);
          }
        }
        S_RTU(icmp1, icmp1) -= bias;
      }
    }
    // Compute the weight matrix
    // weight = sum(phi) sum(phi)^* - sum(phi phi^*)
    dvector sum1(3, 0.);
    dmatrix sum2(3, 3, 0.);
    for (int icmp1 = 0; icmp1 < 3; icmp1++) {
      for (int ist = 0; ist < num_sta; ist++) {
        sum1(icmp1) += (w_spec[idx_RTU[icmp1]][ibuf][ist]);
      }
    }
    for (int icmp1 = 0; icmp1 < 3; icmp1++) {
      for (int icmp2 = 0; icmp2 < 3; icmp2++) {
        for (int k = STATION::if1; k <= STATION::if2; ++k) {
          S_RTU(icmp1, icmp2) += (conj(phi(icmp1, k)) * phi(icmp2, k));
        }
        for (int ist = 0; ist < num_sta; ist++) {
          sum2(icmp1, icmp2) += (w_spec[idx_RTU[icmp1]][ibuf][ist] *
                                 w_spec[idx_RTU[icmp2]][ibuf][ist]);
        }
        weight(icmp1, icmp2) += (sum1[icmp1] * sum1[icmp2]);
        weight(icmp1, icmp2) -= (sum2(icmp1, icmp2));
      }
    }
  }
  // Compute the matrix S_RTU
  for (int icmp1 = 0; icmp1 < 3; icmp1++) {
    for (int icmp2 = 0; icmp2 < 3; icmp2++) {
      S_RTU(icmp1, icmp2) =
          (S_RTU(icmp1, icmp2) / weight(icmp1, icmp2) /
           (double)(STATION::if2 - STATION::if1 + 1)); // * (STATION::df));
    }
  }
  return (S_RTU);
}

static double est_fmax(const PARAM prm, const array3c &buf_spec,
                       const array2d &w_spec, const dvector &dx,
                       const dvector &dy, const int num_ss) {
  const int num_sta = buf_spec.shape()[1];
  cvector phi(STATION::nfreq, 0.);
  dvector psd(STATION::nfreq, 0.);

  double p = prm.p;
  double ex = cos(prm.θ);
  double ey = sin(prm.θ);

  dvector tau(num_sta, 0);
  for (int i = 0; i < num_sta; i++) {
    double η = (ex * dx[i] + ey * dy[i]) / 6371.;
    double ζ = (-ey * dx[i] + ex * dy[i]) / 6371.;
    double cotΔ = cos(prm.Δ) / sin(prm.Δ);
    double dp_Δ = prm.dp_Δ;
    double l = (-η + ζ * ζ * cotΔ / 2.) * 6371.;
    tau[i] = l * (p + dp_Δ * l / 2);
  }

  for (int ibuf = 0; ibuf < num_ss; ibuf++) {
    for (int k = STATION::if1; k < STATION::nfreq; ++k)
      phi[k] = 0;
    for (int i = 0; i < num_sta; i++) {
      double phase = tau[i] * 2. * M_PI * STATION::df;
      double dc = cos(phase);
      double ds = sin(phase);
      double cp0 = cos(phase * STATION::if1);
      double sp0 = sin(phase * STATION::if1);

      for (int k = STATION::if1; k < STATION::nfreq; ++k) {
        double cp1 = cp0 * dc - sp0 * ds;
        double sp1 = sp0 * dc + cp0 * ds;
        phi[k] += (buf_spec[ibuf][i][k] * std::complex<double>(cp0, sp0) *
                   w_spec[ibuf][i]); //*(double)flag_ss[ist]);
        cp0 = cp1;
        sp0 = sp1;
      }
    }
    for (int k = STATION::if1; k < STATION::nfreq; ++k)
      psd[k] += real(conj(phi[k]) * phi[k]);
  }
  auto itr = max_element(psd.begin(), psd.end());
  size_t ifmax = distance(psd.begin(), itr);

  return (ifmax * STATION::df);
}

/**
 * @brief Calculates the Hessian matrix of the objective function.
 *
 * @param prm The parameters to be updated.
 * @param buf_spec 3D array of complex numbers representing the spectral data.
 * @param w_spec 2D array of weights for the spectral data.
 * @param dx Vector of x-coordinates of the stations.
 * @param dy Vector of y-coordinates of the stations.
 * @param dS Vector of the first derivatives of the objective function.
 * @param ddS Matrix of the second derivatives of the objective function.
 * @param num_ss Number of slant stacks.
 */
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

int eval_deri(std::vector<STATION> &sta0, double th, double &median0,
              double &median1) {
  double deri, integ_cur;
  SPCTRM specE, specN, specZ;

  for (int ist = 0; ist < (int)sta0.size(); ist++) { // Calculation of spectra
    sta0[ist].print_spec(specE, specN, specZ);

    tmp_integ0[ist] = specZ.integ[0];
    tmp_integ1[ist] = specZ.integ[1];
  }

  std::sort(tmp_integ0, tmp_integ0 + sta0.size());
  std::sort(tmp_integ1, tmp_integ1 + sta0.size());

  integ_cur = tmp_integ1[sta0.size() / 2]; // For Test
  median0 = tmp_integ0[sta0.size() / 2];
  median1 = tmp_integ1[sta0.size() / 2];

  if (integ_cur != 0 && integ_pre != 0)
    deri = fabs(log(integ_cur / integ_pre));
  else
    deri = 100;
  // fprintf(stderr,"###Der %g %g %g\n",deri,integ_old,integ_cur);
  integ_pre = integ_cur;
  if (integ_old == 0 && deri < log(th))
    integ_old = integ_cur;
  // if(deri<th&&integ_cur<integ_old*2&&integ_cur>integ_old/2){
  if (deri < log(th) &&
      log(integ_cur) < log(integ_old) + log(th) * (count_gap + 1) &&
      log(integ_cur) > log(integ_old) - log(th) * (count_gap + 1)) {
    integ_old = integ_cur;
    count_gap = 0;
    return (1);
  } else {
    count_gap++;
    return (0);
  }
}

int search_max(const array3d &ssRTU, dvector &px, dvector &py,
               const int event_number, double &mad,
               std::vector<double> &quartile, const int icmp) {
  const int number = (2 * ipmax + 1) * (2 * ipmax + 1);
  const int imask = int(1.E-2 / dp + 0.5);
  std::vector<int> dpy_ary(imask + 1);
  for (int dpx = 0; dpx <= imask; dpx++) {
    dpy_ary[dpx] = round(sqrt(imask * imask - dpx * dpx));
  }
  int num_mask = 0;
  for (int dpx = -imask; dpx <= imask; dpx++) {
    for (int dpy = -dpy_ary[abs(dpx)]; dpy <= dpy_ary[abs(dpx)]; dpy++) {
      num_mask++;
    }
  }
  array2b mask(boost::extents[range2d(-ipmax - imask, ipmax + imask + 1)]
                             [range2d(-ipmax - imask, ipmax + imask + 1)]);
  std::fill_n(mask.data(), mask.num_elements(), true);

  struct ss_data {
    double val;
    int ipx, ipy;
  };
  std::vector<ss_data> ss1(number);
  dvector sz(number), s_tmp(number);

  int count = 0;
  for (int i = -ipmax; i <= ipmax; i++) {
    for (int j = -ipmax; j <= ipmax; j++) {
      ss1[count].val = ssRTU[icmp][i][j];
      ss1[count].ipx = i;
      ss1[count].ipy = j;
      s_tmp[count] = ssRTU[icmp][i][j];
      count++;
      if (sqrt(i * i + j * j) * dp > ipmax * dp)
        mask[i][j] = false;
    }
  }
  // Calculation of MAD
  mad = boost::math::statistics::median_absolute_deviation(
      s_tmp); // This function changes the order of s_tmp
  std::sort(s_tmp.begin(), s_tmp.end());
  // Calculation of quartiles
  quartile[0] = s_tmp[(int)(number / 4)];     // 1st quartile
  quartile[1] = s_tmp[(int)(number / 2)];     // Median
  quartile[2] = s_tmp[(int)(number * 3 / 4)]; // 3rd quartile

  std::cerr << "#$" << mad << " " << quartile[0] << " " << icmp << std::endl;
  std::sort(ss1.begin(), ss1.end(), [](const ss_data &a, const ss_data &b) {
    return (a.val) > (b.val);
  });

  // for (int i = 0; i < number; i++)
  //   ss1[i].val -= s_tmp[0]; // Subtract the minimum
  //  ss1[i].val -= quartile[0]; // Subtract the 1st quartile to reduce the
  //  bias.

  int idx = 0;
  for (int i = 0; i < number; i++) {
    // Because ss[i].val - base < ss[i].val - min(ss1) < mad*5,
    // ss[i+j].val - base must be smaller than mad*5 for j>=i.
    // For this reason, we can break the loop.
    if (ss1[i].val - ss1.back().val < mad * 7)
      break;
    if (mask[ss1[i].ipx][ss1[i].ipy] == true) {
      // Estimate the base value by the average of the surrounding values
      double base = 0.;
      int num_base = 0;
      for (int dpx = -imask; dpx <= imask; dpx++) {
        int k = ss1[i].ipx + dpx; // clang-format off
				base += ssRTU[icmp][k][ss1[i].ipy - dpy_ary[abs(dpx)]];
				base += ssRTU[icmp][k][ss1[i].ipy + dpy_ary[abs(dpx)]];
				num_base += 2;
			}	
			for (int dpy = -dpy_ary[imask]+1; dpy <= dpy_ary[imask]-1; dpy++) {
				base += ssRTU[icmp][ss1[i].ipx - imask][dpy];
				base += ssRTU[icmp][ss1[i].ipx + imask][dpy];
				num_base += 2;
			}
			base /= (double)(num_base);
			//Here, we roughly select the peaks by the difference from the minumum
      if (ss1[i].val - base < mad * 7)
          continue;

      int num_peak = 0;
      for (int dpx = -imask; dpx <= imask; dpx++) {
        int k = ss1[i].ipx + dpx; // clang-format off
				if (k < -ipmax || k > ipmax) continue;
        for (int l = std::max(ss1[i].ipy - dpy_ary[abs(dpx)],-ipmax); 
                l <= std::min(ss1[i].ipy + dpy_ary[abs(dpx)], ipmax); l++) {                             // clang-format on
          if (ss1[i].val > ssRTU[icmp][k][l]) // quartile[0])
            num_peak++;
        }
      }
      if (num_peak != num_mask - 1)
        break;
      px[idx] = ss1[i].ipx * dp;
      py[idx] = ss1[i].ipy * dp;
      for (int dpx = -imask; dpx <= imask; dpx++) {
        int k = ss1[i].ipx + dpx;
        for (int l = ss1[i].ipy - dpy_ary[abs(dpx)];
             l <= ss1[i].ipy + dpy_ary[abs(dpx)]; l++) {
          mask[k][l] = false;
        }
      }
      idx++;
    }
    if (idx >= event_number) {
      break;
    }
  }
  // std::cout << idx << " " << idx_peaks << std::endl;
  return (idx);
}
