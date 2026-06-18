#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <boost/format.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <Eigen/Dense>
#define BOOST_DISABLE_ASSERTS

using namespace std;
using namespace GeographicLib;
using namespace boost::posix_time;
using namespace boost::gregorian;

namespace fs = boost::filesystem;
namespace ublas = boost::numeric::ublas; //Consider atd::vector

static const int yr0 = 2004; //2015;
static double dp = 2.5E-3;   //.05E-2;
static int ipmax = (int)(1.5E-1 / dp);
static double px0 = 0;      //-0.007;
static double py0 = 0;      // 0.047;
static double rad0 = 1500.; // The radius of the initial beamforming
static double rad1 = 1500.; // The radius of the array is 500 km.;
static const Geodesic &geod = Geodesic::WGS84();

typedef boost::multi_array<complex<double>, 3> array3c;
typedef array3c::extent_range range3c;
typedef boost::multi_array<double, 3> array3d;
typedef array3d::extent_range range3d;
typedef boost::multi_array<double, 2> array2d;
typedef array2d::extent_range range2d;
typedef boost::multi_array<bool, 2> array2b;
typedef array2b::extent_range range2b;

typedef ublas::vector<double> dvector;
typedef ublas::vector<complex<double>> cvector;
typedef ublas::matrix<double> dmatrix;
typedef ublas::matrix<complex<double>> cmatrix;

#include <math.h>
#include "station_info.h"
#include "util.h"
#include "calTT.h"
#include <omp.h>

struct PARAM
{
  double p, theta, D, dp_l, fmax; //p [s/km], theta [rad], dp_l [s/km/km]
  double dp, dD, dS, epsilon;
};

static int slant_stack(vector<STATION> &sta0, array2d &ssZ);
static double cal_S(const PARAM prm, const array3c &buf_specZ, const array2d &w_specZ,
                    const dvector &dx, const dvector &dy, const int num_ss);
static double est_fmax(const PARAM prm, const array3c &buf_specZ, const array2d &w_specZ,
                       const dvector &dx, const dvector &dy, const int num_ss);
static double cal_matS(const PARAM prm, const array3c &buf_specZ, const array2d &w_specZ, const dvector &dx, const dvector &dy,
                       dvector &dS, dmatrix &ddS, const int num_ss);
static int est_dist(PARAM &prm, const array3c &buf_specZ, const array2d &w_specZ, const dvector &dx, const dvector &dy, const int num_ss);
static int read_station_correction(vector<STATION> &sta0);

static int search_max(const array2d &ssZ, dvector &px, dvector &py, const int lmax, double &mad);

//For definition of transients
int eval_deri(vector<STATION> &sta0, double th, double &median0, double &median1);
int double_cmp(const double *a, const double *b)
{
  if (*a < *b)
    return (-1);
  else if (*a > *b)
    return (1);
  else
    return (0);
}
static double tmp_integ0[2000], tmp_integ1[2000], integ_pre = 0, integ_old = 0;
vector<ptime> CMT_ptime;
vector<double> CMT_slat,CMT_slon,CMT_sdep,CMT_moment;

int main(int argc, char *argv[])
{
  STATION::dt_msec = 500;    //Sampling interval in millisecond
  STATION::len = 1024 * 2;   //2^x
  STATION::stride = 667 * 2; //667*2*128 = (86400*2-1024*2) //824;
  STATION::df = 1. / (STATION::len * STATION::dt_msec * 1E-3);
  STATION::nfreq = (int)(2.6E-1 / STATION::df); //STATION::len/2;
  STATION::npts = 86400 * 1000 / (STATION::dt_msec);
  STATION::init_Freq();
  STATION::nl_h = 50;
  STATION::nl_v = 50;
  //STATION::flag_ary=0;
  STATION::if1 = (int)(0.1 / STATION::df);
  STATION::if2 = (int)(0.2 / STATION::df);

  //Cache setting
  // Make a list of file access properties
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  size_t rdcc_nslots = 1049;             // Hash entry
  size_t rdcc_nbytes = 2147483648;       // 2GB cache
  double rdcc_w0 = 0.75;                 // Balance of writing 
  H5Pset_cache(fapl, 0, rdcc_nslots, rdcc_nbytes, rdcc_w0);

  if(argc != 3){
    cerr << "Usage:   ./bin/cal_ccf_eq <hinet-root> <cmt-catalog>" << endl;
    H5Pclose(fapl);
    return(1);
  }
  string dir_Hinet = argv[1];
  if(!dir_Hinet.empty() && dir_Hinet[dir_Hinet.size() - 1] != '/')
    dir_Hinet += "/";
  const string cmt_catalog = argv[2];
  if(init_CMT(2004, cmt_catalog) != 0){
    H5Pclose(fapl);
    return(1);
  }
  get_CMT(CMT_ptime,CMT_slat,CMT_slon,CMT_sdep,CMT_moment);
  cout << "#" << STATION::if1 * STATION::df << " " << STATION::if2 * STATION::df << endl;

  //Definitions of variables for calculations of slant stack.
  array2d ssZ(boost::extents[range2d(-ipmax, ipmax + 1)][range2d(-ipmax, ipmax + 1)]);

  int sta_num = 0;
  int count0 = 0;
  for (int days = 0; days < 365 * 14.67; days += 1)
  {                                                //*14.67
                                                   //for(int days=0;days<1;days+=1){
    date d(date(yr0, 4, 1) + date_duration(days)); //4
    date d0(d.year(), 1, 1);
    string h5file = dir_Hinet +
                    (boost::format("%d/%02d%02d/%d%03d0000.h5") % d.year() % (int)(d.month()) % d.day() % d.year() % ((d - d0).days() + 1)).str();
    bool flag_eq = 0;
    ptime t1 = ptime(d, time_duration(0, 5, 0));
    ptime t2 = ptime(d, time_duration(23, 55, 0));
    for (std::size_t i = 0; i < CMT_ptime.size(); ++i)
      if (t1 < CMT_ptime[i] && CMT_ptime[i] < t2)
        flag_eq = 1;

    const fs::path path(h5file);
    if (fs::exists(path) && flag_eq == 1)
    {
      cerr << t1 << endl;

      vector<STATION> sta0;
      init_station(h5file, sta0, rad0, rad1,fapl);
      read_station_correction(sta0);

      //Clear sac information
      for (int i = 0; i < (int)sta0.size(); i++)
        sta0[i].clear_sac();
      //fill_n(ssZ.data(), ssZ.num_elements(), 0.);

      sta_num = read_h5(sta0, h5file,fapl);
      cerr << "#Read data " << h5file << " " << sta_num << endl;
      if (sta_num > 0)
        count0 += slant_stack(sta0, ssZ);
    }
  }
  H5Pclose(fapl);
  return(0);
}

//////
static int slant_stack(vector<STATION> &sta0, array2d &ssZ)
{
  ptime t1, t2;
  SPCTRM specE, specN, specZ;
  int count0 = 0;
  double median0, median1;
  //Determine reference time t0
  ptime t0;
  bool flag = 0;
  for (int i = 0; i < (int)sta0.size(); i++)
  {
    if (sta0[i].print_sta_num() == 3)
    {
      if (flag == 0)
        t0 = sta0[i].print_sacE().ts;
      if (t0 > sta0[i].print_sacE().ts)
      {
        t0 = sta0[i].print_sacE().ts;
        flag = 1;
      }
    }
  }

  //Definition for buffer of specZ.
  const int ibuf = (STATION::npts - STATION::len) / STATION::stride + 1;
  array3c buf_specZ(boost::extents[ibuf][(int)sta0.size()][range3c(STATION::if1, STATION::nfreq)]);
  array2d w_specZ(boost::extents[ibuf][(int)sta0.size()]);
  array2c buf_aryZ(boost::extents[(int)sta0.size()][range3c(STATION::if1, STATION::if2 + 1)]);
  array3c sZ(boost::extents[range3c(-ipmax, ipmax + 1)][range3c(-ipmax, ipmax + 1)][range3c(STATION::if1, STATION::if2 + 1)]);
  dvector dx_ary((int)sta0.size()), dy_ary((int)sta0.size());

  long len = STATION::len * STATION::dt_msec;
  dvector msZ(STATION::if2 + 1, 0.);

  int num_ss = 0;
  vector<int> flag_nl((int)sta0.size(), 0);
  vector<int> flag_ss((int)sta0.size(), 0);

  int count_all = 0;
  for (int ist = 0; ist < (int)sta0.size(); ist++)
    if (sta0[ist].print_rad() < rad1)
      count_all++;

  for (std::size_t iev = 0; iev < CMT_ptime.size(); ++iev)
    if (t0 + minutes(5) < CMT_ptime[iev] && CMT_ptime[iev] + milliseconds(len) < t0 + minutes(1440))
    {
      fill_n(ssZ.data(), ssZ.num_elements(), 0.);
      //{cerr<<CMT_ptime[i]<<endl;}}exit(0);
      //for(long int_t=0;int_t+len<= 3600*24*1000 ;int_t += stride){
      t1 = CMT_ptime[iev] - minutes(2);
      t2 = t1 + milliseconds(len);
      for (int k = STATION::if1; k <= STATION::if2; ++k)
        msZ[k] = 0;
      fill_n(sZ.data(), sZ.num_elements(), 0.);

#pragma omp parallel for private(specE, specN, specZ)
      for (int ist = 0; ist < (int)sta0.size(); ist++)
      {
        flag_nl[ist] = 0;
        if (sta0[ist].print_sta_num() && sta0[ist].cal_spec(t1))
        {
          sta0[ist].print_spec(specE, specN, specZ);
          if (1.0 < specZ.integ[0] && specZ.integ[0] < 4E10 &&
              1.0 < specZ.integ[1] && specZ.integ[1] < 4E9 &&
              1.0 < specZ.integ[2] && specZ.integ[2] < 4E9)
          {
            //For debug
            /*if(ist%20==0){
              for(int isg=(t1-t0).total_seconds()*2;isg<(t1-t0).total_seconds()*2+STATION::len;isg++){
                cout << sta0[ist].print_lon()<<" "<<sta0[ist].print_lat()<<" "<<isg*.5<<" "<< sta0[ist].print_sacZ().sgram[isg] << endl;
              }
              cout << endl;
            }*/
            flag_nl[ist] = 1;
          }
        }
      }

      int count_ss = 0;
      int count_est = 0;
      for (int ist = 0; ist < (int)sta0.size(); ist++)
      {
        if (flag_nl[ist] == 1 && sta0[ist].print_rad() < rad1)
        {
          flag_ss[ist] = 1;
          count_ss++;
        }
        else
          flag_ss[ist] = 0;
        if (flag_nl[ist] == 1)
          count_est++;
      }
      int flag_deri = eval_deri(sta0, 0.03, median0, median1);
      cerr << count_ss / (double)count_all << " " << count_est / (double)sta0.size() << " " << flag_deri << " " << median0 << " " << median1 << endl;
      if (count_all * .8 && count_est > 0.8 * (int)(sta0.size()))
      {
        //Copy of specZ for further loop of the parameter search
        int num_ary = 0;
        for (int ist = 0; ist < (int)sta0.size(); ist++)
        {
          if (flag_nl[ist] == 1)
            w_specZ[num_ss][ist] = 1.;
          else
            w_specZ[num_ss][ist] = 0.;

          sta0[ist].print_spec(specE, specN, specZ);
          ///Travel time correction
          
          double tau = sta0[ist].print_dtau();
          double phase = tau * 2. * M_PI * STATION::df;
          double dc = cos(phase);
          double ds = sin(phase);
          double cp0 = cos(phase * STATION::if1);
          double sp0 = sin(phase * STATION::if1);

          for (int k = STATION::if1; k < STATION::nfreq; ++k)
          {
            double cp1 = cp0 * dc - sp0 * ds;
            double sp1 = sp0 * dc + cp0 * ds;
            specZ.spec[k] = (specZ.spec[k] * complex<double>(cp0, sp0)); 
            buf_specZ[num_ss][ist][k] = specZ.spec[k];
            cp0 = cp1;
            sp0 = sp1;
          }
          
          if (flag_ss[ist] == 1)
          {
            for (int k = STATION::if1; k <= STATION::if2; ++k)
              buf_aryZ[num_ary][k] = specZ.spec[k];
            dx_ary[num_ary] = sta0[ist].print_dx();
            dy_ary[num_ary] = sta0[ist].print_dy();

            num_ary++;
          }
        }
        num_ss++;
#pragma omp parallel for //collapse(2)//private(specE,specN,specZ)
        for (int ipx = -ipmax; ipx <= ipmax; ipx++)
        {
          //for(int ist=0;ist<count_ss;++ist){
          for (int ist = 0; ist < num_ary; ++ist)
          {
            for (int ipy = -ipmax; ipy <= ipmax; ipy++)
            {
              double px = px0 + ipx * dp;
              double py = py0 + ipy * dp;

              double tau = -(px * dx_ary[ist] + py * dy_ary[ist]);
              double phase = tau * 2. * M_PI * STATION::df;
              double dc = cos(phase);
              double ds = sin(phase);
              double cp0 = cos(phase * STATION::if1);
              double sp0 = sin(phase * STATION::if1);

              for (int k = STATION::if1; k <= STATION::if2; ++k)
              {
                double cp1 = cp0 * dc - sp0 * ds;
                double sp1 = sp0 * dc + cp0 * ds;
                sZ[ipx][ipy][k] += (buf_aryZ[ist][k] * complex<double>(cp0, sp0));
                cp0 = cp1;
                sp0 = sp1;
              }
            }
          }
        }

        for (int ist = 0; ist < (int)sta0.size(); ist++)
        {
          if (flag_ss[ist] == 1)
          {
            sta0[ist].print_spec(specE, specN, specZ);
            for (int k = STATION::if1; k <= STATION::if2; ++k)
              msZ[k] += norm(specZ.spec[k]);
          }
        }
        cerr << t1 << " " << median0 << " " << median1 << " " << flag_deri << " " << count_ss << endl;

        double normZ;
        for (int ipx = -ipmax; ipx <= ipmax; ipx++)
        {
          for (int ipy = -ipmax; ipy <= ipmax; ipy++)
          {
            normZ = 0;
            for (int k = STATION::if1; k <= STATION::if2; ++k)
              normZ += norm(sZ[ipx][ipy][k] / (double)(count_ss));
          }
        }
        count0++;
        for (int ipx = -ipmax; ipx <= ipmax; ipx++)
        {
          for (int ipy = -ipmax; ipy <= ipmax; ipy++)
          {
            for (int k = STATION::if1; k <= STATION::if2; ++k)
            {
              ssZ[ipx][ipy] += ((norm(sZ[ipx][ipy][k]) - msZ[k]) / count_ss / (STATION::if2 - STATION::if1 + 1));
            }
          }
        }
      }
      if (num_ss == 0)
	return (0);
      double mad;
      dvector px(10), py(10);
      dvector dx((int)sta0.size()), dy((int)sta0.size());
      for (int ist = 0; ist < (int)sta0.size(); ist++)
	{
	  dx[ist] = sta0[ist].print_dx();
	  dy[ist] = sta0[ist].print_dy();
	}
      int l2 = search_max(ssZ, px, py, 10, mad);
      for (int i = 0; i < l2; i++)
	{
	  double s12;
	  double max = ssZ[round(px[i] / dp)][round(py[i] / dp)];
	  cerr << "#(px,py) =(" << px[i] << " " << py[i] << ")"
	       << ", num_ss=" << num_ss << " " << max << endl;
	  PARAM prm;
	  prm.p = sqrt(px[i] * px[i] + py[i] * py[i]);
	  prm.theta = atan2(py[i], px[i]);
	  prm.D = M_PI / 2.;
	  prm.dp_l = 0.;
	  int flag2 = est_dist(prm, buf_specZ, w_specZ, dx, dy, num_ss);
	  if (flag2 > 0)
	    {
	      double lat_s, lon_s;
	      geod.ArcDirect(STATION::lat_ary, STATION::lon_ary, 90 - prm.theta / M_PI * 180., prm.D / M_PI * 180., lat_s, lon_s);
	      double arc = geod.Inverse(STATION::lat_ary,STATION::lon_ary,CMT_slat[iev],CMT_slon[iev],s12);//
	      cout << t1.date() << " " << prm.D / M_PI * 180 << " "
		   << lat_s << " " << lon_s << " "
		   << cos(prm.theta) * prm.p << " " << sin(prm.theta) * prm.p << " " << prm.p / prm.dp_l << " "
		   << " " << prm.dp << " " << -prm.dD / M_PI * 180 << " " << prm.dS << " "
		   << max << " " << mad << " " << num_ss << " " << prm.epsilon << " " << prm.fmax << " "
		   << STATION::lat_ary << " " << STATION::lon_ary<<" "
		   << CMT_ptime[iev]<<" "<<CMT_slat[iev]<<" "<<CMT_slon[iev]<<" "<<CMT_sdep[iev]<<" "<<CMT_moment[iev] <<" "<<arc<<" "<<i<<endl;
	    }
	}
    }
  return (count0);
}
int est_dist(PARAM &prm, const array3c &buf_specZ, const array2d &w_specZ,
             const dvector &dx, const dvector &dy, const int num_ss)
{
  double S0 = cal_S(prm, buf_specZ, w_specZ, dx, dy, num_ss);
  double maxS = 0;
  PARAM prm0 = prm, prm_init, prm_tmp;

  //ublas::permutation_matrix<> pm(4);
  dvector dS(4, 0.), p0(4), p1(4);
  dmatrix ddS(4, 4, 0.);

  prm.D = -1;
  prm.dp_l = 0;

  {
    const int ddeg = 5;
    dvector S(180 / ddeg);
    vector<PARAM> prm1(180 / ddeg);
#pragma omp parallel for
    for (int ideg = 0; ideg < 180 / ddeg - 1; ideg++)
    {
      prm1[ideg] = prm;
      prm1[ideg].D = (ideg + 1) * ddeg / 180. * M_PI;
      S[ideg] = cal_S(prm1[ideg], buf_specZ, w_specZ, dx, dy, num_ss);
    }

    for (int ideg = 0; ideg < 180 / ddeg - 1; ideg++)
    {
      if (maxS < S[ideg])
      {
        maxS = S[ideg];
        prm = prm1[ideg];
      }
    }
    prm_init = prm;
    cerr << "#deg max=" << prm.D / M_PI * 180 << " " << prm.D << endl;
  }
  {
    maxS = 0.;
    dvector S(40);
    vector<PARAM> prm1(40);
    double dp_l0,dp_l1;
    if(prm.p>0.04){
      dp_l0 = -2E-5;
      dp_l1 =  1E-5;
    }
    else{
      dp_l0 = -2E-5;
      dp_l1 =  2E-5;
    }
#pragma omp parallel for
    for (int i = 0; i < 40; i++)
    {
      prm1[i] = prm;
      //prm1[i].dp_l = (i - 20) * .04 / (30. * 111) / 30.;
      prm1[i].dp_l = (dp_l1-dp_l0)*i/40. +dp_l0;
      S[i] = cal_S(prm1[i], buf_specZ, w_specZ, dx, dy, num_ss);
    }

    for (int i = 0; i < 40; i++)
    {
      if (maxS < S[i])
      {
        maxS = S[i];
        prm = prm1[i];
      }
    }
    cerr << "#dp_l max=" << prm.dp_l << endl;
  }
  int num_loop = 0;
  int num_eig = 0;
  if (prm.D > 0)
  {
    double Sinit = S0;

    Eigen::Matrix4f S_tmp;
    Eigen::Vector4f d, iLambda, dprm, W;
    W(0) = 0.06;
    W(1) = M_PI / 2;
    W(2) = M_PI / 2;
    W(3) = .04 / (30. * 111);
    double epsilon = 0;

    prm0 = prm;
    S0 = cal_S(prm0, buf_specZ, w_specZ, dx, dy, num_ss);

    double S1 = S0 / 2;

    for (int i = 0; i < 100; ++i)
    {
      cal_matS(prm0, buf_specZ, w_specZ, dx, dy, dS, ddS, num_ss);
      for (int k = 0; k < 4; k++)
      {
        for (int l = 0; l < 4; l++)
        {
          ddS(k, l) *= (W(l) * W(k));
          S_tmp(k, l) = ddS(k, l);
        }
        dS(k) *= W(k);
        d(k) = dS(k);
      }
      //
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix4f> eigensolver(S_tmp);
      Eigen::Vector4f Lambda = eigensolver.eigenvalues();
      //if(Lambda(2)>0)return(-1);
      Eigen::Matrix4f Q = eigensolver.eigenvectors();
      //cerr <<"#Eigen\n"<< eigensolver.eigenvalues() << endl;
      //cerr <<"#Test\n"<< Q*Lambda.asDiagonal()*Q.transpose() << endl;
      Eigen::Vector4f dd = Q.transpose() * d;

      //cerr << "#"<<prm0.dp_l<<endl;
      num_eig = 0;
      for (int k = 0; k < 4; k++)
      {
        //cerr << dd(k)/Lambda(k) <<" "<<Lambda(k)<<endl;
       // if (abs(dd(k) / Lambda(k))>0.5||Lambda(k) > 0)
       //   iLambda(k) = 0.;
       // else iLambda(k) = dd(k) / Lambda(k);
        if(Lambda(k) > 0){
	  iLambda(k) = -dd(k) / Lambda(k)/2.;
	}
        else{
	  iLambda(k) =  dd(k) / Lambda(k);
	  num_eig++;
        }
      }
      dprm = W.asDiagonal() * Q * iLambda;
      //
      //for(int idx=0;idx<4;idx++) ddS(idx,idx) += ddS(0,0)*5E-3;
      //lu_factorize(ddS,pm);
      //lu_substitute(ddS,pm,dS);//dS => dpo
      //for(int k=0; k<4; k++) dS(k) *= W(k);
      //if(i==0){
      //prm0.dp = dS(0)/prm0.p;////sqrt((dS[0]*dS[0]+dS[1]*dS[1]))/prm0.p;
      //prm0.dD = dS(2);
      //}

      // double r = 1.0;
      bool flag_loop = 0;
      for (double r = 1.0; r > 1E-2; r *= 0.8)
      {
        prm_tmp = prm0;
        prm_tmp.p -= r * dprm(0);     //dS(0);
        prm_tmp.theta -= r * dprm(1); //dS(1);
        prm_tmp.D -= r * dprm(2);     //dS(2);
        prm_tmp.dp_l -= r * dprm(3);  //dS(3);

        S1 = cal_S(prm_tmp, buf_specZ, w_specZ, dx, dy, num_ss);
        epsilon = ((S1 - S0) / S0);
        if (S1 > S0)
        {
          flag_loop = 1;
          break;
        }
      }
      if (flag_loop == 1)
      {
        S0 = S1;
        prm0 = prm_tmp;
      }
      else
      {
        return (0);
      }
      //cout <<"##"<< S0/Sinit<< " "<<epsilon<<" "<<dprm(2)*180./M_PI<<endl;
      //if(epsilon<0) return(-1);
      if (abs(epsilon) < 1E-9 && dprm(2) * 180. / M_PI < 0.1)
      {
        num_loop = i+1;
        break;
      }
      //if(i==3)r=.6;
    }
    if(num_loop == 0 || num_eig !=4 ) return(0);
    double fmax = est_fmax(prm0, buf_specZ, w_specZ, dx, dy, num_ss);
    prm = prm0;
    prm.dS = S0 / Sinit;
    prm.epsilon = epsilon;
    prm.fmax = fmax;
    prm.dp = sqrt(pow(prm.p * cos(prm.theta) - prm_init.p * cos(prm_init.theta), 2) +
                  pow(prm.p * sin(prm.theta) - prm_init.p * sin(prm_init.theta), 2));
    prm.dD = prm.D - prm_init.D;
        return (num_loop);
  }

  return (0);
}

static double cal_S(const PARAM prm, const array3c &buf_specZ, const array2d &w_specZ, const dvector &dx, const dvector &dy, const int num_ss)
{
  //const int num_ss  = buf_specZ.shape()[0];
  const int num_sta = buf_specZ.shape()[1];
  //const int STATION::if1 = buf_specZ.index_bases()[2];
  //const int STATION::if2 = buf_specZ.shape()[2]+buf_specZ.index_bases()[2]-1;
  //cout << num_ss<<" "<<num_sta<<" "<<STATION::if1<<" "<<STATION::if2<<endl;
  cvector phi(STATION::if2 + 1);

  double p = prm.p;
  double ex = cos(prm.theta);
  double ey = sin(prm.theta);

  dvector tau(num_sta, 0);
  for (int i = 0; i < num_sta; i++)
  {
    double eta = (ex * dx[i] + ey * dy[i]) / 6371.;
    double zeta = (-ey * dx[i] + ex * dy[i]) / 6371.;
    double cotDelta = cos(prm.D) / sin(prm.D);
    double dp_l = prm.dp_l;
    double l = (-eta + zeta * zeta * cotDelta / 2. + eta * zeta * zeta * (1. / 6. + cotDelta * cotDelta / 2.)) * 6371.;

    tau[i] = l * (p + dp_l * l / 2);
  }

  double S = 0;
  for (int ibuf = 0; ibuf < num_ss; ibuf++)
  {
    for (int k = STATION::if1; k <= STATION::if2; ++k)
      phi[k] = 0;
    for (int i = 0; i < num_sta; i++)
    {
      double phase = tau[i] * 2. * M_PI * STATION::df;
      double dc = cos(phase);
      double ds = sin(phase);
      double cp0 = cos(phase * STATION::if1);
      double sp0 = sin(phase * STATION::if1);

      for (int k = STATION::if1; k <= STATION::if2; ++k)
      {
        double cp1 = cp0 * dc - sp0 * ds;
        double sp1 = sp0 * dc + cp0 * ds;
        phi[k] += (buf_specZ[ibuf][i][k] * complex<double>(cp0, sp0) * w_specZ[ibuf][i]); //*(double)flag_ss[ist]);
        cp0 = cp1;
        sp0 = sp1;
      }
    }
    for (int k = STATION::if1; k <= STATION::if2; ++k)
      S += real(phi[k] * conj(phi[k]));
  }
  return (S);
}
static double est_fmax(const PARAM prm, const array3c &buf_specZ, const array2d &w_specZ, const dvector &dx, const dvector &dy, const int num_ss)
{
  const int num_sta = buf_specZ.shape()[1];
  //const int STATION::if1 = buf_specZ.index_bases()[2];
  //const int STATION::if2 = buf_specZ.shape()[2]+buf_specZ.index_bases()[2]-1;
  cvector phi(STATION::nfreq);
  dvector psd(STATION::nfreq, 0.);

  double p = prm.p;
  double ex = cos(prm.theta);
  double ey = sin(prm.theta);

  dvector tau(num_sta, 0);
  for (int i = 0; i < num_sta; i++)
  {
    double eta = (ex * dx[i] + ey * dy[i]) / 6371.;
    double zeta = (-ey * dx[i] + ex * dy[i]) / 6371.;
    double cotDelta = cos(prm.D) / sin(prm.D);
    double dp_l = prm.dp_l;
    double l = (-eta + zeta * zeta * cotDelta / 2.) * 6371.;
    tau[i] = l * (p + dp_l * l / 2);
  }

  for (int ibuf = 0; ibuf < num_ss; ibuf++)
  {
    for (int k = STATION::if1; k < STATION::nfreq; ++k)
      phi[k] = 0;
    for (int i = 0; i < num_sta; i++)
    {
      double phase = tau[i] * 2. * M_PI * STATION::df;
      double dc = cos(phase);
      double ds = sin(phase);
      double cp0 = cos(phase * STATION::if1);
      double sp0 = sin(phase * STATION::if1);

      for (int k = STATION::if1; k < STATION::nfreq; ++k)
      {
        double cp1 = cp0 * dc - sp0 * ds;
        double sp1 = sp0 * dc + cp0 * ds;
        phi[k] += (buf_specZ[ibuf][i][k] * complex<double>(cp0, sp0) * w_specZ[ibuf][i]); //*(double)flag_ss[ist]);
        cp0 = cp1;
        sp0 = sp1;
      }
    }
    for (int k = STATION::if1; k < STATION::nfreq; ++k)
      psd[k] += real(phi[k] * conj(phi[k]));
  }
  auto itr = max_element(psd.begin(), psd.end());
  size_t ifmax = distance(psd.begin(), itr);

  return (ifmax * STATION::df);
}

static double cal_matS(const PARAM prm, const array3c &buf_specZ, const array2d &w_specZ, const dvector &dx, const dvector &dy,
                       dvector &dS, dmatrix &ddS, const int num_ss)
{
  //const int num_ss  = buf_specZ.shape()[0];
  const int num_sta = buf_specZ.shape()[1];
  //const int STATION::if1 = buf_specZ.index_bases()[2];
  //const int STATION::if2 = buf_specZ.shape()[2]+buf_specZ.index_bases()[2]-1;

  dvector tau(num_sta);
  array2d dtau(boost::extents[num_sta][4]);
  array3d ddtau(boost::extents[num_sta][4][4]);

  cvector phi(STATION::if2 + 1);
  cmatrix dphi(4, STATION::if2 + 1);
  array3c ddphi(boost::extents[4][4][STATION::if2 + 1]);

  const double p = prm.p;
  const double ex = cos(prm.theta);
  const double ey = sin(prm.theta);
  const double dp_l = prm.dp_l;
  const double cotDelta = cos(prm.D) / sin(prm.D);
  const double sinD = sin(prm.D);

  for (int i = 0; i < num_sta; i++)
  {
    double eta = (ex * dx[i] + ey * dy[i]) / 6371.;
    double zeta = (-ey * dx[i] + ex * dy[i]) / 6371.;
    //double l = (-eta+zeta*zeta*cotDelta/2.) *6371.;
    double l = (-eta + zeta * zeta * cotDelta / 2. + eta * zeta * zeta * (1. / 6. + cotDelta * cotDelta / 2.)) * 6371.;

    //double dl_t = -zeta*(1+eta*cotDelta)*6371.;
    //double dl_D = -pow(zeta/sinD,2)/2*6371.;
    //double ddl_t2 = (eta+(pow(eta,2)-pow(zeta,2))*cotDelta)*6371.;
    //double ddl_tD = zeta*eta/(pow(sinD,2))*6371.;
    //double ddl_D2 = (pow(zeta/sinD,2)*cotDelta)*6371.;
    double dl_t = (-zeta * (1 + eta * cotDelta) + (1 / 6. + cotDelta * cotDelta / 2.) * zeta * (zeta * zeta - 2 * eta * eta)) * 6371.;
    double dl_D = (-pow(zeta / sinD, 2) / 2 - eta * zeta * zeta * cotDelta / (sinD * sinD)) * 6371.;
    double ddl_t2 = ((eta + (pow(eta, 2) - pow(zeta, 2)) * cotDelta) + 1 / 6. * (1 + 3 * cotDelta * cotDelta) * (2 * eta * eta * eta - 7 * eta * zeta * zeta)) * 6371.;
    double ddl_tD = (zeta * eta / (pow(sinD, 2)) + zeta / (sinD * sinD) * cotDelta * (2 * eta * eta - zeta * zeta)) * 6371.;
    double ddl_D2 = (pow(zeta / sinD, 2) * (cotDelta + eta * (3 / (sinD * sinD) - 2))) * 6371.;

    tau[i] = l * (p + dp_l * l / 2);

    //0: p, 1: θ, 2: Δ
    dtau[i][0] = l;
    dtau[i][1] = p * dl_t + dp_l * l * dl_t;
    dtau[i][2] = p * dl_D + dp_l * l * dl_D;
    dtau[i][3] = pow(l, 2) / 2.;

    ddtau[i][0][0] = 0;
    ddtau[i][0][1] = dl_t;
    ddtau[i][0][2] = dl_D;
    ddtau[i][0][3] = 0.;
    ddtau[i][1][0] = ddtau[i][0][1];
    ddtau[i][1][1] = p * ddl_t2 + dp_l * (pow(dl_t, 2) + l * ddl_t2);
    ddtau[i][1][2] = p * ddl_tD + dp_l * (dl_t * dl_D + l * ddl_tD);
    ddtau[i][1][3] = l * dl_t;
    ddtau[i][2][0] = ddtau[i][0][2];
    ddtau[i][2][1] = ddtau[i][1][2];
    ddtau[i][2][2] = p * ddl_D2 + dp_l * (pow(dl_D, 2) + l * ddl_D2);
    ddtau[i][2][3] = l * dl_D;
    ddtau[i][3][0] = 0.;
    ddtau[i][3][1] = ddtau[i][1][3];
    ddtau[i][3][2] = ddtau[i][2][3];
    ddtau[i][3][3] = 0.;
  }
  double S = 0;
  dS.clear();
  ddS.clear();

  for (int ibuf = 0; ibuf < num_ss; ibuf++)
  {
    phi.clear();
    dphi.clear();
    fill_n(ddphi.data(), ddphi.num_elements(), 0.);

    for (int i = 0; i < num_sta; i++)
    {
      double phase = tau[i] * 2. * M_PI * STATION::df;
      double dc = cos(phase);
      double ds = sin(phase);
      double cp0 = cos(phase * STATION::if1);
      double sp0 = sin(phase * STATION::if1);
      ///from python
      for (int k = STATION::if1; k <= STATION::if2; ++k)
      {
        double omega = k * STATION::df * 2. * M_PI;
        double cp1 = cp0 * dc - sp0 * ds;
        double sp1 = sp0 * dc + cp0 * ds;
        complex<double> amp = buf_specZ[ibuf][i][k] * complex<double>(cp0, sp0) * w_specZ[ibuf][i];
        phi[k] += amp;
        for (int m = 0; m < 4; m++)
        {
          dphi(m, k) += (amp * dtau[i][m] * complex<double>(0, 1) * omega);
          for (int n = m; n < 4; n++)
          {
            ddphi[m][n][k] += (-pow(omega, 2) * dtau[i][m] * dtau[i][n] + complex<double>(0, 1) * omega * ddtau[i][m][n]) * amp;
          }
        }
        cp0 = cp1;
        sp0 = sp1;
      }
    }
    for (int k = STATION::if1; k <= STATION::if2; ++k)
      S += real(phi[k] * conj(phi[k]));

    for (int k = STATION::if1; k <= STATION::if2; ++k)
    {
      for (int m = 0; m < 4; m++)
      {
        dS[m] += 2 * real(conj(phi[k]) * dphi(m, k));
        for (int n = m; n < 4; n++)
        {
          ddS(m, n) += 2 * real(conj(dphi(m, k)) * dphi(n, k) + conj(phi[k]) * ddphi[m][n][k]);
        }
      }
    }
  }
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < m; n++)
      ddS(m, n) = ddS(n, m);

  return (S);
}

int eval_deri(vector<STATION> &sta0, double th, double &median0, double &median1)
{
  double deri, integ_cur;
  SPCTRM specE, specN, specZ;

  for (int ist = 0; ist < (int)sta0.size(); ist++)
  { //Calculation of spectra
    sta0[ist].print_spec(specE, specN, specZ);

    tmp_integ0[ist] = specZ.integ[0];
    tmp_integ1[ist] = specZ.integ[1];
  }

  qsort(tmp_integ0, sta0.size(), sizeof(double), (int (*)(const void *, const void *))double_cmp);
  qsort(tmp_integ1, sta0.size(), sizeof(double), (int (*)(const void *, const void *))double_cmp);

  integ_cur = tmp_integ0[sta0.size() / 2];
  median0 = tmp_integ0[sta0.size() / 2];
  median1 = tmp_integ1[sta0.size() / 2];

  if (integ_cur != 0 && integ_pre != 0)
    deri = fabs(log10(integ_cur / integ_pre));
  else
    deri = 100;
  //fprintf(stderr,"###Der %g %g %g\n",deri,integ_old,integ_cur);
  integ_pre = integ_cur;
  if (integ_old == 0 && deri < th)
    integ_old = integ_cur;
  if (deri < th && integ_cur < integ_old * 2 && integ_cur > integ_old / 2)
  {
    integ_old = integ_cur;
    return (1);
  }
  else
    return (0);
}

int read_station_correction(vector<STATION> &sta0)
{
  array3d alpha(boost::extents[60][range2d(1280, 1460)][range2d(300, 460)]);
  dvector dt(sta0.size(), 0.);
  double dep, lon, lat, vs, vp, rho;

  ifstream ifs("./vel_Nishida2008.dat");
  while (ifs)
  {
    ifs >> dep >> lon >> lat >> vs >> vs >> vp >> rho;
    int idp = (int)(dep + .5);
    int ilon = (int)(lon / .1 + .5);
    int ilat = (int)(lat / .1 + .5);
    alpha[idp][ilon][ilat] = vp;
  }
  ifs.close();

  for (int ist = 0; ist < (int)sta0.size(); ist++)
  {
    double dtau = 0;
    int ilon = (int)(sta0[ist].print_lon() / 0.1);
    lat = sta0[ist].print_lat();
    lat = lat - (11.55 / 60.) * sin(2 * lat * M_PI / 180); //geographic => geocentirc
    int ilat = (int)(lat / 0.1);
    double stel = sta0[ist].print_ele() / 1E3;
    if (stel > 0 && alpha[0][ilon][ilat] > 0)
    {
      dtau = stel / alpha[0][ilon][ilat];
      for (int i = 0; i < 60; i++)
        dtau += 1. / alpha[i][ilon][ilat];
    }
    else if (stel <= 0 && alpha[0][ilon][ilat] > 0)
    {
      //stel = -2.1
      int idp0 = (int)(-stel);
      dtau = (1 + stel + idp0) / alpha[idp0][ilon][ilat];
      for (int i = idp0 + 1; i < 60; i++)
        dtau += 1. / alpha[i][ilon][ilat];
    }
    else
      dtau = 9.1;

    dtau -= 9.;
    sta0[ist].set_dtau(dtau);
    //cout << sta0[ist].print_lon()<<" "<<sta0[ist].print_lat()<<" "<<dtau<<endl;
  } //exit(0);
  return (0);
}

int search_max(const array2d &ssZ, dvector &px, dvector &py, const int lmax, double &mad)
{
  const int number = (2 * ipmax + 1) * (2 * ipmax + 1);
  const int imask = 10;
  array2b mask(boost::extents[range2d(-ipmax - imask, ipmax + imask + 1)][range2d(-ipmax - imask, ipmax + imask + 1)]);
  fill_n(mask.data(), mask.num_elements(), true);

  struct ss_data
  {
    double val;
    int ipx, ipy;
  };
  vector<ss_data> ss1(number);
  dvector sz(number);

  int count = 0;
  for (int i = -ipmax; i <= ipmax; i++)
  {
    for (int j = -ipmax; j <= ipmax; j++)
    {
      ss1[count].val = ssZ[i][j];
      ss1[count].ipx = i;
      ss1[count].ipy = j;
      sz[count] = ssZ[i][j];
      count++;
    }
  }
  //Calculation of MAD
  sort(sz.begin(), sz.end());
  double med = sz[number / 2];
  for (int i = 0; i < number; i++)
    sz[i] = fabs(sz[i] - med);
  sort(sz.begin(), sz.end());
  mad = sz[number / 2];

  sort(ss1.begin(), ss1.end(), [](const ss_data &a, const ss_data &b) { return (a.val) > (b.val); });

  count = 0;
  int idx = 0;
  for (int i = 0; i < number; i++)
  {
    if (mask[ss1[i].ipx][ss1[i].ipy] == true)
    {
      if (ss1[i].val < mad * 5)
        break;
      px[idx] = ss1[i].ipx * dp;
      py[idx] = ss1[i].ipy * dp;
      for (int k = ss1[i].ipx - imask; k <= ss1[i].ipx + imask; k++)
      {
        for (int l = ss1[i].ipy - imask; l <= ss1[i].ipy + imask; l++)
        {
          mask[k][l] = false;
        }
      }
      idx++;
    }
    count++;
    if (idx >= lmax)
      break;
  }
  return (idx);
}
