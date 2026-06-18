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
#include <string>
#include <vector>

//using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;

namespace ublas = boost::numeric::ublas;//Consider atd::vector

#include <complex>
#include <math.h>

const double D2R = M_PI/180.;

typedef boost::multi_array< std::complex<double>, 4> array4c;
typedef array4c::extent_range range4c;

typedef boost::multi_array< double, 3> array3d;
typedef array3d::extent_range range3d;

typedef boost::multi_array< std::complex<double>, 3> array3c;
typedef array3c::extent_range range3c;

typedef boost::multi_array< double, 2> array2d;
typedef array2d::extent_range range2d;

typedef boost::multi_array< std::complex<double>, 2> array2c;
typedef array2c::extent_range range2c;

typedef boost::multi_array< bool, 2> array2b;
typedef array2b::extent_range range2b;

typedef ublas::vector<double> dvector;
typedef ublas::vector<std::complex<double> > cvector;
typedef ublas::matrix<double> dmatrix;
typedef ublas::matrix<std::complex<double> > cmatrix;

struct GeoDistance {
  double dist_km;
  double az_deg;
  double baz_deg;
  double gcarc_deg;
};

GeoDistance calc_geodesic(double evla, double evlo, double stla, double stlo);
double mk_NLNM(double f);
int rtr(double *data,int w);
int CMT(ptime p0,ptime p1);
int init_CMT(int year0, const std::string &catalog_path);
double cos_win(double f1, double f2, double f3, double f4, double f);
int Legendre(const int l,const double x,double &f0,double &f1,double &f2,double &f3);
double Ylm(int l,int m, double theta,double phi);//Real spherical harmonics
inline int idx(int i,int j);
double eval_mid(double theta1,double phi1,double theta2,double phi2
             ,double *theta3, double *phi3);
int eval_pole(double lon1,double lat1,double lon2,double lat2
	      ,double &lon3, double &lat3);

int lp_filt(float *x, float *y, int npts,float nf);
int bp_filt(float *x, float *y, int npts,float nf);
int hp_filt(float *x, float *y, int npts,float nf);
  
double const rho_w   = 1.02E3;
double const alpha_w = 1450.; 
double const rho_c   = 2.6E3;
double const alpha_c = 5800.; 

int get_CMT(std::vector<ptime> &CMT_ptime,std::vector<double> &CMT_slat,std::vector<double> &CMT_slon,std::vector<double> &CMT_sdep, std::vector<double> &CMT_moment);
