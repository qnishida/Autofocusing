#include <complex>
#include <string>
#include <fftw3.h>
#include <vector>
#include <hdf5.h>
#include <hdf5_hl.h>

//using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;

struct SAC_data{
  float *sgram;
  int npts;
  ptime ts,te;
  double Dt;
  double cmpaz;
  //int pos;
};

struct SPCTRM{
  std::complex<double> *spec;
  //double *amp;
  std::vector<double> integ = std::vector<double>(5); 
  //std::vector<double> integ(5) was wrong because it was interpreted as a definition of a function
  //double integ[5];
};

class STATION{
 public:
  std::string print_net(){return(net);};
  std::string print_sta(){return(sta);};
  double print_dx() const {return(dx);};
  double print_dy() const {return(dy);};
  double print_rad() const {return(rad);};
  double print_lon() const {return(lon);};
  double print_lat() const {return(lat);};
  double print_ele() const {return(elev);};
  double print_dtau() const {return(dtau);};
  int print_spec(SPCTRM &specE2, SPCTRM &specN2, SPCTRM &specZ2){specE2=specE;specN2=specN;specZ2=specZ;return(0);};
  SAC_data print_sacE(){return(sacE);};
  SAC_data print_sacN(){return(sacN);};
  SAC_data print_sacZ(){return(sacZ);};

  int set_loc(double dx2, double dy2, double rad2){dx=dx2;dy=dy2;rad=rad2;return(0);};
  int set_dtau(double correction){dtau=correction;return(0);};
  int set_station(std::string net0, std::string sta0,double lon,double lat, double elev);
  int set_SAC_data(std::string cmp, SAC_data sac_buf);
  int cal_spec(ptime t1);

  int clear_sac();
  int print_sta_num(){return(num_sac);};
  
  STATION();
  STATION(const STATION& obj);
  virtual ~STATION();
  void operator=(const STATION& src);
  static int dt_msec;  // in millisec 
  static int len; // Length in number of segments (2^n)
  static int stride;
  static double df;
  static int nfreq;
  static int if1,if2;
  static int npts;
  static int init_Freq();
  static double *NLNM;
  static double nl_h,nl_v;
  //static int flag_ary;
  static double lon_ary,lat_ary;
 private:
  static double *freq;
  static double *taper;
  double *in,*out;
  fftw_plan plan;
  int cal_spec_1st(ptime t0,SAC_data sac0,SPCTRM &spec0);
  
  std::string net,sta;
  double lon,lat,elev;
  double dx,dy,dtau,rad;//dx, dy are distance from the array center, dtau is the crustal correction and rad is the radius from the center

  int countE,countN,countZ;
  SPCTRM specE,specN,specZ;

  SAC_data sacE,sacN,sacZ;
  int num_sac;
};


int get_station_num(std::vector<STATION> &sta0,std::string sta,std::string net);
void init_station(const std::string h5file, std::vector<STATION> &sta0,const double rad0, const double rad1, hid_t fapl);
int read_h5(std::vector<STATION> &sta0,std::string h5_file, hid_t fapl);
