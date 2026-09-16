//#define  BOOST_FILESYSTEM_VERSION 2
#include <string>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <unordered_map>
#include <cstdlib>
#include <complex>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/algorithm/string.hpp>//split, join, etc.
#include <omp.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <fftw3.h>
#include <math.h>

#include "station_info.h"
#include "util.h"

//using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;

#define ENZ 1

void init_station(const std::string file_hd5,std::vector<STATION> &sta0, hid_t fapl);//fapl: cache

/////Class library
//int STATION::flag_ary=0;
double STATION::lon_ary, STATION::lat_ary;
bool STATION::horizontal_only = false;
int STATION::dt_msec;
int STATION::len;
int STATION::stride;
double STATION::df;
int STATION::nfreq;
int STATION::if1, STATION::if2;
double *STATION::freq;
int STATION::npts;
double *STATION::NLNM;
double STATION::nl_h,STATION::nl_v;
//static
double *STATION::taper;

static inline double window(int i){
double facm=(STATION::len-1.)/2.;
double facp=2./(STATION::len+1);
double win;
win=1.0-pow(((i-1.)-facm)*facp,2.);
return(win);
}

///Constructor destructor
void STATION::operator=(const STATION &src){
  if(this== & src)return;

  net=src.net; sta=src.sta;
  lon=src.lon; lat=src.lat;elev = src.elev;
  countE=src.countE;  countN=src.countN;  countZ=src.countZ;
  dtau = src.dtau;
  sacE.ts=src.sacE.ts; sacN.ts=src.sacN.ts; sacZ.ts=src.sacZ.ts;
  sacE.te=src.sacE.te; sacN.te=src.sacN.te; sacZ.te=src.sacZ.te;
  sacE.Dt=src.sacE.Dt; sacN.Dt=src.sacN.Dt; sacZ.Dt=src.sacZ.Dt;
  sacE.cmpaz=src.sacE.cmpaz; sacN.cmpaz=src.sacN.cmpaz; sacZ.cmpaz=src.sacZ.cmpaz;
  sacE.npts=src.sacE.npts; sacN.npts=src.sacN.npts; sacZ.npts=src.sacZ.npts;

  num_sac =src.num_sac;

  for(int i=0;i<5;i++){
    specE.integ[i]=src.specE.integ[i];
    specN.integ[i]=src.specN.integ[i];
    specZ.integ[i]=src.specZ.integ[i];
  }
}

STATION::STATION(const STATION &src)
{
  
  net=src.net; sta=src.sta;
  lon=src.lon; lat=src.lat; elev = src.elev;
  dtau = src.dtau;
  countE=-1;  countN=-1;  countZ=-1;

  sacE.ts=src.sacE.ts; sacN.ts=src.sacN.ts; sacZ.ts=src.sacZ.ts;
  sacE.te=src.sacE.te; sacN.te=src.sacN.te; sacZ.te=src.sacZ.te;
  sacE.Dt=src.sacE.Dt; sacN.Dt=src.sacN.Dt; sacZ.Dt=src.sacZ.Dt;
  sacE.cmpaz=src.sacE.cmpaz; sacN.cmpaz=src.sacN.cmpaz; sacZ.cmpaz=src.sacZ.cmpaz;
  sacE.npts=src.sacE.npts; sacN.npts=src.sacN.npts; sacZ.npts=src.sacZ.npts;

  num_sac =src.num_sac;

  for(int i=0;i<5;i++){
    specE.integ[i]=src.specE.integ[i];
    specN.integ[i]=src.specE.integ[i];
    specZ.integ[i]=src.specE.integ[i];
  }
  
  specE.spec = new std::complex<double> [STATION::nfreq];
  specN.spec = new std::complex<double> [STATION::nfreq];
  specZ.spec = new std::complex<double> [STATION::nfreq];

//specE.amp = new double [STATION::nfreq];
//specN.amp = new double [STATION::nfreq];
//specZ.amp = new double [STATION::nfreq];

  sacE.sgram = new float [STATION::npts];
  sacN.sgram = new float [STATION::npts];
  sacZ.sgram = new float [STATION::npts];

  in   = new double [STATION::len+2];//
  out  = new double [STATION::len+2];//
  plan = fftw_plan_r2r_1d(STATION::len,in, out,FFTW_R2HC,FFTW_MEASURE);//
}
STATION::STATION() 
{
  countE=-1;
  countN=-1;
  countZ=-1;
  num_sac=0;
  dtau=0;

  specE.spec = new std::complex<double> [STATION::nfreq];
  specN.spec = new std::complex<double> [STATION::nfreq];
  specZ.spec = new std::complex<double> [STATION::nfreq];

//specE.amp = new double [STATION::nfreq];
//specN.amp = new double [STATION::nfreq];
//specZ.amp = new double [STATION::nfreq];
  
  sacE.sgram = new float [STATION::npts];
  sacN.sgram = new float [STATION::npts];
  sacZ.sgram = new float [STATION::npts];

  in   = new double [STATION::len+2];//
  out  = new double [STATION::len+2];//
  plan = fftw_plan_r2r_1d(STATION::len,in, out,FFTW_R2HC,FFTW_MEASURE);
}
STATION::~STATION(){
  delete[] sacE.sgram;
  delete[] sacN.sgram;
  delete[] sacZ.sgram;

  delete[] specE.spec;
  delete[] specN.spec;
  delete[] specZ.spec;

//delete[] specE.amp;
//delete[] specN.amp;
//delete[] specZ.amp;

  delete[] in;
  delete[] out;
  fftw_destroy_plan(plan);
}
int STATION::init_Freq()
{
  STATION::freq = new double[STATION::nfreq];
  for(int i=0;i<nfreq;i++)STATION::freq[i]=i*STATION::df;
    
//STATION::in   = new double [STATION::len+2];
//STATION::out  = new double [STATION::len+2];
  STATION::taper= new double [STATION::len+2];
  STATION::NLNM = new double [STATION::nfreq];

  //STATION::plan = fftw_plan_r2r_1d(STATION::len,in, out,FFTW_R2HC,FFTW_MEASURE);

  double sumw=0;
  for(int i=0;i<STATION::len;i++)  sumw+=pow(window(i),2.0);   /* Welch  */
  for(int i=0;i<STATION::len;i++)
    taper[i]=window(i)/sqrt(sumw/(STATION::dt_msec*1E-3))*M_SQRT2;   /* Welch  */
  for(int i=1;i<STATION::nfreq;++i){
    //std::err << i*STATION::df<<" "<<mk_NLNM(i*STATION::df)<<std::endl;
    NLNM[i]=mk_NLNM(i*STATION::df)/pow(2.*M_PI*i*STATION::df,2);
    }
  NLNM[0] = 0;
  return(0);
}
/////Set data
int STATION::set_station(std::string net_in, std::string sta_in,double lat_in,double lon_in,double stel)
{
  net=net_in;
  sta=sta_in;
  lon=lon_in;
  lat=lat_in;
  elev = stel;
  
  countN=-1;
  countE=-1;
  countZ=-1;
  
  num_sac =0;
  return(0);
}


////Calculation of spec
int STATION::cal_spec(ptime t1){
  ptime t2= t1 + milliseconds(STATION::len*STATION::dt_msec);

  if(num_sac == (horizontal_only ? 2 : 3) &&
     sacE.ts<=t1 && t2<sacE.te && sacN.ts<=t1 && t2<sacN.te &&
     (horizontal_only || (sacZ.ts<=t1 && t2<sacZ.te))){
#if ENZ
    cal_spec_1st(t1,sacE,specE);
    cal_spec_1st(t1,sacN,specN);
#endif
    if (!horizontal_only)
      cal_spec_1st(t1,sacZ,specZ);
    else {
      std::fill_n(specZ.spec, nfreq, std::complex<double>(0., 0.));
      std::fill(specZ.integ.begin(), specZ.integ.end(), 0.);
    }

#if ENZ
    //Correction for azimuth
    double dcs = cos(sacN.cmpaz/180.*M_PI);
    double dsn = sin(sacN.cmpaz/180.*M_PI);
    
    std::complex <double> tmpE,tmpN;
    for(int i=1;i<STATION::nfreq;i++){
      tmpN=dcs*specN.spec[i]-dsn*specE.spec[i];
      tmpE=dsn*specN.spec[i]+dcs*specE.spec[i];
      specN.spec[i]=tmpN;
      specE.spec[i]=tmpE;
    }


    ///Cal. of amp
    /*
    for(int i=1;i<STATION::nfreq;i++){
      specE.amp[i]=sqrt(norm(specE.spec[i]));
      specN.amp[i]=sqrt(norm(specN.spec[i]));
      specZ.amp[i]=sqrt(norm(specZ.spec[i]));
    }*/
#endif
    //cout <<t1 <<" "<< specZ.integ[0] <<" "<<specZ.integ[1]<<" "<<specN.integ[0]<< std::endl;
    return(1);
  }
  return(0);
}

int STATION::cal_spec_1st(ptime t1,SAC_data sac0,SPCTRM &spec0){//double Dt,SPCTRM &spec0){

  int     i0 =  (t1-sac0.ts).total_milliseconds()/STATION::dt_msec;
  double dt2 = ((t1-sac0.ts).total_milliseconds()-i0*STATION::dt_msec)*1E-3;

  double c0=1.,s0=0.,c1,s1;
  double dc=cos(2.*M_PI*STATION::df*(-sac0.Dt+dt2)),ds=sin(2.*M_PI*STATION::df*(-sac0.Dt+dt2));

  for(int i=0;i<STATION::len;i++) in[i] = sac0.sgram[i+i0];
  rtr(in,STATION::len);
  for(int i=0;i<STATION::len;i++) in[i]=in[i]*STATION::taper[i];

  fftw_execute(plan);
  
  for(int i=1;i<STATION::nfreq;i++){
    c1=c0*dc-s0*ds; //Additional theorem of cos
    s1=c0*ds+s0*dc; //Additional theorem of sin
    spec0.spec[i] =    std::complex<double>(out[i]*c1-out[STATION::len-i]*s1,
				      +out[i]*s1+out[STATION::len-i]*c1);

    c0=c1;
    s0=s1;
  }

  spec0.integ[0]=0.;
  int count=0;
  for(int i=(int)(4E-2/STATION::df);i*STATION::df<1E-1;++i){
    spec0.integ[0] += norm(spec0.spec[i])/STATION::NLNM[i];
    count++;
  }
  spec0.integ[0] /= (double)(count);

  spec0.integ[1]=0.;
  count=0;
  for(int i=(int)(1E-1/STATION::df);i*STATION::df<2E-1;++i){
    spec0.integ[1] += norm(spec0.spec[i])/STATION::NLNM[i];
    count++;
  }
  spec0.integ[1] /= (double)(count);

  spec0.integ[2]=0.;
  count=0;
  // QC retains its original upper bound even when analysis needs more FFT bins.
  const int qc_end = std::min(STATION::nfreq, (int)(2.6E-1/STATION::df));
  for(int i=(int)(2E-1/STATION::df);i<qc_end;++i){
    spec0.integ[2] += norm(spec0.spec[i])/STATION::NLNM[i];
    count++;
  }
  spec0.integ[2] /= (double)(count);

  return(0);
}

int STATION::set_SAC_data(std::string cmp, SAC_data sac_buf){
  SAC_data *sac_tmp;
  int npts2 = std::min(STATION::npts,sac_buf.npts);
  if(npts2>STATION::npts/2){
    if     (cmp=="E"){sacE.cmpaz = 90.; sac_tmp = &sacE;}
    else if(cmp=="N"){sacN.cmpaz = 0.; sac_tmp = &sacN;}
    else if(cmp=="U"){sacZ.cmpaz = 0.; sac_tmp = &sacZ;}
    else{
      std::cout << "Error of component name" << std::endl;
      return(-1);
    }
    sac_tmp->cmpaz= sac_buf.cmpaz;
    sac_tmp->npts = sac_buf.npts;
    sac_tmp->ts = sac_buf.ts;
    sac_tmp->te = sac_buf.te;
    sac_tmp->Dt = sac_buf.Dt;
    for(int k=0;k<sac_buf.npts;k++)sac_tmp->sgram[k]=sac_buf.sgram[k];
    if(sac_tmp->npts > 0) num_sac++;
  }
  return(1);
}

int STATION::clear_sac(){
  num_sac=0;
  sacE.npts = sacN.npts = sacZ.npts = 0;
  std::fill_n(specE.spec, nfreq, std::complex<double>(0., 0.));
  std::fill_n(specN.spec, nfreq, std::complex<double>(0., 0.));
  std::fill_n(specZ.spec, nfreq, std::complex<double>(0., 0.));
  std::fill(specE.integ.begin(), specE.integ.end(), 0.);
  std::fill(specN.integ.begin(), specN.integ.end(), 0.);
  std::fill(specZ.integ.begin(), specZ.integ.end(), 0.);
  return(0);
}

///// Global functions
static std::vector<std::string> component_names(hid_t file_id, const std::string &station) {
  auto has = [&](const std::string &component) {
    return H5LTpath_valid(file_id, ("/" + station + "/" + component + "/sgram").c_str(), 1) > 0;
  };
  if (STATION::horizontal_only) {
    if (has("LE") && has("LN")) return {"LE", "LN"};
    if (has("E") && has("N")) return {"E", "N"};
  } else if (has("U") && has("E") && has("N")) {
    return {"U", "E", "N"};
  }
  return {};
}

namespace {
struct H5File {
  hid_t id;
  H5File(const std::string &path, hid_t fapl)
      : id(H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl)) {}
  ~H5File() { if (id >= 0) H5Fclose(id); }
  H5File(const H5File &) = delete;
  H5File &operator=(const H5File &) = delete;
};
struct StationMetadata {
  std::string name;
  std::vector<std::string> components;
  float latitude = 0, longitude = 0, elevation = 0;
};
std::vector<StationMetadata> station_metadata(hid_t file, bool locations) {
  hsize_t count = 0;
  if (H5Gget_num_objs(file, &count) < 0)
    throw std::runtime_error("Cannot enumerate HDF5 stations");
  std::vector<StationMetadata> result;
  result.reserve(count);
  for (hsize_t i = 0; i < count; ++i) {
    if (H5Gget_objtype_by_idx(file, i) != H5G_GROUP) continue;
    char name[100] = {};
    if (H5Gget_objname_by_idx(file, i, name, sizeof(name)) < 0)
      throw std::runtime_error("Cannot read HDF5 station name");
    StationMetadata item;
    item.name = name;
    item.components = component_names(file, item.name);
    if (locations && (!STATION::horizontal_only || !item.components.empty())) {
      const std::string group = "/" + item.name;
      if (H5LTget_attribute_float(file, group.c_str(), "stlo", &item.longitude) < 0 ||
          H5LTget_attribute_float(file, group.c_str(), "stla", &item.latitude) < 0 ||
          H5LTget_attribute_float(file, group.c_str(), "stel", &item.elevation) < 0 ||
          !std::isfinite(item.longitude) || !std::isfinite(item.latitude) ||
          !std::isfinite(item.elevation))
        throw std::runtime_error("Invalid station coordinates: " + item.name);
    }
    result.push_back(std::move(item));
  }
  return result;
}
void initialize_stations(const std::string &h5_file, std::vector<STATION> &sta0,
                         double rad0, double rad1,
                         const std::vector<StationMetadata> &metadata) {
  sta0.reserve(sta0.size() + metadata.size());
  for (const auto &item : metadata) {
    if (STATION::horizontal_only && item.components.empty()) continue;
    sta0.emplace_back();
    sta0.back().set_station("Hi-net", item.name, item.latitude, item.longitude, item.elevation);
  }
  //Selection of stations and definition of station locations in Cartesian coordinate.
  //if(STATION::flag_ary==0){
  double arrlo=0.,arrla=0.;
  int iarr=0;
  //Initial location of the center of the array
  //double center_lon = 133., center_lat = 34.5;
  double center_lon = 138.5, center_lat = 36.5;
  for(std::vector<STATION>::iterator it = sta0.begin();it!=sta0.end();){
    const GeoDistance geo = calc_geodesic(center_lat, center_lon,
                                          it->print_lat(), it->print_lon());
    
    //cout << it->print_lat()<<" "<< it->print_lon()<<" "<< geo.dist_km<<" "<<rad0<<std::endl;
    if(geo.dist_km > rad0) it=sta0.erase(it);//Station selection
    else{
	if(geo.dist_km<rad1){
	  arrlo +=it->print_lon();
	  arrla +=it->print_lat();
	  iarr++;}
	it++;
    }
  }
  if(iarr == 0){
    std::ostringstream oss;
    oss << "No stations remain inside rad1=" << rad1
        << " km while initializing " << h5_file;
    throw std::runtime_error(oss.str());
  }
  STATION::lon_ary = arrlo/iarr;
  STATION::lat_ary = arrla/iarr;

  //cout<<STATION::lon_ary<<" "<<STATION::lat_ary<<std::endl;exit(0);
  for(int i=0;i < (int)sta0.size();i++){//Location of stations in Cartesian coordinate
    const GeoDistance geo = calc_geodesic(STATION::lat_ary, STATION::lon_ary,
                                          sta0[i].print_lat(), sta0[i].print_lon());
    
    double dx = geo.dist_km*sin(geo.az_deg/180.*M_PI);
    double dy = geo.dist_km*cos(geo.az_deg/180.*M_PI);
    sta0[i].set_loc(dx,dy,sqrt(dx*dx+dy*dy));
  }
}//exit(0);

} // namespace

void init_station(const std::string h5_file, std::vector<STATION> &sta0,
                  const double rad0, const double rad1, hid_t fapl) {
  H5File file(h5_file, fapl);
  if (file.id < 0) throw std::runtime_error("Cannot open HDF5: " + h5_file);
  initialize_stations(h5_file, sta0, rad0, rad1, station_metadata(file.id, true));
}

int get_station_num(std::vector<STATION> &sta0,std::string sta,std::string net){
  int i=-1;

  //std::cerr << sta<<" "<<net<<std::endl;
  for(i=0;i<(int)sta0.size();i++){
    //std::cerr << sta0[i].print_sta()<<std::endl;
    //std::cerr << sta0[i].print_sta() <<"#"<<sta<<"#"<< (sta0[i].print_sta()==sta) <<std::endl;
    //std::cerr << sta0[i].print_net() <<"#"<<net<<"#"<< (sta0[i].print_net()==net) <<std::endl;
    if(sta0[i].print_sta()==sta&&sta0[i].print_net()==net) return(i);
  }
  return(-1);
}

static int read_station_waveforms(std::vector<STATION> &sta0, hid_t file_id,
                                  const std::vector<StationMetadata> &metadata,
                                  double &filter_seconds) {
  struct FilterJob { float *data; int count; float cutoff; };
  std::vector<FilterJob> filter_jobs;
  filter_jobs.reserve(sta0.size() * (STATION::horizontal_only ? 2 : 3));
  int number=0;
  for (auto &station : sta0) station.clear_sac();
  std::vector<int> ibuf(STATION::npts);
  std::vector<float> samples(STATION::npts);
  SAC_data sac_buf{};
  sac_buf.sgram = samples.data();
  std::unordered_map<std::string, int> index;
  for (size_t i = 0; i < sta0.size(); ++i)
    if (sta0[i].print_net() == "Hi-net") index.emplace(sta0[i].print_sta(), int(i));

  for (const auto &item : metadata) {
    const auto &stnm = item.name;
    const auto &cmps = item.components;
    if (!cmps.empty()) {
      const auto found = index.find(stnm);
      const int istnm = found == index.end() ? -1 : found->second;
      if(istnm !=-1){
	//components
	for(size_t j = 0;j<cmps.size();j++){
          int year=0,jday=0,hour=0,min=0,sec=0,msec=0,sr=0;
          float a0=0,cmpaz=0;
          const std::string group = "/"+stnm+"/"+cmps[j];
          const std::string time_group = group+"/time";
          if (H5LTget_attribute_float(file_id, group.c_str(), "sensitivity", &a0) < 0 ||
              H5LTget_attribute_float(file_id, group.c_str(), "cmpaz", &cmpaz) < 0 ||
              H5LTget_attribute_int(file_id, group.c_str(), "npts", &sac_buf.npts) < 0 ||
              H5LTget_attribute_int(file_id, group.c_str(), "sr", &sr) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "year", &year) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "jday", &jday) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "hour", &hour) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "min", &min) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "sec", &sec) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "msec", &msec) < 0 ||
              !std::isfinite(a0) || a0 == 0 || !std::isfinite(cmpaz)) {
            fprintf(stderr, "Invalid waveform metadata: %s\n", group.c_str());
            continue;
          }
          if (STATION::horizontal_only) {
            hsize_t unit_dims[1];
            H5T_class_t unit_class;
            size_t unit_size=0;
            if (H5LTget_attribute_info(file_id, group.c_str(), "unit", unit_dims,
                                      &unit_class, &unit_size) < 0 ||
                unit_class != H5T_STRING || unit_size > 32) {
              fprintf(stderr, "Missing/unsupported velocity unit: %s\n", group.c_str());
              continue;
            }
            char unit[33] = {};
            if (H5LTget_attribute_string(file_id, group.c_str(), "unit", unit) < 0 ||
                std::string(unit) != "nm/s") {
              fprintf(stderr, "Expected nm/s in %s\n", group.c_str());
              continue;
            }
          }
	  date d(date(year,1,1)+date_duration(jday-1));
	  sac_buf.ts = ptime(d,time_duration(hour,min,sec))+milliseconds(msec);
	  sac_buf.te = sac_buf.ts+milliseconds(sac_buf.npts*STATION::dt_msec);
	  sac_buf.Dt = 0;
	  sac_buf.cmpaz = cmpaz;
	  
	  int rank = 0;
          hsize_t dims[1] = {0};
          const std::string dataset = "/"+stnm+"/"+cmps[j]+"/sgram";
          if (H5LTget_dataset_ndims(file_id, dataset.c_str(), &rank) < 0 || rank != 1 ||
              H5LTget_dataset_info(file_id, dataset.c_str(), dims, nullptr, nullptr) < 0 ||
              dims[0] != static_cast<hsize_t>(sac_buf.npts)) {
            fprintf(stderr, "Invalid waveform extent: %s\n", dataset.c_str());
            continue;
          }
          if(1000 != sr*STATION::dt_msec || sac_buf.npts <= 0 || sac_buf.npts> STATION::npts){
	    fprintf(stderr,"Sampling rate (%d) or npts (%d) is wrong,\n",sr,sac_buf.npts);
	    sac_buf.npts = 0;
	  }
	  else{
	    if (H5LTread_dataset_int(file_id, dataset.c_str(), ibuf.data()) < 0) continue;
	    for(int k=0;k<sac_buf.npts;k++) sac_buf.sgram[k]=ibuf[k]*a0*1E-9;
            const std::string component = cmps[j].substr(cmps[j].size()-1);
            if (sta0[istnm].set_SAC_data(component, sac_buf) == 1) {
              const SAC_data stored = component == "E" ? sta0[istnm].print_sacE()
                  : component == "N" ? sta0[istnm].print_sacN() : sta0[istnm].print_sacZ();
              if (stored.npts > 0)
                filter_jobs.push_back({stored.sgram, stored.npts, float(3E-2/(sr*1.))});
            }
	  }
	}
        if (STATION::horizontal_only && sta0[istnm].print_sta_num() == 2) {
          const double az_delta = std::fmod(sta0[istnm].print_sacE().cmpaz -
              sta0[istnm].print_sacN().cmpaz + 720., 360.);
          if (std::abs(az_delta - 90.) > 0.01) {
            fprintf(stderr, "Non-orthogonal LE/LN or E/N axes: %s\n", stnm.c_str());
            sta0[istnm].clear_sac();
          }
        }
        if (sta0[istnm].print_sta_num() == (STATION::horizontal_only ? 2 : 3)) ++number;
      }
    }
  }

  const auto filter_start = std::chrono::steady_clock::now();
  // All HDF5 calls and metadata updates finish before worker threads start.
  // Jobs retain the original lengths even if subsequent QC clears a station.
#pragma omp parallel for schedule(static) if(filter_jobs.size() > 1)
  for (int i = 0; i < static_cast<int>(filter_jobs.size()); ++i) {
    const auto &job = filter_jobs[i];
    hp_filt(job.data, job.data, job.count, job.cutoff);
  }
  filter_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - filter_start).count();
  return(number);
}    
  
//init: location of the center of an array
//    : ibuf[86400*sr]?

int read_h5(std::vector<STATION> &sta0, std::string h5_file, hid_t fapl) {
  for (auto &station : sta0) station.clear_sac();
  H5File file(h5_file, fapl);
  if (file.id < 0) {
    fprintf(stderr, "Cannot open file %s\n", h5_file.c_str());
    return -1;
  }
  double filter_seconds = 0;
  return read_station_waveforms(sta0, file.id, station_metadata(file.id, false), filter_seconds);
}

int load_h5(std::vector<STATION> &sta0, const std::string &h5_file,
            double rad0, double rad1, hid_t fapl) {
  const auto start = std::chrono::steady_clock::now();
  H5File file(h5_file, fapl);
  if (file.id < 0) throw std::runtime_error("Cannot open HDF5: " + h5_file);
  const auto metadata = station_metadata(file.id, true);
  initialize_stations(h5_file, sta0, rad0, rad1, metadata);
  const auto initialized = std::chrono::steady_clock::now();
  double filter_seconds = 0;
  const int count = read_station_waveforms(sta0, file.id, metadata, filter_seconds);
  const double init_seconds = std::chrono::duration<double>(initialized - start).count();
  const double read_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - initialized).count();
  if (std::getenv("AUTOFOCUSING_PROFILE"))
    std::cerr << "#LOAD_PROFILE init_s=" << init_seconds
              << " read_decode_copy_s=" << read_seconds - filter_seconds
              << " filter_s=" << filter_seconds
              << " total_s=" << init_seconds + read_seconds << '\n';
  return count;
}
