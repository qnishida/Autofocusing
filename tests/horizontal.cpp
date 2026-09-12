// Include the driver to exercise the same internal rotation/quality/matrix code.
#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main
#include <cstdlib>

static void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
static void configure() {
  STATION::dt_msec=500; STATION::len=2048; STATION::stride=928;
  STATION::df=1./1024; STATION::nfreq=int(.26/STATION::df);
  STATION::npts=172800; STATION::if1=int(.1/STATION::df);
  STATION::if2=int(.25/STATION::df); STATION::init_Freq();
}
static void attr(hid_t f, const std::string &p, const char *name, int value) {
  require(H5LTset_attribute_int(f,p.c_str(),name,&value,1)>=0,"write integer attribute");
}
static void attr(hid_t f, const std::string &p, const char *name, float value) {
  require(H5LTset_attribute_float(f,p.c_str(),name,&value,1)>=0,"write float attribute");
}
static void group(hid_t f, const std::string &p) {
  hid_t g=H5Gcreate2(f,p.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  require(g>=0,"create group"); H5Gclose(g);
}
static std::vector<int> waveform(int amplitude) {
  std::vector<int> data(STATION::npts);
  for(size_t i=0;i<data.size();++i)
    data[i]=int(amplitude*(sin(2*M_PI*.0625*i*.5)+
        sin(2*M_PI*.125*i*.5)+sin(2*M_PI*.21875*i*.5)));
  return data;
}
static void component(hid_t f,const std::string &p,int amplitude,float az,
                      int sr=2,const char *unit="nm/s") {
  group(f,p); group(f,p+"/time");
  attr(f,p,"sr",sr); attr(f,p,"npts",STATION::npts); attr(f,p,"gap",0);
  attr(f,p,"sensitivity",.1f); attr(f,p,"cmpaz",az);
  H5LTset_attribute_string(f,p.c_str(),"unit",unit);
  const auto t=p+"/time";
  attr(f,t,"year",2004); attr(f,t,"jday",1); attr(f,t,"hour",0);
  attr(f,t,"min",0); attr(f,t,"sec",0); attr(f,t,"msec",0);
  auto data=waveform(amplitude); hsize_t dims[]={data.size()};
  require(H5LTmake_dataset_int(f,(p+"/sgram").c_str(),1,dims,data.data())>=0,"write waveform");
}
static std::vector<STATION> load(const std::string &file,const std::string &name,
                               bool horizontal,int expected) {
  STATION::horizontal_only=horizontal;
  std::vector<STATION> stations(1);
  stations[0].set_station("Hi-net",name,36.,138.,0.);
  stations[0].set_loc(0,0,0);
  require(read_h5(stations,file,H5P_DEFAULT)==expected,"unexpected loaded station count");
  return stations;
}
static void synthetic(const std::string &file) {
  require((fs::path("a") / "b").string()=="a/b", "Boost filesystem ABI mismatch; use a consistent C++ toolchain");
  hid_t f=H5Fcreate(file.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  for(const auto &name : {"TILT","ENU","MISSING","SKEW","RATE","UNIT"}) {
    const std::string p=std::string("/")+name;
    group(f,p); attr(f,p,"stlo",138.f); attr(f,p,"stla",36.f); attr(f,p,"stel",0.f);
  }
  component(f,"/TILT/LE",2000,120); component(f,"/TILT/LN",1000,30);
  component(f,"/ENU/E",2000,120); component(f,"/ENU/N",1000,30); component(f,"/ENU/U",1500,0);
  component(f,"/MISSING/LE",2000,90);
  component(f,"/SKEW/LE",2000,100); component(f,"/SKEW/LN",1000,0);
  component(f,"/RATE/LE",2000,90,4); component(f,"/RATE/LN",1000,0);
  component(f,"/UNIT/LE",2000,90,2,"rad"); component(f,"/UNIT/LN",1000,0);
  H5Fclose(f);
  auto tilt=load(file,"TILT",true,1);
  const ptime t(date(2004,1,1));
  require(tilt[0].cal_spec(t)==1,"horizontal spectrum requires no U");
  SPCTRM e,n,z; tilt[0].print_spec(e,n,z);
  require(std::abs(e.spec[128])>0,"nonzero horizontal spectrum");
  for(int k=1;k<STATION::nfreq;++k) require(z.spec[k]==std::complex<double>(0,0),"vertical spectrum must be zero internally");
  std::vector<float> scaled(STATION::npts); auto raw=waveform(2000);
  for(size_t k=0;k<scaled.size();++k) scaled[k]=raw[k]*.1f*1E-9;
  hp_filt(scaled.data(),scaled.data(),scaled.size(),.03/2);
  auto sac=tilt[0].print_sacE();
  require(std::abs(sac.sgram[500]-scaled[500])<1E-15,"counts to m/s and filter");
  // Independent zero-azimuth reference, then analytic 30-degree rotation.
  STATION ref; ref.clear_sac(); sac.cmpaz=90; ref.set_SAC_data("E",sac);
  sac=tilt[0].print_sacN(); sac.cmpaz=0; ref.set_SAC_data("N",sac);
  require(ref.cal_spec(t)==1,"reference spectrum"); SPCTRM re,rn,rz; ref.print_spec(re,rn,rz);
  for(int k=STATION::if1;k<=STATION::if2;++k) {
    require(std::abs(e.spec[k]-(cos(M_PI/6)*re.spec[k]+.5*rn.spec[k]))<1E-15,"east azimuth rotation");
    require(std::abs(n.spec[k]-(cos(M_PI/6)*rn.spec[k]-.5*re.spec[k]))<1E-15,"north azimuth rotation");
  }
  auto en=load(file,"ENU",true,1); require(en[0].cal_spec(t)==1,"E/N fallback");
  SPCTRM ee,nn,zz; en[0].print_spec(ee,nn,zz);
  for(int k=1;k<STATION::nfreq;++k) require(e.spec[k]==ee.spec[k] && n.spec[k]==nn.spec[k],"LE/LN and E/N equivalence");
  auto three=load(file,"ENU",false,1); require(three[0].cal_spec(t)==1,"3c spectrum");
  three[0].print_spec(ee,nn,zz);
  for(int k=1;k<STATION::nfreq;++k) require(e.spec[k]==ee.spec[k] && n.spec[k]==nn.spec[k],"3c and horizontal spectra equivalence");
  require(std::abs(zz.spec[128])>0,"3c vertical preserved");
  load(file,"TILT",false,0);
  for(const auto &name : {"MISSING","SKEW","RATE","UNIT"}) {
    auto bad=load(file,name,true,0); require(bad[0].cal_spec(t)==0,"reject incomplete/invalid station");
  }
  STATION::horizontal_only=true;
  integ_pre=integ_old=0; count_gap=0; double m0,m1;
  require(eval_deri(tilt,1.12,m0,m1)==0,"first quality window has no history");
  require(eval_deri(tilt,1.12,m0,m1)==1 && m0>0 && m1>0,"horizontal stability quality without U");
  array4c spec(boost::extents[5][1][2][range4c(STATION::if1,STATION::nfreq)]);
  array3d weights(boost::extents[5][1][2]);
  std::fill_n(spec.data(),spec.num_elements(),std::complex<double>(0,0));
  std::fill_n(weights.data(),weights.num_elements(),0.);
  for(int i=0;i<2;++i) {
    weights[3][0][i]=weights[4][0][i]=1;
    for(int k=STATION::if1;k<STATION::nfreq;++k) spec[3][0][i][k]=1.;
  }
  PARAM prm{}; prm.Δ=M_PI/2.; dvector dx(2,0),dy(2,0);
  auto matrix=cal_S_matrix(prm,spec,weights,dx,dy,1);
  require(std::abs(matrix(0,0).real()-1)<1E-12 && matrix(1,1)==std::complex<double>(0,0),"R/T beam power");
  require(std::isnan(matrix(2,2).real()) && std::isnan(matrix(0,2).imag()),"missing U matrix entries");
  std::vector<STATION> array(2);
  for(auto &station:array) station.set_station("Hi-net","TEST",36,138,0);
  STATION::lat_ary=36; STATION::lon_ary=138;
  prm.θ=M_PI/2; prm.Δ=.2;
  for(int i=0;i<2;++i) {
    weights[0][0][i]=weights[1][0][i]=1;
    for(int k=STATION::if1;k<=STATION::if2;++k) {
      spec[0][0][i][k]=1; spec[1][0][i][k]=2;
    }
  }
  rotate_EN_RT(array,spec,weights,prm);
  require(std::abs(spec[3][0][0][128]-std::complex<double>(2,0))<1e-12 &&
          std::abs(spec[4][0][0][128]-std::complex<double>(1,0))<1e-12,
          "geographic E/N to R/T rotation");
  std::cout << "PASS: loading, units, rotation, missing/invalid components, 3c equivalence, horizontal quality and beam matrix\n";
}
static void real_data(const std::string &file) {
  STATION::horizontal_only=true; std::vector<STATION> sta;
  init_station(file,sta,1800,1800,H5P_DEFAULT);
  int loaded=read_h5(sta,file,H5P_DEFAULT), valid=0, quality=0;
  for(auto &s:sta) {
    if(s.print_sta_num()!=2 || !s.cal_spec(s.print_sacE().ts)) continue;
    ++valid; SPCTRM e,n,z; s.print_spec(e,n,z);
    bool good=true;
    for(int band=0;band<3;++band) {
      double upper=band==0?4E6:4E5;
      good &= e.integ[band]>1 && e.integ[band]<upper && n.integ[band]>1 && n.integ[band]<upper;
    }
    quality+=good;
  }
  require(loaded>0 && valid>0,"real HDF5 load/spectrum");
  std::cout << file << " selected=" << sta.size() << " loaded=" << loaded
            << " first-window-spectra=" << valid << " amplitude-QC=" << quality << '\n';
}
#ifndef AUTOFOCUSING_IO_TEST
int main(int argc,char **argv) {
  try {
    configure();
    if(argc==3 && std::string(argv[1])=="--real") real_data(argv[2]);
    else { require(argc==2,"provide temporary fixture path or --real HDF5"); synthetic(argv[1]); }
  } catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}

#endif
