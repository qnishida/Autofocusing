#include "util.h"
#include "station_info.h"
#include <iostream>
#include <iomanip>
int main(int argc,char **argv) {
  if(argc!=2) return 1;
  STATION::dt_msec=500; STATION::len=2048; STATION::stride=928;
  STATION::df=1./1024; STATION::nfreq=int(.26/STATION::df);
  STATION::npts=172800; STATION::if1=int(.1/STATION::df); STATION::if2=int(.25/STATION::df);
  STATION::init_Freq(); std::vector<STATION> sta(1);
  sta[0].set_station("Hi-net","ENU",36.,138.,0.); sta[0].clear_sac();
  if(read_h5(sta,argv[1],H5P_DEFAULT)!=1 || !sta[0].cal_spec(ptime(date(2004,1,1)))) return 2;
  SPCTRM e,n,z; sta[0].print_spec(e,n,z);
  std::cout << std::setprecision(17);
  for(const auto &s:{e,n,z}) {
    for(int k=1;k<STATION::nfreq;++k) std::cout << s.spec[k].real() << ' ' << s.spec[k].imag() << '\n';
    for(int k=0;k<3;++k) std::cout << s.integ[k] << '\n';
  }
}
