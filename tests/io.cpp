// Reuse the established fixture and compare with the frozen pre-I/O loader.
#define AUTOFOCUSING_IO_TEST
#include "horizontal.cpp"
#include "io_reference.h"
#include <cstring>
#include <iomanip>
#include <sys/resource.h>
#ifdef __APPLE__
#include <libproc.h>
#endif

static int load_current(const std::string &file, std::vector<STATION> &sta, hid_t fapl) {
  init_station(file,sta,1800,1800,fapl);
  return read_h5(sta,file,fapl);
}
static int load_reference(const std::string &file, std::vector<STATION> &sta, hid_t fapl) {
  io_reference::init_station(file,sta,1800,1800,fapl);
  return io_reference::read_h5(sta,file,fapl);
}
static void compare(std::vector<STATION> &a,std::vector<STATION> &b) {
  require(a.size()==b.size(),"station count differs");
  for(size_t i=0;i<a.size();++i) {
    require(a[i].print_sta()==b[i].print_sta() && a[i].print_net()==b[i].print_net(),"station order differs");
    require(a[i].print_dx()==b[i].print_dx() && a[i].print_dy()==b[i].print_dy() && a[i].print_rad()==b[i].print_rad(),"station coordinates differ");
    require(a[i].print_sta_num()==b[i].print_sta_num(),"component acceptance differs");
    SAC_data x[]={a[i].print_sacE(),a[i].print_sacN(),a[i].print_sacZ()};
    SAC_data y[]={b[i].print_sacE(),b[i].print_sacN(),b[i].print_sacZ()};
    for(int c=0;c<3;++c) {
      require(x[c].npts==y[c].npts,"waveform length differs");
      if(!x[c].npts) continue;
      require(x[c].ts==y[c].ts && x[c].te==y[c].te && x[c].Dt==y[c].Dt && x[c].cmpaz==y[c].cmpaz,"waveform metadata differs");
      require(std::memcmp(x[c].sgram,y[c].sgram,x[c].npts*sizeof(float))==0,"waveform bits differ");
    }
  }
}
static void equivalence(const std::string &file,bool horizontal) {
  STATION::horizontal_only=horizontal;
  std::vector<STATION> reference;
  int count=load_reference(file,reference,H5P_DEFAULT);
  const double lat=STATION::lat_ary,lon=STATION::lon_ary;
  for(int threads:{1,4,16}) {
    omp_set_num_threads(threads);
    std::vector<STATION> current;
    require(load_current(file,current,H5P_DEFAULT)==count,"loaded count differs");
    require(lat==STATION::lat_ary && lon==STATION::lon_ary,"array center differs");
    compare(reference,current);
    std::cout<<"PASS mode="<<(horizontal?"horizontal":"3c")<<" threads="<<threads<<" stations="<<count<<" full_waveforms=bitwise_equal\n";
  }
}
static long long disk_bytes() {
#ifdef __APPLE__
  rusage_info_v4 r{};
  if(!proc_pid_rusage(getpid(),RUSAGE_INFO_V4,reinterpret_cast<rusage_info_t *>(&r))) return r.ri_diskio_bytesread;
#endif
  return -1;
}
static void benchmark(const std::string &file,bool reference,int threads,int repeats) {
  omp_set_num_threads(threads);STATION::horizontal_only=true;
  hid_t fapl=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_cache(fapl,0,1049,2147483648ULL,.75);
  for(int i=0;i<repeats;++i) {
    std::vector<STATION> stations;
    const auto before=disk_bytes();
    const auto start=std::chrono::steady_clock::now();
    struct rusage r0{},r1{};getrusage(RUSAGE_SELF,&r0);
    int count=reference?load_reference(file,stations,fapl):load_current(file,stations,fapl);
    double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    auto after=disk_bytes();getrusage(RUSAGE_SELF,&r1);
    auto cpu=[](const rusage &r){return r.ru_utime.tv_sec+r.ru_stime.tv_sec+(r.ru_utime.tv_usec+r.ru_stime.tv_usec)*1e-6;};
    require(count>0,"benchmark loaded no stations");
    std::cout<<std::setprecision(17)<<"{\"reference\":"<<(reference?"true":"false")<<",\"threads\":"<<threads<<",\"repeat\":"<<i
      <<",\"stations\":"<<count<<",\"wall_s\":"<<wall<<",\"cpu_s\":"<<cpu(r1)-cpu(r0)<<",\"disk_read_bytes\":"<<(before<0||after<0?-1:after-before)
      <<",\"max_rss_native\":"<<r1.ru_maxrss<<"}\n";
  }
  H5Pclose(fapl);
}
int main(int argc,char **argv) try {
  configure();
  if(argc==6 && std::string(argv[1])=="--benchmark") {
    std::string mode=argv[3];require(mode=="reference"||mode=="current","benchmark mode");
    int threads=std::stoi(argv[4]),repeats=std::stoi(argv[5]);require(threads>0&&repeats>0,"positive benchmark settings");
    benchmark(argv[2],mode=="reference",threads,repeats);
  } else if(argc==3 && std::string(argv[1])=="--real") {
    equivalence(argv[2],true);
  } else {
    require(argc==2,"provide fixture path, --real FILE, or --benchmark FILE reference|current THREADS REPEATS");
    synthetic(argv[1]);
    equivalence(argv[1],true);equivalence(argv[1],false);
    std::vector<STATION> empty;
    require(read_h5(empty,std::string(argv[1])+".missing",H5P_DEFAULT)==-1,"missing input handling");
  }
  return 0;
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
