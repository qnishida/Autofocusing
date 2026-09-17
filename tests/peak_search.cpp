// Parse Boost's checked array access before the production assertion-disable macro.
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <boost/multi_array.hpp>
#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main

using Point = std::pair<int,int>;
static array3d grid() {
  array3d g(boost::extents[3][range2d(-ipmax,ipmax+1)][range2d(-ipmax,ipmax+1)]);
  for(int c=0;c<3;++c) for(int x=-ipmax;x<=ipmax;++x) for(int y=-ipmax;y<=ipmax;++y)
    g[c][x][y]=.001*(x+ipmax)+.00001*(y+ipmax);
  return g;
}
static std::vector<Point> find(const array3d &g,int capacity=4,int component=2) {
  dvector px(capacity),py(capacity);double mad;std::vector<double> q(3);
  int n=search_max(g,px,py,capacity,mad,q,component);
  std::vector<Point> result;
  for(int i=0;i<n;++i) result.emplace_back(std::lround(px[i]/dp),std::lround(py[i]/dp));
  return result;
}
static void check(const char *name,const array3d &g,std::vector<Point> expected,int capacity=4,int component=2) {
  auto actual=find(g,capacity,component);
  if(actual!=expected) {
    std::cerr<<"FAIL "<<name<<" actual:";
    for(auto xy:actual) std::cerr<<' '<<xy.first<<','<<xy.second;
    std::cerr<<'\n';throw std::runtime_error(name);
  }
  std::cout<<"PASS "<<name<<'\n';
}
int main(int argc, char **argv) try {
  dp=.005;ipmax=33;
  if(argc==2 && std::string(argv[1])=="--boundary-only") {
    auto g=grid();g[2][33][0]=1000;
    check("checked boundary access",g,{{33,0}});
    return 0;
  }
  for(int c=0;c<3;++c) {
    auto g=grid();g[c][0][0]=1000;g[c][10][10]=500;
    check("interior isolated peaks",g,{{0,0},{10,10}},4,c);
    check("candidate limit",g,{{0,0}},1,c);
  }
  {
    auto g=grid();g[2][0][0]=1000;g[2][0][1]=1000;g[2][10][10]=500;
    check("rejected plateau does not hide later peak",g,{{10,10}});
  }
  for(auto xy:std::vector<Point>{{33,0},{-33,0},{0,33},{0,-33},{32,1},{-32,1},{1,32},{1,-32}}) {
    auto g=grid();g[2][xy.first][xy.second]=1000;g[2][0][0]=500;
    check("partial neighborhood retains edge seed",g,{xy,{0,0}});
  }
  {
    auto g=grid();g[2][33][33]=2000;g[2][0][0]=500;
    check("outside-circle center excluded",g,{{0,0}});
  }
  {
    auto g=grid();g[2][25][23]=2000;g[2][23][23]=1000;g[2][0][0]=500;
    check("outside-circle grid remains part of neighborhood",g,{{0,0}});
  }
  {
    auto g=grid();g[2][33][0]=.18;
    check("boundary background uses actual reference count",g,{});
    g[2][33][0]=.20;
    check("boundary contrast above unchanged threshold",g,{{33,0}});
  }
  {
    auto g=grid();g[2][0][0]=.14;
    check("interior contrast below threshold",g,{});
    g[2][0][0]=.17;
    check("interior contrast above threshold",g,{{0,0}});
  }
  {
    auto g=grid();g[2][0][0]=1000;g[2][1][0]=900;g[2][10][10]=500;
    check("nearby sidelobe suppressed",g,{{0,0},{10,10}});
  }
  check("no qualifying peak",grid(),{});
  dp=.005;ipmax=80;
  { auto g=grid();g[2][80][0]=1000;check("wide configured search edge",g,{{80,0}}); }
  dp=.0025;ipmax=66;
  { auto g=grid();g[2][-66][0]=1000;check("finer step edge",g,{{-66,0}}); }
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n';return 1; }
