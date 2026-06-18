#include <iostream>
#include <fstream>
#include <vector>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

//#include <fftw3.h>
#include "util.h"

#include <GeographicLib/Geodesic.hpp>

#include <math.h>

#define MAX_FILT    100
int rtr2(float *data,int w,int width);
int buthip(float *h,int *m,float *gn, int *n,float fp,float fs,float ap,float as);
int butpas(float *h,int *m,float *gn, int *n,float fl,float fh,float fs,float ap,float as);
int butlop(float *h,int *m,float *gn, int *n,float fp,float fs,float ap,float as);
int tandem(float *x,float *y,int n,float *h,int m,int nml,float *uv);
int recfil(float *x,float *y,int n,float *h,int nml,float *uv);


using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace GeographicLib;

static const Geodesic& geod = Geodesic::WGS84();

static double normalize_azimuth(double azimuth);
static double NLNM(double a1, double a2,double t);
static vector<ptime> ts,te;
static vector<double> CMT_lat2,CMT_lon2,CMT_dep2,CMT_moment2;

static double normalize_azimuth(double azimuth)
{
  azimuth = fmod(azimuth, 360.);
  if(azimuth < 0.)
    azimuth += 360.;
  return(azimuth);
}

double mk_NLNM(double f)
{
  double t=1/f,psd=0;

  if(t>0.1&&t<=0.17)
    {
      psd=NLNM(-162.36,5.64,t);
    }                               
  if(t>0.17&&t<=0.40)
    {
      psd=NLNM(-166.7,0.00,t);
    }                               
  if(t>0.40&&t<=0.80)
    {
      psd=NLNM(-170,-8.30,t);
    }                               
  if(t>0.80&&t<=1.24)
    {
      psd=NLNM(-166.4,28.9,t);
    }                               
  if(t>1.24&&t<=2.40)
    {
      psd=NLNM(-168.6,52.48,t);
    }                               
  if(t>2.40&&t<=4.30)
    {
      psd=NLNM(-159.98,29.81,t);
    }                               
  if(t>4.30&&t<=5.00)
    {
      psd=NLNM(-141.1,0.,t);
    }                               
  if(t>5.00&&t<=6.00)
    {
      psd=NLNM(-71.36,-99.77,t);
    }                               
  if(t>6.0&&t<=10)
    {
      psd=NLNM(-97.26,-66.49,t);
    }                               
  if(t>10&&t<=12)
    {
      psd=NLNM(-132.18,-31.57,t);
    }                               
  if(t>12&&t<=15.6)
    {
      psd=NLNM(-205.27,36.16,t);
    }                               
  if(t>15.6&&t<=21.90)
    {
      psd=NLNM(-37.65,-104.33,t);
    }                               
  if(t>21.90&&t<=31.6)
    {
      psd=NLNM(-114.37,-47.1,t);
    }                               
  if(t>31.6&&t<=45)
    {
      psd=NLNM(-160.58,-16.28,t);
    }                               
  if(t>45.&&t<=70)
    {
      psd=NLNM(-187.5,0.00,t);
    }                               
  if(t>70.&&t<=101)
    {
      psd=NLNM(-216.47,15.7,t);
    }                               
  if(t>101&&t<=154)
    {
      psd=NLNM(-185,0.0,t);
    }                               
  if(t>154&&t<=328)
    {
      psd=NLNM(-168.34,-7.61,t);
    }                               
  if(t>328&&t<=600)
    {
      psd=NLNM(-217.43,11.9,t);
    }                               
  if(t>600.&&t<=1E4)
    {
      psd=NLNM(-258.28,26.6,t);
    }                               
  if(t>1E4&&t<=1E6)
    {
      psd=NLNM(-346.88,48.75,t);
    }         
  return(psd);
}

static double NLNM(double a1, double a2,double t)
{
  double tmp;
  tmp=a1+a2*log10(t);
  tmp=pow(10,(tmp/10.));
  return(tmp);
}

////
int rtr(double *data,int w)
{
  double det,a,b;
  int i;
  double x1,x2,xy,y1;
  x1=0;x2=0;xy=0;y1=0;
  for (i=0;i<w;i++)
    {
      x1+=i;
      x2+=(i*i);
      xy+=(i*data[i]);
      y1+=(data[i]);
    }
  det=w*x2-x1*x1;
  a=(w*xy-x1*y1)/det;
  b=(-x1*xy+x2*y1)/det;
              
  for (i=0;i<w;i++)
    {
      data[i]-=(a*i+b);
    }

#if 0
  double sum=0;
  for (i=0;i<w;i++)
    sum+=data[i];
  for (i=0;i<w;i++)
    data[i]-=(sum/(double)(w));
#endif

  return(0);
}

int init_CMT(int year0, const std::string &catalog_path){
  int year,jday,hour,min;
  double sec,moment,slon,slat,sdep,arc,s12,moment_mod;
  string line;
  ifstream ifs(catalog_path.c_str());
  if(!ifs){
    cerr << "Cannot open CMT catalog: " << catalog_path << endl;
    return(-1);
  }
  ts.clear();
  te.clear();
  CMT_lat2.clear();
  CMT_lon2.clear();
  CMT_dep2.clear();
  CMT_moment2.clear();
  while(getline(ifs,line)){
    istringstream iss(line);
    iss >> year>>jday>>hour>>min>>sec>>slon>>slat>>sdep>>moment;
    if(!iss)
      continue;
    year+=1900;
    arc = geod.Inverse(36.8246,137.604,slat,slon,s12);//
    if(arc > 10||arc>170) moment_mod = moment / sin(arc*M_PI/180.)*sin(60*M_PI/180.);
    else moment_mod = moment/sin(10*M_PI/180.)*sin(60*M_PI/180.); 
    if(year0<=year&&moment_mod>1E25) //1E25
      {
	      date d(date(year,1,1)+date_duration(jday-1));
	      ptime t0 = ptime(d,time_duration(hour+9,min,sec+arc*6.0));//Convert from UT to JST
        //cout <<t0 <<" "<< arc <<" "<< arc*6.0/60<<endl;

	      ts.push_back(t0);
	      CMT_lat2.push_back(slat);
	      CMT_lon2.push_back(slon);
	      CMT_dep2.push_back(sdep);
	      CMT_moment2.push_back(moment);
	      
	      //double DT = 360*(log(moment*1E-23)/(2*M_PI*5E-3));//360
	      double DT = 1000*(log(moment*1E-23)/(2*M_PI*1E-1));
	      ptime t1 = t0 + seconds((int)(DT));
	      te.push_back(t1);
      }
  }
  //exit(0);
  if(ts.empty()){
    cerr << "No CMT events loaded from catalog: " << catalog_path << endl;
    return(-1);
  }
  return(0);
}

int CMT(ptime t0,ptime t1){
  int count=0;
  for(int i=0;i< (int)ts.size();++i)
    if(t0< te[i] && t1>= ts[i])count++;
  return(count);
}

int get_CMT(vector<ptime> &CMT_ptime,vector<double> &CMT_slat,vector<double> &CMT_slon,vector<double> &CMT_sdep,vector<double> &CMT_moment){
  CMT_ptime.resize(ts.size());
  copy(ts.begin(),ts.end(),CMT_ptime.begin());

  CMT_slat.resize(CMT_lat2.size());
  copy(CMT_lat2.begin(),CMT_lat2.end(),CMT_slat.begin());

  CMT_slon.resize(CMT_lon2.size());
  copy(CMT_lon2.begin(),CMT_lon2.end(),CMT_slon.begin());

  CMT_sdep.resize(CMT_dep2.size());
  copy(CMT_dep2.begin(),CMT_dep2.end(),CMT_sdep.begin());

  CMT_moment.resize(CMT_moment2.size());
  copy(CMT_moment2.begin(),CMT_moment2.end(),CMT_moment.begin());
  
  return(1);
}


GeoDistance calc_geodesic(double evla, double evlo, double stla, double stlo)
{
  double s12, azi1, azi2;
  const double a12 = geod.Inverse(evla, evlo, stla, stlo, s12, azi1, azi2);

  GeoDistance distance;
  distance.dist_km = s12 / 1000.;
  distance.az_deg = normalize_azimuth(azi1);
  distance.baz_deg = normalize_azimuth(azi2 + 180.);
  distance.gcarc_deg = a12;
  return(distance);
}

double cos_win(double f1, double f2, double f3, double f4, double f)
{
  double win;
  if(f<f1)
    {
      win=0;
    }
  else if(f>=f1&&f<f2)
    {
      win=-0.5*cos(M_PI*(f-f1)/(f2-f1))+.5;
    }
  else if(f>=f2&&f<f3)
    {
      win=1;
    }
  else if(f>=f3&&f<f4)
    {
      win= 0.5*cos(M_PI*(f-f3)/(f4-f3))+0.5;
    }
  else
    {
      win=0;
    }
  return(win);
}


double Ylm(int l,int m, double theta,double phi){//Real spherical harmonics
  namespace sf = boost::math;
  double x;
  if(m<0 )     x = M_SQRT2* sf::spherical_harmonic_r(l, -m, theta,phi);
  else if(m==0)x =          sf::spherical_harmonic_r(l,  0, theta,phi);
  else         x = M_SQRT2* sf::spherical_harmonic_i(l,  m, theta,phi);
  return(x);
}

int Legendre(const int l,const double x,double &f0,double &f1,double &f2,double &f3){
  namespace sf = boost::math;
  double P0,P1;
  if(l==0){
    f0=1.E0;
    f1=0.;
    f2=0.;
    f3=0.;
  }
  else if(x!=0&&x<M_PI-1E-10){
    double cs=cos(x);
    P0  = sf::legendre_p(l,0,cs);
    P1  = sf::legendre_p(l,1,cs);
    
    f0 = P0;
    f1 = (-cs/sin(x)*P1-(l+1)*l*P0)/(l*(l+1.));
    f2 = P1/sqrt(l*(l+1.));
    f3 = P1/sin(x)/(l*(l+1.));
    //cerr <<l<<" "<< f1 << sin(x)<<" "<<sin(x)*P1<<endl;
  }
  else if(x==0)
    {
      f0= 1.E0;
      f1=-.5E0;
      f2= 0.;
      f3=-.5E0;
    }
  else if(x>=M_PI-1E-10)
    {
      if(l%2==0){
	f0= 1.E0;
	f1=-5E-1;
	f3= 5E-1;
      }
      else{
	f0=-1.0E0;
	f1= 5E-1;
	f3=-5E-1;
      }
      //cerr <<"#"<<l<<" "<< f1 <<endl;
      f2= 0.;
    }

  return(0);
}

inline int idx(int i,int j){
  return((i*(i+1))/2+j+1);
}
double eval_mid(double theta1,double phi1,double theta2,double phi2
             ,double *theta3, double *phi3)
{
  double x1,y1,z1;
  double x2,y2,z2;
  double x3,y3,z3;
  double xt,yt,zt,xp,yp,zp,azimuth;

  x1=sin(theta1)*cos(phi1);
  y1=sin(theta1)*sin(phi1);
  z1=cos(theta1);
  
  x2=sin(theta2)*cos(phi2);
  y2=sin(theta2)*sin(phi2);
  z2=cos(theta2);

  x3=(x1+x2)/2.;
  y3=(y1+y2)/2.;
  z3=(z1+z2)/2.;

  *theta3=acos(z3/sqrt(x3*x3+y3*y3+z3*z3));
  *phi3=atan2(y3,x3);

  xt=cos(*theta3)*cos(*phi3);
  yt=cos(*theta3)*sin(*phi3);
  zt=-sin(*theta3);

  xp=-sin(*phi3);
  yp= cos(*phi3);
  zp= 0;

  azimuth=atan2( xt*(x2-x1)+yt*(y2-y1)+zt*(z2-z1)
                ,xp*(x2-x1)+yp*(y2-y1)+zp*(z2-z1));

  *theta3/= D2R;
  *phi3  /= D2R;
  azimuth/= D2R; 

  return(azimuth);
}


int eval_pole(double lon1,double lat1,double lon2,double lat2
	      ,double &lon3, double &lat3)
{
  double x1,y1,z1;
  double x2,y2,z2;
  double x3,y3,z3,nrm;
  double theta1,phi1,theta2,phi2,theta3,phi3;

  theta1= M_PI/2.-(lat1-11.55/60.*sin(2.*lat1*M_PI/180.))/180.*M_PI;
  theta2= M_PI/2.-(lat2-11.55/60.*sin(2.*lat2*M_PI/180.))/180.*M_PI;
  phi1=lon1*D2R;
  phi2=lon2*D2R;

  x1=sin(theta1)*cos(phi1);
  y1=sin(theta1)*sin(phi1);
  z1=cos(theta1);
  
  x2=sin(theta2)*cos(phi2);
  y2=sin(theta2)*sin(phi2);
  z2=cos(theta2);

  x3=y1*z2-z1*y2;
  y3=z1*x2-x1*z2;
  z3=x1*y2-y1*x2;
  
  nrm=sqrt(x3*x3+y3*y3+z3*z3);
  if(nrm==0)
    {
      return(0);
    }

  theta3=acos(z3/nrm);
  phi3=atan2(y3,x3);

  lat3 = (M_PI/2.-theta3)/D2R;// /M_PI*180.;
  lat3= lat3+11.55/60.*sin(2.*lat3*D2R);
  lon3= phi3/D2R;
  return(1);
}


int lp_filt(float *x, float *y, int npts,float nf)
{
    float h[1000],gn,uv[MAX_FILT*4];
    int m,n,i;

    rtr2(x,npts,int(10./nf));
    for(i=0;i<MAX_FILT*4;i++)
      uv[i]=0;

    butlop(h,&m,&gn,&n,nf*.8,nf*1.2,.5,5.);
    tandem(x,y,npts,h,m,1,uv);

    //cout <<"#"<<m<<" "<<gn<<" "<<nf<<endl;
    //for(i=0;i<npts;i++)	cout<<i*10.<<" "<<x[i]<< " "<<y[i]*gn <<endl;
    //exit(0);
    //printf("# %g\n",gn);
    for(i=0;i<npts;i++)	y[i]=y[i]*gn;
    return(0);
}

int hp_filt(float *x, float *y, int npts,float nf)
{
    float h[1000],gn,uv[MAX_FILT*4];
    int m,n,i;

    rtr2(x,npts,int(10./nf));
    for(i=0;i<MAX_FILT*4;i++)
      uv[i]=0;

    buthip(h,&m,&gn,&n,nf*.8,nf*1.2,.5,5.);
    tandem(x,y,npts,h,m,1,uv);

    //cout <<"#"<<m<<" "<<gn<<" "<<nf<<endl;
    //for(i=0;i<npts;i++)	cout<<i*10.<<" "<<x[i]<< " "<<y[i]*gn <<endl;
    //exit(0);
    //printf("# %g\n",gn);
    for(i=0;i<npts;i++)	y[i]=y[i]*gn;
    return(0);
}

int bp_filt(float *x, float *y, int npts,float nf)
{
    float h[1000],gn,uv[MAX_FILT*4];
    int m,n,i;

    rtr2(x,npts,(int)(10./nf));
    for(i=0;i<MAX_FILT*4;i++)   uv[i]=0;

    butpas(h,&m,&gn,&n,nf*.8,nf*2.,nf*3.,.5,5.);
    tandem(x,y,npts,h,m,1,uv);

    //cout <<"#"<<m<<" "<<gn<<" "<<nf<<endl;
    //for(i=0;i<npts;i++)	cout<<i*10.<<" "<<x[i]<< " "<<y[i]*gn <<endl;
    //exit(0);
    for(i=0;i<npts;i++)  y[i]=y[i]*gn;
    return(0);
}


/*
+   BUTTERWORTH BAND PASS FILTER COEFFICIENT
+
+   ARGUMENTS
+   H : FILTER COEFFICIENTS
+   M : ORDER OF FILTER
+   GN  : GAIN FACTOR
+   N : ORDER OF BUTTERWORTH FUNCTION
+   FL  : LOW  FREQUENCY CUT-OFF  (NON-DIMENSIONAL)
+   FH  : HIGH FREQUENCY CUT-OFF
+   FS  : STOP BAND FREQUENCY
+   AP  : MAX. ATTENUATION IN PASS BAND
+   AS  : MIN. ATTENUATION IN STOP BAND
+
+   M. SAITO  (7/I/76)
*/
int butpas(float *h,int *m,float *gn, int *n,float fl,float fh,float fs,float ap,float as){
  float wl,wh,ws,clh,op,ww,ts,os,pa,sa,cc,c,dp,g,fj,rr,tt,
    re,ri,a,wpc,wmc;
  int k,l,j,i;
  struct {
    float r;
    float c;
  } oj,aa,cq,r[2];
  if(fabs(fl)<fabs(fh)) wl=fabs(fl)*M_PI;
  else wl=fabs(fh)*M_PI;
  if(fabs(fl)>fabs(fh)) wh=fabs(fl)*M_PI;
  else wh=fabs(fh)*M_PI;
  ws=fabs(fs)*M_PI;
  if(wl==0.0 || wl==wh || wh>=M_PI_2 || ws==0.0 || ws>=M_PI_2 ||(ws-wl)*(ws-wh)<=0.0){
      fprintf(stderr,
	      "? (butpas) invalid input : fl=%14.6e fh=%14.6e fs=%14.6e ?\n",
	      fl,fh,fs);
      *m=0;
      *gn=1.0;
      return 1;
   }
   /****  DETERMINE N & C */
   clh=1.0/(cos(wl)*cos(wh));
   op=sin(wh-wl)*clh;
   ww=tan(wl)*tan(wh);
   ts=tan(ws);
   os=fabs(ts-ww/ts);
   if(fabs(ap)<fabs(as)) pa=fabs(ap);
   else pa=fabs(as);
   if(fabs(ap)>fabs(as)) sa=fabs(ap);
   else sa=fabs(as);
   if(pa==0.0) pa=0.5;
   if(sa==0.0) sa=5.0;
   if((*n=(int)(fabs(log(pa/sa)/log(op/os))+0.5))<2) *n=2;
   cc=exp(log(pa*sa)/(float)(*n))/(op*os);
   c=sqrt(cc);
   ww=ww*cc;
   
   dp=M_PI_2/(float)(*n);
   k=(*n)/2;
   *m=k*2;
   l=0;
   g=fj=1.0;
   
   for(j=0;j<k;j++){
      oj.r=cos(dp*fj)*0.5;
      oj.c=sin(dp*fj)*0.5;
      fj=fj+2.0;
      aa.r=oj.r*oj.r-oj.c*oj.c+ww;
      aa.c=2.0*oj.r*oj.c;
      rr=sqrt(aa.r*aa.r+aa.c*aa.c);
      tt=atan(aa.c/aa.r);
      cq.r=sqrt(rr)*cos(tt/2.0);
      cq.c=sqrt(rr)*sin(tt/2.0);
      r[0].r=oj.r+cq.r;
      r[0].c=oj.c+cq.c;
      r[1].r=oj.r-cq.r;
      r[1].c=oj.c-cq.c;
      g=g*cc;
      
      for(i=0;i<2;i++){
	 re=r[i].r*r[i].r;
	 ri=r[i].c;
	 a=1.0/((c+ri)*(c+ri)+re);
	 g=g*a;
	 h[l  ]=0.0;
	 h[l+1]=(-1.0);
	 h[l+2]=2.0*((ri-c)*(ri+c)+re)*a;
	 h[l+3]=((ri-c)*(ri-c)+re)*a;
	 l=l+4;
      }
   }
   /****  EXIT */
   *gn=g;
   if(*n==(*m)) return 0;
   /****  FOR ODD N */
   *m=(*m)+1;
   wpc=  cc *cos(wh-wl)*clh;
   wmc=(-cc)*cos(wh+wl)*clh;
   a=1.0/(wpc+c);
   *gn=g*c*a;
   h[l  ]=0.0;
   h[l+1]=(-1.0);
   h[l+2]=2.0*wmc*a;
   h[l+3]=(wpc-c)*a;
   return 0;
}

/*
+   BUTTERWORTH LOW PASS FILTER COEFFICIENT
+
+   ARGUMENTS
+   H : FILTER COEFFICIENTS
+   M : ORDER OF FILTER  (M=(N+1)/2)
+   GN  : GAIN FACTOR
+   N : ORDER OF BUTTERWORTH FUNCTION
+   FP  : PASS BAND FREQUENCY  (NON-DIMENSIONAL)
+   FS  : STOP BAND FREQUENCY
+   AP  : MAX. ATTENUATION IN PASS BAND
+   AS  : MIN. ATTENUATION IN STOP BAND
+
+   M. SAITO  (17/XII/75)
*/
int butlop(float *h,int *m,float *gn, int *n,float fp,float fs,float ap,float as)
{
   float wp,ws,tp,ts,pa,sa,cc,c,dp,g,fj,c2,sj,tj,a;
   int k,j;
   if(fabs(fp)<fabs(fs)) wp=fabs(fp)*M_PI;
   else wp=fabs(fs)*M_PI;
   if(fabs(fp)>fabs(fs)) ws=fabs(fp)*M_PI;
   else ws=fabs(fs)*M_PI;
   if(wp==0.0 || wp==ws || ws>=M_PI_2){
       fprintf(stderr,"? (butlop) invalid input : fp=%14.6e fs=%14.6e ?\n",
	       fp,fs);
       return 1;
   }
   /****  DETERMINE N & C */
   tp=tan(wp);
   ts=tan(ws);
   if(fabs(ap)<fabs(as)) pa=fabs(ap);
   else pa=fabs(as);
   if(fabs(ap)>fabs(as)) sa=fabs(ap);
   else sa=fabs(as);
   if(pa==0.0) pa=0.5;
   if(sa==0.0) sa=5.0;
   if((*n=(int)(fabs(log(pa/sa)/log(tp/ts))+0.5))<2) *n=2;
   cc=exp(log(pa*sa)/(float)(*n))/(tp*ts);
   c=sqrt(cc);
   
   dp=M_PI_2/(float)(*n);
   *m=(*n)/2;
   k=(*m)*4;
   g=fj=1.0;
   c2=2.0*(1.0-c)*(1.0+c);
   
   for(j=0;j<k;j+=4){
      sj=pow(cos(dp*fj),2.0);
      tj=sin(dp*fj);
      fj=fj+2.0;
      a=1.0/(pow(c+tj,2.0)+sj);
      g=g*a;
      h[j  ]=2.0;
      h[j+1]=1.0;
      h[j+2]=c2*a;
      h[j+3]=(pow(c-tj,2.0)+sj)*a;
   }
   /****  EXIT */
   *gn=g;
   if(*n%2==0) return 0;
   /****  FOR ODD N */
   *m=(*m)+1;
   *gn=g/(1.0+c);
   h[k  ]=1.0;
   h[k+1]=0.0;
   h[k+2]=(1.0-c)/(1.0+c);
   h[k+3]=0.0;
   return 0;
}

int buthip(float *h,int *m,float *gn,int *n,float fp,float fs,float ap,float as)
{
   float wp,ws,tp,ts,pa,sa,cc,c,dp,g,fj,c2,sj,tj,a;
   int k,j;
   if(fabs(fp)>fabs(fs)) wp=fabs(fp)*M_PI;
   else wp=fabs(fs)*M_PI;
   if(fabs(fp)<fabs(fs)) ws=fabs(fp)*M_PI;
   else ws=fabs(fs)*M_PI;
   if(wp==0.0 || wp==ws || wp>=M_PI/2.){
      fprintf(stderr,"? (buthip) invalid input : fp=%14.6e fs=%14.6e ?\n",
              fp,fs);
      return 1;
   }
   /****  DETERMINE N & C */
   tp=tan(wp);
   ts=tan(ws);
   if(fabs(ap)<fabs(as)) pa=fabs(ap);
   else pa=fabs(as);
   if(fabs(ap)>fabs(as)) sa=fabs(ap);
   else sa=fabs(as);
   if(pa==0.0) pa=0.5;
   if(sa==0.0) sa=5.0;
   if((*n=(int)(fabs(log(sa/pa)/log(tp/ts))+0.5))<2) *n=2;
   cc=exp(log(pa*sa)/(float)(*n))*(tp*ts);
   c=sqrt(cc);

   dp=M_PI/2./(float)(*n);
   *m=(*n)/2;
   k=(*m)*4;
   g=fj=1.0;
   c2=(-2.0)*(1.0-c)*(1.0+c);
   
   for(j=0;j<k;j+=4){
      sj=pow(cos(dp*fj),2.0);
      tj=sin(dp*fj);
      fj=fj+2.0;
      a=1.0/(pow(c+tj,2.0)+sj);
      g=g*a;
      h[j  ]=(-2.0);
      h[j+1]=1.0;
      h[j+2]=c2*a;
      h[j+3]=(pow(c-tj,2.0)+sj)*a;
   }
   /****  EXIT */
   *gn=g;
   if(*n%2==0) return 0;
   /****  FOR ODD N */
   *m=(*m)+1;
   *gn=g/(c+1.0);
   h[k  ]=(-1.0);
   h[k+1]=0.0;
   h[k+2]=(c-1.0)/(c+1.0);
   h[k+3]=0.0;
   return 0;
}
/*
+   RECURSIVE FILTERING IN SERIES
+
+   ARGUMENTS
+   X : INPUT TIME SERIES
+   Y : OUTPUT TIME SERIES  (MAY BE EQUIVALENT TO X)
+   N : LENGTH OF X & Y
+   H : COEFFICIENTS OF FILTER
+   M : ORDER OF FILTER
+   NML : >0 ; FOR NORMAL  DIRECTION FILTERING
+       <0 ;   REVERSE DIRECTION FILTERING
+   uv  : past data and results saved
+
+   SUBROUTINE REQUIRED : RECFIL
+
+   M. SAITO  (6/XII/75)
*/
int tandem(float *x,float *y,int n,float *h,int m,int nml,float *uv)
{
    int i;
    if(n<=0 || m<=0){
	fprintf(stderr,"? (tandem) invalid input : n=%d m=%d ?\n",n,m);
	return 1;
    }
    /****  1-ST CALL */
    recfil(x,y,n,h,nml,uv);
    /****  2-ND AND AFTER */
    if(m>1) for(i=1;i<m;i++) recfil(y,y,n,&h[i*4],nml,&uv[i*4]);
    return 0;
}

/*
+   RECURSIVE FILTERING : F(Z) = (1+A*Z+AA*Z**2)/(1+B*Z+BB*Z**2)
+
+   ARGUMENTS
+   X : INPUT TIME SERIES
+   Y : OUTPUT TIME SERIES  (MAY BE EQUIVALENT TO X)
+   N : LENGTH OF X & Y
+   H : FILTER COEFFICIENTS ; H(1)=A, H(2)=AA, H(3)=B, H(4)=BB
+   NML : >0 ; FOR NORMAL  DIRECTION FILTERING
+       <0 ; FOR REVERSE DIRECTION FILTERING
+   uv  : past data and results saved
+
+   M. SAITO  (6/XII/75)
*/
int recfil(float *x,float *y,int n,float *h,int nml,float *uv)
{
   int i,j,jd;
   float a,aa,b,bb,u1,u2,u3,v1,v2,v3;
   if(n<=0){
       fprintf(stderr,"? (recfil) invalid input : n=%d ?\n",n);
       return 1;
   }
   if(nml>=0){
      j=0;
      jd=1;
   }
   else{
      j=n-1;
      jd=(-1);
   }
   a =h[0];
   aa=h[1];
   b =h[2];
   bb=h[3];
   u1=uv[0];
   u2=uv[1];
   v1=uv[2];
   v2=uv[3];
   /****  FILTERING */
   for(i=0;i<n;i++){
      u3=u2;
      u2=u1;
      u1=x[j];
      v3=v2;
      v2=v1;
      v1=u1+a*u2+aa*u3-b*v2-bb*v3;
      y[j]=v1;
      j+=jd;
   }
   uv[0]=u1;
   uv[1]=u2;
   uv[2]=v1;
   uv[3]=v2;
   return 0;
}

int rtr2(float *data,int w, int width)
{
  int i;
  double x1,x2;
  if(width > w/2){fprintf(stderr,"Erro of rtr2\n");return(-1);}
  x1=0;x2=0;
  for (i=0;i<width;i++)
    {
      x1+=(data[i]);
      x2+=(data[w-i-1]);
    }

  x1/=(w/10.);
  x2/=(w/10.);

  for (i=0;i<w;i++)
    {
      data[i]-=
        x1+i*(x2-x1)/(w*1.);
    }
  return(0);
}
