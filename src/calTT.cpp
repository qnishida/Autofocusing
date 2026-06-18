#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;
const double D2R = M_PI/180.;

static double *ttP, *pP, *RP2,*RS2;
static double *ttS, *pS;
static double *ttSS,*pSS;

int init_ttTable(){
  ttP = new double [1000];
  pP  = new double [1000];
  RP2  = new double [1000];//Square of geometrical spreading factor R of P-wave
  RS2  = new double [1000];//Square of geometrical spreading factor R of S-wave
  ttS = new double [1000];
  pS  = new double [1000];
  ttSS= new double [2000];
  pSS = new double [2000];

  double theta,tt0,p0,dp0;
  ifstream ifs("TravelTime/IASP91/TravelTimeS");
  while(ifs){
    ifs >> theta >> tt0 >> p0 >> dp0;
    int ith=(int)(theta/0.1+.5);
    ttS[ith]=tt0;
    pS[ith]=p0/D2R;//in sec/radian
    RS2[ith]=sin(theta*D2R)/pS[ith]/fabs(dp0);//in sec/radian^2
  }
  ifs.close();

  //Data Angular distance: Travel time: Slowness [s/deg]: dpdth [s/rad/rad]
  ifs.open("TravelTime/IASP91/TravelTimeP");
  while(ifs){
    ifs >> theta >> tt0 >> p0 >> dp0;
    int ith=(int)(theta/0.1+.5);
    ttP[ith]=tt0;
    pP[ith]=p0/D2R;//in sec/radian
    RP2[ith]=sin(theta*D2R)/pP[ith]/fabs(dp0);//in radian^3/s^2
  }
  ifs.close();

  ifs.open("TravelTime/IASP91/TravelTimeSS");
  while(ifs){
    ifs >> theta >> tt0 >> p0;
    int ith=(int)(theta/0.1+.5);
    ttSS[ith]=tt0;
    pSS[ith]=p0/D2R;//in sec/radian
  }
  return 0;
}

double cal_RP2(double theta){//in degree
  if(theta > 45 && theta <98){
    int ith =(int)(theta/0.1);
    double dth = theta-ith*0.1;
    double RP2_0 =RP2 [ith ];
    double RP2_1 =RP2 [ith+1];
    return((RP2_0*(0.1-dth)+RP2_1*dth)/0.1);
  }
  return 0.;
}

double cal_RS2(double theta){//in degree
  if(theta > 45 && theta <95){
    int ith =(int)(theta/0.1);
    double dth = theta-ith*0.1;
    double RS2_0 =RS2 [ith ];
    double RS2_1 =RS2 [ith+1];
    return((RS2_0*(0.1-dth)+RS2_1*dth)/0.1);
  }
  return 0.;
}


double p2gcarc(double p, int flag){//Input: in s/rad, Out: s/degree
  if(flag==0){//P-wave
    for(int i=(int)(28/.1+.5);i*.1<95.;i++){
      if(pP[i]>=p && pP[i+1]< p){
	return((pP[i]-p)/(pP[i]-pP[i+1])*0.1+i*0.1 );
      }
    }
  }
  else if(flag==1){//S-wave
    for(int i=(int)(28/.1+.5);i*.1<95.;i++){
      if(pS[i]>=p && pS[i+1]< p){
	return((pS[i]-p)/(pS[i]-pS[i+1])*0.1+i);
      }
    }
  }
  else
    return(0);
  return 0.;
}


double cal_Slowness(double theta, int flag){//Input: in degree, Out: s/degree
  if(flag==0){//P-wave
    if(theta > 45 && theta <98){
      int ith =(int)(theta/0.1);
      double dth = theta-ith*0.1;
      //double ttP0=ttP[ith];
      double pP0 =pP [ith ];
      double pP1 =pP [ith+1];
      return((pP0*(0.1-dth)+pP1*dth)/0.1);
    }
  }
  else if(flag==1){//S-wave
    if(theta > 45 && theta <98){
      int ith =(int)(theta/0.1);
      double dth = theta-ith*0.1;
      //double ttS0=ttS[ith];
      double pS0 =pS [ith   ];
      double pS1 =pS [ith + 1];
      return((pS0*(0.1-dth)+pS1*dth)/0.1);
    }
  }
  else if(flag==2){//SS-wave
    if(theta > 80 && theta <150){
      int ith =(int)(theta/0.1);
      double dth = theta-ith*0.1;
      //double ttSS0=ttSS[ith];
      double pSS0 =pSS [ith   ];
      double pSS1 =pSS [ith +1];
      return((pSS0*(0.1-dth)+pSS1*dth)/0.1);
    }
  }
  else
    return(0);
  return 0.;
}

double cal_TT(double theta, int flag){//in Degree
  if(flag==0){//P-wave
    if(theta > 45 && theta <98){
      int ith =(int)(theta/0.1);
      double dth = theta-ith*0.1;
      //double ttP0=ttP[ith];
      double ttP0 =ttP [ith ];
      double ttP1 =ttP [ith+1];
      return((ttP0*(0.1-dth)+ttP1*dth)/0.1);
    }
  }
  else if(flag==1){//S-wave
    if(theta > 45 && theta <98){
      int ith =(int)(theta/0.1);
      double dth = theta-ith*0.1;
      //double ttS0=ttS[ith];
      double ttS0 =ttS [ith   ];
      double ttS1 =ttS [ith + 1];
      return((ttS0*(0.1-dth)+ttS1*dth)/0.1);
    }
  }
  else if(flag==2){//SS-wave
    if(theta > 80 && theta <150){
      int ith =(int)(theta/0.1);
      double dth = theta-ith*0.1;
      //double ttSS0=ttSS[ith];
      double ttSS0 = ttSS [ith   ];
      double ttSS1 = ttSS [ith +1];
      return((ttSS0*(0.1-dth)+ttSS1*dth)/0.1);
    }
  }
  else
    return(0);
  return 0.;
}
