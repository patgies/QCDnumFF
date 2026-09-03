/*
 * ---------------------------------------------------------------------
 * Enjoy the speed of QCDNUM and ZMSTF
 * ---------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <random>
#include <time.h>
#include "QCDNUM/QCDNUM.h"
using namespace std;

#include "../testjobs/inputpdfCxx.h"

int main() {
  int    iord = 3, nfin = 0;                                //NNLO, VFNS
  double as0 = 0.364, r20 = 2.0;                          //input alphas
  double xmin[] = {1.e-5, 0.2, 0.4, 0.6, 0.75};                 //x-grid
  int    iwt[] = {1, 2, 4, 8, 16};                              //x-grid
  int    ngx = 5, nxin = 100, iosp = 3;                         //x-grid
  double qlim[] = {2e0, 1e4}, wt[] = {1e0, 1e0};              //mu2-grid
  int    ngq = 2, nqin = 50;                                  //mu2-grid
  double q2c = 3, q2b = 25, q0 = 2;                   //thresholds, mu20

  double def[] =                             //input flavour composition
  // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
    { 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,      // 1=dval
      0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,      // 2=uval
      0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,      // 3=sval
      0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,      // 4=dbar
      0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,      // 5=ubar
      0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // 6=sbar
      0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,      // 7=cval
      0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // 8=cbar
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // 9=zero
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      //10=zero
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      //11=zero
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};     //12=zero
    
  double xx[1000], q2[1000], ff[1000];                          //output
  double proton[] = {4.,1.,4.,1.,4.,1.,0.,1.,4.,1.,4.,1.,4.};   //proton
  for(int i = 0; i<13; i++) proton[i] /= 9;                   //proton:9
  
//Generate 1000 random x-mu2 points in the HERA kinematic range
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  double roots = 300;                                  //HERA COM energy
  double shera = roots*roots;
  double xlog1 = log(xmin[0]);
  double xlog2 = log(0.99);
  double qlog1 = log(qlim[0]);
  double qlog2 = log(qlim[ngq-1]);
  int  ntot  = 0;
  while(ntot < 1000) {
    double rval  = distribution(generator);
    double xlog  = xlog1 + rval*(xlog2-xlog1);
    double xxxx  = exp(xlog);
    rval         = distribution(generator);
    double qlog  = qlog1 + rval*(qlog2-qlog1);
    double qqqq  = exp(qlog);
    if(qqqq <= xxxx*shera) {
      xx[ntot] = xxxx;
      q2[ntot] = qqqq;
      ntot += 1;
    }
  }

//Setup evolution and structure function calculation
  int nx, nq, id1, id2, nw, nfout, ierr ; double eps;
  int lun = 6 ; string outfile = " ";
  string unpfile = "../weights/unpolarised.wgt";
  string zmsfile = "../weights/zmstf.wgt";
    
  QCDNUM::qcinit(lun,outfile);                              //initialize
  QCDNUM::gxmake(xmin,iwt,ngx,nxin,nx,iosp);                    //x-grid
  QCDNUM::gqmake(qlim,wt,ngq,nqin,nq);                        //mu2-grid
    

  QCDNUM::wtfile(1,unpfile);                             //weights unpol

  int lunw = QCDNUM::nxtlun(10);  
  QCDNUM::zmreadw(lunw,zmsfile,nw,ierr);                 //weights zmstf
  if(ierr) { QCDNUM::zmfillw(nw);
             QCDNUM::zmdumpw(lunw,zmsfile); }

  QCDNUM::setord(iord);                                  //LO, NLO, NNLO
  QCDNUM::setalf(as0,r20);                                //input alphas
  int iqc  = QCDNUM::iqfrmq(q2c);                      //charm threshold
  int iqb  = QCDNUM::iqfrmq(q2b);                     //bottom threshold
  QCDNUM::setcbt(nfin,iqc,iqb,999);             //thresholds in the VFNS
  int iq0  = QCDNUM::iqfrmq(q0);                           //start scale

//Evolution and structure function timing loop
  cout << " Wait: 1000 evols and 2.10^6 stfs will take ... " << endl;
  const clock_t tim1 = clock();
  for(int iter=0; iter < 1000; iter++ ) {
    QCDNUM::evolfg(1,func,def,iq0,eps);                         //evolve
    QCDNUM::zmstfun(1,proton,xx,q2,ff,1000,1);               //FL proton
    QCDNUM::zmstfun(2,proton,xx,q2,ff,1000,1);               //F2 proton
    }
  const clock_t tim2 = clock();

  cout << setprecision(4) << fixed;
  cout << " CPU : " << float(tim2 - tim1)/CLOCKS_PER_SEC << " secs";
  cout << endl;
  
  return 0;
}
