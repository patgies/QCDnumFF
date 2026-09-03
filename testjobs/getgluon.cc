/*
 * ---------------------------------------------------------------------
 * Use different interpolation routines to get the gluon
 * ---------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

#include "../testjobs/inputpdfCxx.h"

int main() {
  int    ityp = 1, iord = 3, nfin = 0;                //unpol, NLO, VFNS
  double as0 = 0.364, r20 = 2.0;                          //input alphas
  double xmin[] = {1.e-5, 0.2, 0.4, 0.6, 0.75};                 //x-grid
  int    iwt[] = {1, 2, 4, 8, 16}, ng = 5;                      //x-grid
  int    nxin = 100, iosp = 3;                                  //x-grid
  int    nqin = 60;                                           //mu2-grid
  double qlim[] = {2e0, 1e4}, wt[] = {1e0, 1e0};              //mu2-grid
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
    
  double xx = 1e-4, qq = 1e2;                                   //output
  double xa[] = {1e-4}, qa[] = {1e2}, ga[1];                    //output
  double pdf[13], coef[13];                                     //output
    
  int nx, nq, id1, id2, nw, nfout, ierr ; double eps;
  int lun = 6 ; string outfile = " ";
  string unpfile = "../weights/unpolarised.wgt";
    
  QCDNUM::qcinit(lun,outfile);                              //initialize
  QCDNUM::gxmake(xmin,iwt,ng,nxin,nx,iosp);                     //x-grid
  QCDNUM::gqmake(qlim,wt,2,nqin,nq);                          //mu2-grid
  QCDNUM::wtfile(1,unpfile);                               //get weights
  QCDNUM::setalf(as0,r20);                                //input alphas
  int iqc  = QCDNUM::iqfrmq(q2c);                      //charm threshold
  int iqb  = QCDNUM::iqfrmq(q2b);                     //bottom threshold
  QCDNUM::setcbt(nfin,iqc,iqb,999);             //thresholds in the VFNS
  int iq0  = QCDNUM::iqfrmq(q0);                           //start scale
  QCDNUM::setord(iord);                                      //set order
  QCDNUM::evolfg(1,func,def,iq0,eps);                           //evolve
  
  double gluon;
  gluon = QCDNUM::bvalxq(1,0,xx,qq,1);
  cout << scientific << setprecision(4);
  cout << " BVALXQ gluon =    " << gluon << endl;
  gluon = QCDNUM::fvalxq(1,0,xx,qq,1);
  cout << " FVALXQ gluon =    " << gluon << endl;
  QCDNUM::allfxq(1,xx,qq,pdf,0,1);
  cout << " ALLFXQ gluon =    " << pdf[6] << endl;
  gluon = QCDNUM::sumfxq(1,coef,0,xx,qq,1);
  cout << " SUMFXQ gluon =    " << gluon << endl;
  QCDNUM::fflist(1,coef,0,xa,qa,ga,1,1);
  cout << " FFLIST gluon =    " << ga[0] << endl;
  QCDNUM::ftable(1,coef,0,xa,1,qa,1,ga,1);
  cout << " FTABLE gluon =    " << ga[0] << endl << endl;

  int ix = QCDNUM::ixfrmx(xx);
  int iq = QCDNUM::iqfrmq(qq);
  gluon = QCDNUM::bvalij(1,0,ix,iq,1);
  cout << " BVALIJ gluon =    " << gluon << endl;
  gluon = QCDNUM::fvalij(1,0,ix,iq,1);
  cout << " FVALIJ gluon =    " << gluon << endl;
  QCDNUM::allfij(1,ix,iq,pdf,0,1);
  cout << " ALLFIJ gluon =    " << pdf[6] << endl;
  gluon = QCDNUM::sumfij(1,coef,0,ix,iq,1);
  cout << " SUMFIJ gluon =    " << gluon << endl;

  return 0;
}
