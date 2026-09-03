/*
 * ---------------------------------------------------------------------
 * 23-01-17  Basic QCDNUM example job in C++
 * ---------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

//----------------------------------------------------------------------
double xupv(double x) {
  double au = 5.107200;
  double pd = au * pow(x,0.8) * pow((1-x),3);
  return pd;
}
//----------------------------------------------------------------------
double xdnv(double x) {
  double ad = 3.064320;
  double pd = ad * pow(x,0.8) * pow((1-x),4);
  return pd;
}
//----------------------------------------------------------------------
double xglu(double x) {
  double ag = 1.7;
  double pd = ag * pow(x,-0.1) * pow((1-x),5);
  return pd;
}
//----------------------------------------------------------------------
double xdbar(double x) {
  double adbar = 0.1939875;
  double pd = adbar * pow(x,-0.1) * pow((1-x),6);
  return pd;
}
//----------------------------------------------------------------------
double xubar(double x) {
  double pd = xdbar(x) * (1-x);
  return pd;
}
//----------------------------------------------------------------------
double xsbar(double x) {
  double pd = 0.2 * ( xdbar(x) + xubar(x) );
  return pd;
}
//----------------------------------------------------------------------
double func(int* ipdf, double* x) {
  int i = *ipdf;
  double xb = *x;
  double f = 0;
  if(i ==  0) f = xglu(xb);
  if(i ==  1) f = xdnv(xb);
  if(i ==  2) f = xupv(xb);
  if(i ==  3) f = 0;
  if(i ==  4) f = xdbar(xb);
  if(i ==  5) f = xubar(xb);
  if(i ==  6) f = xsbar(xb);
  if(i ==  7) f = 0;
  if(i ==  8) f = 0;
  if(i ==  9) f = 0;
  if(i == 10) f = 0;
  if(i == 11) f = 0;
  if(i == 12) f = 0;
  return f;
}
//----------------------------------------------------------------------
int main() {
  int    ityp = 1, iord = 3, nfin = 0;                //unpol, NLO, VFNS
  double as0 = 0.364, r20 = 2.0;                          //input alphas
  double xmin[] = {1.e-4};                                      //x-grid
  int    iwt[] = {1}, ng = 1, nxin = 100, iosp = 3;             //x-grid
  int    nqin = 60;                                           //mu2-grid
  double qq[] = { 2e0, 1e4}, wt[] = { 1e0, 1e0};              //mu2-grid
  double q2c = 3, q2b = 25, q0 = 2;                   //thresholds, mu20
  double x = 1e-3, q = 1e3, qmz2 = 8315.25, pdf[13];            //output
    
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
    
  int nx, nq, id1, id2, nw, nfout, ierr ; double eps;
  int lun = 6 ; string outfile = " ";
    
  QCDNUM::qcinit(lun,outfile);                              //initialize
  QCDNUM::gxmake(xmin,iwt,ng,nxin,nx,iosp);                     //x-grid
  QCDNUM::gqmake(qq,wt,2,nqin,nq);                            //mu2-grid
  QCDNUM::wtfile(1,"../weights/unpolarised.wgt");    //calculate weights
  QCDNUM::setord(iord);                                  //LO, NLO, NNLO
  QCDNUM::setalf(as0,r20);                                //input alphas
  int iqc  = QCDNUM::iqfrmq(q2c);                      //charm threshold
  int iqb  = QCDNUM::iqfrmq(q2b);                     //bottom threshold
  QCDNUM::setcbt(nfin,iqc,iqb,999);             //thresholds in the VFNS
  int iq0  = QCDNUM::iqfrmq(q0);                           //start scale
  QCDNUM::evolfg(1,func,def,iq0,eps);                 //evolve all pdf's
  QCDNUM::allfxq(1,x,q,pdf,0,1);                 //interpolate all pdf's
    
  double csea = 2 * pdf[2];                         //charm sea at x,mu2
  double asmz = QCDNUM::asfunc(qmz2,nfout,ierr);           //alphas(mz2)

//Format with C++
  string s = " ";
  cout << scientific << setprecision(4);
  cout << "x, q, CharmSea = " << x << s << q << s << csea << endl;
  cout << "as(mz2)        = " << asmz << endl;

//    Format with MBUTIL
//    string s = " ";
//    string tx, tq, tcs, tas;
//    int len;
//    MBUTIL::smb_dtoch(x,5,tx,len);
//    MBUTIL::smb_dtoch(q,5,tq,len);
//    MBUTIL::smb_dtoch(csea,5,tcs,len);
//    MBUTIL::smb_dtoch(asmz,5,tas,len);
//    cout << "x, q, Charmsea = " << tx << s << tq << s << tcs << endl;
//    cout << "as(mz2)        = " << tas << endl;

  return 0;
}
