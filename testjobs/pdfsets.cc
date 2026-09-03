/*
 * ---------------------------------------------------------------------
 * Play around with pdf sets in memory
 *
 * set 1 : unpolarised NNLO
 * set 2 : polarised   NLO
 * set 5 : unpolarised LO
 * set 6 : unpolarised NLO
 * set 7 : unpolarised NNLO
 * set 8 : unpolarised NNLO + 2 extra pdfs (singlet and proton)
 *
 * This tests if PDFCPY, PDFEXT, FTABLE, FFLIST are working correctly
 * ---------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

#include "../testjobs/inputpdfCxx.h"

/*----------------------------------------------------------------------
 * Return q,qbar,g for ipdf = [-6:6], proton singlet [7], nonsinglet [8]
 *----------------------------------------------------------------------
 */
 double fimport(int* ii, double* xx, double* qq, bool* fst) {
  static double proton[13];
  int ipdf = *ii;
  double x = *xx;
  double q = *qq;
  bool first = *fst;
    if(first) {
    for(int i=1; i<=13; i++) { QCDNUM::qstore("Read",i,proton[i-1]); }
    }
  double f = 0;
  if(ipdf >= -6 && ipdf <= 6) f = QCDNUM::fvalxq(1,ipdf,x,q,1);
  if(ipdf == 7)               f = QCDNUM::sumfxq(1,proton,2,x,q,1);
  if(ipdf == 8)               f = QCDNUM::sumfxq(1,proton,3,x,q,1);
  return f;
}

/*----------------------------------------------------------------------
 * Pointer function
 *----------------------------------------------------------------------
 */
 inline int kij(int i, int j, int n) { return i-1 + n*(j-1); }

/*----------------------------------------------------------------------
 * Main program
 *----------------------------------------------------------------------
 */

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
    
  double proton[] = {4, 1, 4, 1, 4, 1, 0, 1, 4, 1, 4, 1, 4};    //proton
    
  double xx[] = {1e-5,1e-4,1e-3,1e-2,1e-1};                     //output
  double qq[] = {1e1,1e2,1e3};                                  //output
  double xf1[15],xf2[15],xf8[15],psi[15],pns[15],epar[13];      //output
    
  int nx, nq, id1, id2, nw, nfout, ierr ; double eps;
  int lun = 6 ; string outfile = " ";
  string unpfile = "../weights/unpolarised.wgt";
  string polfile = "../weights/polarised.wgt";
    
  for(int i=0; i<=13; i++) {proton[i] /= 9;}        //divide proton by 9
    
  QCDNUM::qcinit(lun,outfile);                              //initialize
    
  for(int i=0 ; i<=13; i++) {
      QCDNUM::qstore("Write",i+1,proton[i]); }            //store proton
    
  QCDNUM::gxmake(xmin,iwt,ng,nxin,nx,iosp);                     //x-grid
  QCDNUM::gqmake(qlim,wt,2,nqin,nq);                          //mu2-grid
    
  QCDNUM::wtfile(1,unpfile);                             //weights unpol
  QCDNUM::wtfile(2,polfile);                               //weights pol
    
  QCDNUM::setalf(as0,r20);                                //input alphas
  int iqc  = QCDNUM::iqfrmq(q2c);                      //charm threshold
  int iqb  = QCDNUM::iqfrmq(q2b);                     //bottom threshold
  QCDNUM::setcbt(nfin,iqc,iqb,999);             //thresholds in the VFNS
  int iq0  = QCDNUM::iqfrmq(q0);                           //start scale
    
//Unpolarsed LO, NLO and NNLO in sets 5, 6 and 7 -----------------------
  for(iord=1; iord<=3; iord++) {
    QCDNUM::setord(iord);
    QCDNUM::evolfg(1,func,def,iq0,eps);
    QCDNUM::pdfcpy(1,4+iord);
    }
  
//Print gluon table from set 1 -----------------------------------------
  double gt1[15];
  QCDNUM::ftable(1,proton,0,xx,5,qq,3,gt1,1);
  cout << scientific << setprecision(4);
  string s = "           ";
  cout << " Gluon table IX =  1";
  cout << s << 2 << s << 3 << s << 4 << s << 5 << endl;
  s = "  ";
  for(int j=1; j<=3; j++) {                          //print gluon table
    cout << " IQ =  " << j ;
    for(int i = 1; i<=5; i++) { cout << s << gt1[ kij(i,j,5) ]; }
    cout << endl;
    }
  cout << endl;
 
//Make a gluon list of set 1 -------------------------------------------
  double xl[15], ql[15], gl1[15];
  int k = 0;
  for(int j=0; j<3; j++) {
    for(int i=0; i<5; i++) { xl[k] = xx[i]; ql[k] = qq[j]; k++;}
    }
  QCDNUM::fflist(1,proton,0,xl,ql,gl1,15,1);
    
//Check that table = list ----------------------------------------------
  double dif = -10;
  for(int i=1; i<15; i++) { dif = max( dif, fabs(gl1[i]-gt1[i]) ); }
  cout << " This should be small:  " << dif << endl;
  
//Polarised pdfs at NLO in default set 2  ------------------------------
  QCDNUM::setord(2);
  QCDNUM::evolfg(2,func,def,iq0,eps);
    
//Import NNLO unpol (from set 1) into set 8 with 2 extra pdfs ----------
//Make sure that extpdf uses the params of set 1
  QCDNUM::pushcp();                                          //save pars
  QCDNUM::usepar(1);                                    //use pars set 1
  double offset = 0;
//s/r fimport stores proton singlet and nonsinglet in idpdf 7 and 8
  QCDNUM::extpdf(fimport,8,2,offset,eps);
  QCDNUM::pullcp();                                       //restore pars
    
//Make proton table of set 1 and 8 -------------------------------------
  QCDNUM::ftable(1,proton,1,xx,5,qq,3,xf1,1);         //proton     set 1
  QCDNUM::ftable(8,proton,1,xx,5,qq,3,xf8,1);         //proton     set 8
  QCDNUM::ftable(8,proton,13,xx,5,qq,3,psi,1);        //singlet    set 8
  QCDNUM::ftable(8,proton,14,xx,5,qq,3,pns,1);        //nonsinglet set 8

//See if set 1 is properly imported into set 8  ------------------------
  dif = -10;
  for(int i=1; i<15; i++) { dif = max( dif, fabs(xf8[i]-xf1[i]) ); }
  cout << " This should be small:  " << dif << endl;
    
//See if proton singlet and nonsinglet are properly stored in 8 --------
  dif = -10;
  for(int i=1; i<15; i++) {
      dif = max( dif, fabs(xf8[i]-psi[i]-pns[i]) ); }
  cout << " This too            :  " << dif << endl << endl;
    
//Now print proton table of set 2 --------------------------------------
  QCDNUM::ftable(2,proton,1,xx,5,qq,3,xf2,1);             //proton set 2
  s = "           ";
  cout << " Proton table IX = 1";
  cout << s << 2 << s << 3 << s << 4 << s << 5 << endl;
  s = "  ";
  for(int j=1; j<=3; j++) {                         //print proton table
    cout << " IQ =  " << j ;
    for(int i = 1; i<=5; i++) { cout << s << xf2[ kij(i,j,5) ]; }
    cout << endl;
    }
  cout << endl;
    
//Print some parameters of all sets in memory --------------------------
  cout << " Parameters of    ISET:   IORD    ITYP     NPDF     KEY";
  cout << endl;
  s = "       ";
  for(int iset=1; iset<=8; iset++){
    int npdfs = QCDNUM::nptabs(iset);
    if(npdfs){                                        //check set exists
      QCDNUM::cpypar(epar,13,iset);
      int iorder = epar[0];
      int ipdfs  = epar[12];
      int key    = QCDNUM::keypar(iset);
      cout << s << s << s << iset << s << iorder << s << ipdfs;
      cout << s << npdfs << s << key << endl;
      }
    }
    
  return 0;
}
