/*
 * ---------------------------------------------------------------------
 * LO or NLO F2 with fast convolution engine, compare with ZMSTF
 * With the convolution engine there are various ways to arrive at
 * the same result (depending on taste and programming style). This
 * is illustrated here by two F2 routines: getFtwo and altFtwo.
 * ---------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

#include "../testjobs/inputpdfCxx.h"

//----------------------------------------------------------------------

void   getFtwo(double *w, int *id, double *def, double *x, double *q,
               double *f, int n, int ichk);
void   altFtwo(double *w, int *id, double *def, double *x, double *q,
               double *f, int n, int ichk);
double timesAs(int *ix, int *iq, int *nf, int *ithresh);
void   myF2wgt(double *w, int nw, int *idw);
double dF2Achi(double *qq);
double dF2Const1(double *xx, double *qq, int *nn);
double dF2C2QB(double *xx, double *qq, int *nn);
double dF2C2GA(double *xx, double *qq, int *nn);

//----------------------------------------------------------------------
int main() {
    int    ityp = 1, iord = 2, nfin = 0;              //unpol, NLO, VFNS
    double as0 = 0.364, r20 = 2.0;                        //input alphas
    double xmin[] = {1e-5,0.2,0.4,0.6,0.75};                    //x-grid
    int    iwt[] = {1,2,4,8,16};                                //x-grid
    int    ng = 5, nxin = 100, iosp = 3;                        //x-grid
    int    nqin = 60;                                         //mu2-grid
    double qq[] = {2e0, 1e4}, wt[] = {1e0, 1e0};              //mu2-grid
    double q2c = 3, q2b = 25, q0 = 2;                 //thresholds, mu20
    
    double def[] =                           //input flavour composition
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
    
    double pro[] = { 4./9., 1./9., 4./9., 1./9., 4./9., 1./9., 0.,
                     1./9., 4./9., 1./9., 4./9., 1./9., 4./9. };
    
    double xpr[] = {1e-4,1e-3,1e-2,1e-1,0.5};
    
    int nx, nq, id1, id2, nw, nfout, ierr ; double eps;
    int lun = 6 ; string outfile = " ";
    
    QCDNUM::qcinit(lun,outfile);                            //initialize
    QCDNUM::gxmake(xmin,iwt,ng,nxin,nx,iosp);                   //x-grid
    QCDNUM::gqmake(qq,wt,2,nqin,nq);                          //mu2-grid
    QCDNUM::wtfile(1,"../weights/unpolarised.wgt");  //calculate weights
    QCDNUM::setord(iord);                                //LO, NLO, NNLO
    QCDNUM::setalf(as0,r20);                              //input alphas
    int iqc  = QCDNUM::iqfrmq(q2c);                    //charm threshold
    int iqb  = QCDNUM::iqfrmq(q2b);                   //bottom threshold
    QCDNUM::setcbt(nfin,iqc,iqb,999);           //thresholds in the VFNS
    int iq0  = QCDNUM::iqfrmq(q0);                         //start scale
    QCDNUM::evolfg(1,func,def,iq0,eps);               //evolve all pdf's

    QCDNUM::zmwfile("../weights/zmstf.wgt");            //ZMSTF  weights
    
    int nww = 7000, idw[3];
    double w[nww];
    myF2wgt(w, nww, idw);                               //my own weights
    
    string s = "  ";
    cout << scientific << setprecision(5);
    cout << "     x" << "           Q2" << "        F2(zmstf) "  <<
    "  F2(getftwo) " << " F2(altftwo)" << endl;
    double qpr[1] = {100};
    double ft[1], f2[1], fa[1];
    for(int ix=0; ix<5; ix++) {
        QCDNUM::zmstfun(2,pro,&xpr[ix],&qpr[0],ft,1,1);
        getFtwo(w,idw,pro,&xpr[ix],&qpr[0],f2,1,1);
        altFtwo(w,idw,pro,&xpr[ix],&qpr[0],fa,1,1);
        cout << xpr[ix] << s << qpr[0] << s << ft[0]
                                       << s << f2[0]
                                       << s << fa[0] << endl;
    }

    return 0;
}

//----------------------------------------------------------------------

void getFtwo(double *w, int *id, double *def, double *x, double *q,
             double *f, int n, int ichk) {

//   This routine uses the fast convolution engine with the
//   perturbative expansion capability of fastFxK. This gives
//   compact code with QCDNUM caring about LO or NLO and
//   multiplication by as(mu2)

    int  idw[] = {0,0,0,0};
      
    int iord;
    QCDNUM::getord(iord);
    if(iord == 3) {
        cout << "GETFTWO: cannot handle NNLO'" << endl;
        exit(1);
    }
//  Setup fast buffers
    QCDNUM::fastini(x,q,n,ichk);
//  F2 quarks
    QCDNUM::fastsns(1,def,7,1);        // 1 = all quarks
    idw[0] = id[0] ;                   // LO
    idw[1] = id[1] ;                   // NLO
    QCDNUM::fastfxk(w,idw,1,3);        // 3 = F2(quarks)
//  F2 qluon
    QCDNUM::fastsns(1,def,0,1);        // 1 = gluon * <e2>
    idw[0] = 0;                        // LO no gluon contribution
    idw[1] = id[2];                    // NLO
    QCDNUM::fastfxk(w,idw,1,2);        // 2 = F2(gluon)
    QCDNUM::fastcpy(2,3,1);            // 3 = F2(quarks)+F2(gluon)
//  Interpolate
    QCDNUM::fastfxq(3,f,n);
}

//----------------------------------------------------------------------

void altFtwo(double *w, int *id, double *def, double *x, double *q,
             double *f, int n, int ichk) {

//  This routine uses the fast convolution engine but not the
//  perturbative expansion capability of fastFxK. Here we need
//  explicit multiplication with alfas, depending on LO or NLO

    int  idw[] = {0,0,0,0};
    
    int iord;
    QCDNUM::getord(iord);
    if(iord == 3) {
        cout << "ALTFTWO: cannot handle NNLO'" << endl;
        exit(1);
    }
//  Setup fast buffers
    QCDNUM::fastini(x,q,n,ichk);
//  F2 LO
    QCDNUM::fastsns(1,def,7,1);          // 1 = all quarks
    QCDNUM::fastcpy(1,4,0);              // 4 = F2(LO)
//  F2 NLO
    if(iord == 2) {
//      NLO quarks
        idw[0] = id[1];                  // C2Q
        QCDNUM::fastfxk(w,idw,1,3);      // 3 = quarks * C2Q
//      NLO qluon
        QCDNUM::fastsns(1,def,0,1);      // 1 = gluon * <e2>
        idw[0] = id[2];                  // C2G
        QCDNUM::fastfxk(w,idw,1,2);      // 2 = gluon * C2G
//      Add NLO quark and gluon and multiply by alfas
        QCDNUM::fastcpy(2,3,1);          // 3 = q*C2Q + g*C2G
        QCDNUM::fastkin(3,timesAs);      // 3 = F2(NLO)
//      Add LO and NLO
        QCDNUM::fastcpy(3,4,1);          // 4 = F2(LO) + F2(NLO)
    }
//  Interpolate
    QCDNUM::fastfxq(4,f,n);
}

//----------------------------------------------------------------------

double timesAs(int *ix, int *iq, int *nf, int *ithresh) {

//  Returns as/2pi  (with correct renormalisation scale dependence!)

    int ierr;
    int jq = *iq;
    return QCDNUM::altabn(0,jq,-1,ierr);
}

//----------------------------------------------------------------------

void myF2wgt(double *w, int nw, int *idw) {

//  Make F2 weights up to NLO
//  Returns identifiers idw[0] = LO, idw[1] = NLO C2q, idw[2] = NLO C2g
    
    int itypes[] = {2,1,0,0,0,0}, npars = 0, newt = 0;
      
//  Book two type-1 tables and one type-2 table
    int isetf2w, nwords;
    QCDNUM::maketab(w,nw,itypes,npars,newt,isetf2w,nwords);
//  LO weights
    idw[0] = 1000*isetf2w+101;
    QCDNUM::makewtd(w,idw[0],dF2Const1,dF2Achi);
//  NLO weights
    idw[1] = 1000*isetf2w+102;
    QCDNUM::makewtb(w,idw[1],dF2C2QB,dF2Achi,0);
    idw[2] = 1000*isetf2w+201;
    QCDNUM::makewta(w,idw[2],dF2C2GA,dF2Achi);
}

//----------------------------------------------------------------------

double dF2Achi(double *qq) {

//  Returns coefficient 'a' of scaling function chi = a*x
//  For zero-mass coefficient functions chi = x and a = 1
    
    return 1.;
}

//----------------------------------------------------------------------

double dF2Const1(double *xx, double *qq, int *nn) {

//  Coefficient of LO delta(1-x) terms
//  Returns the value 1.D0
    
    return 1.;
}

//----------------------------------------------------------------------

double dF2C2QB(double *xx, double *qq, int *nn) {

//  Returns the value of C2Q(x)     singular
    
    double x    = *xx;
    double c1mx = 1 - x;
    return 3. + (5./3.)*x + ((4./3.)*log(c1mx/x)-1.) * (1.+x*x) / c1mx;
}

//----------------------------------------------------------------------

double dF2C2GA(double *xx, double *qq, int *nn) {
    
//  Returns the value of C2G(x,nf) regular
    
    double x    = *xx;
    int nf      = *nn;
    double c1mx = 1 - x;
    double C2g  = -.5 + 4.*x*c1mx + .5 * (x*x+c1mx*c1mx) * log(c1mx/x);
    return  2.*nf*C2g;
}

