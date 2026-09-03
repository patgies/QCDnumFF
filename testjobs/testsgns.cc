/*
 * -----------------------------------------------------------------------------
 * Compare EVOLSG andf EVSGNS
 * -----------------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

#include "../testjobs/inputpdfCxx.h"

//------------------------------------------------------------------------------

void compa(int jset1, int jset2, int id, int jq, int ixmin);
double func2(int *ipdf, double *x);

//------------------------------------------------------------------------------

int main() {
    
//  LO/NLO/NNLO,  unpol/pol/frag,  starting scale, debug printout
    int  iord = 3, itype = 1, idbug = 0;
    double q0 = 2.0;

//  FFNS,VFNS,MFNS and flavour thresholds
    int   nfix = 6;
    double q2c = 3.0, q2b = 25.0, q2t = 1.e11;

//  Input alphas
    double as0 = 0.364, r20 = 2.0;

//  Renormalisation scale and factorisation scale
    double aar = 1.0, bbr = 0.0, aaq = 1.0, bbq = 0.0;

//  x-grid
    double xmin[5] = {1.e-5, 0.2, 0.4, 0.6, 0.75};
    int iwt[5] = {1, 2, 4, 8, 16};
    int nxsubg = 5, nxin = 100, iosp = 3;

//  mu2-grid
    double qq[2] = {2.0, 1e4}, wt[2] = {1., 1.};
    int nqin = 60;

//  Input parton distributions and flavour composition of these
    double def[] = {
//  tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
//  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
     0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,     //dval
     0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,     //uval
     0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,     //sval
     0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,     //dbar
     0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,     //ubar
     0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,     //sbar
     0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,     //cval
     0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,     //cbar
     0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,     //bval
     0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,     //bbar
    -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,     //tval
     1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };   //tbar
    
//  Weight file
    string fnamw[3] = {
        "../weights/unpolarised.wgt",
        "../weights/polarised.wgt  ",
        "../weights/timelike.wgt   " };

//  Evolution type definitions for EVSGNS: 1=SG, -1=V, 2=NS+, -2=NS-
//                   1   2   3   4   5   6   7   8   9  10  11  12
    int isns[12] = { 1,  2,  2,  2,  2,  2, -1, -2, -2, -2, -2, -2 };

//  ----------------------------------------------------------------------------
//  QCDNUM set-up
//  ----------------------------------------------------------------------------

    int lun    = 6, nx, nq;
    QCDNUM::qcinit(lun," ");                                   // initialisation
    QCDNUM::gxmake(xmin,iwt,nxsubg,nxin,nx,iosp);                      // x-grid
    QCDNUM::gqmake(qq,wt,2,nqin,nq);                                 // mu2-grid
    QCDNUM::wtfile(itype,fnamw[itype-1]);                  // make weight tables
    QCDNUM::setord(iord);                                       // LO, NLO, NNLO
    QCDNUM::setalf(as0,r20);                                     // input alphas
    int iqc = QCDNUM::iqfrmq(q2c);                            // charm threshold
    int iqb = QCDNUM::iqfrmq(q2b);                           // bottom threshold
    if(nfix >= 0) {
        QCDNUM::setcbt(nfix,iqc,iqb,0);                // thresholds in the VFNS
    }
    else {
        QCDNUM::mixfns(abs(nfix),q2c,q2b,0.0);         // thresholds in the MFNS
    }
    QCDNUM::setabr(aar,bbr);                            // renormalisation scale
    int iq0  = QCDNUM::iqfrmq(q0);                             // starting scale
    
//  Print settings
    printf(" xmin = %-10G",xmin[0]); printf("        subgrids = %i \n",nxsubg);
    printf(" qmin = %-10G",qq[0]);   printf("        qmax     = %G \n",qq[1]);
    printf(" nfix = %-10i",nfix);
    printf("        q2c,b,t  = %G %G %G \n",q2c,q2b,q2t);
    printf(" alfa = %-10G",as0);     printf("        r20      = %G \n",r20);
    printf(" aar  = %-10G",aar);     printf("        bbr      = %G \n",bbr);
    printf(" q20  = %-10G \n \n ",q0);
    
//  ----------------------------------------------------------------------------
//  Evolutions
//  ----------------------------------------------------------------------------

//  Evolfg
    double eps;
    int iset1  = 1;                                          // evolve pdf set 1
    int jtype  = 10*iset1+itype;
    QCDNUM::evolfg(jtype,func,def,iq0,eps);

//  Evsgns func2 gets the start values at iq0 from iset1
    double val = double(iset1);
    QCDNUM::qstore("write",1,val);                // pass iset1 via common block
    val = double(iq0);
    QCDNUM::qstore("write",2,val);                  // pass iq0 via common block
    
//  Evsgns
    int iset2 = 2;                                           // evolve pdf set 2
    jtype = 10*iset2+itype;
    QCDNUM::setint("edbg",idbug);
    QCDNUM::evsgns(jtype,func2,isns,12,iq0,eps);
    
//  ----------------------------------------------------------------------------
//  Compare
//  ----------------------------------------------------------------------------

    cout << endl << " Compare EVOLFG and EVSGNS" << endl << endl;
    int id1 = 0;
    int id2 = 12;
    int iqq = 0;
    int ixm = 1;
    for(int id = id1; id <= id2; id++) {
        compa(iset1,iset2,id,iqq,ixm);
    }
    cout << endl;

}

//  ----------------------------------------------------------------------------

void compa(int jset1, int jset2, int id, int jq, int ixmin) {

//  Compare basis pdfs
//
//  jset1,2   (in) : pdf sets to compare
//  id        (in) : basis pdf id [0,12]
//  jq        (in) : 0 = all mu2; != 0 --> compare for this jq only
//  ixmin     (in) : print range [1,ixmin] when jq != 0

    int nx, nq, iosp;
    double xmi, xma, qmi, qma;
    QCDNUM::grpars(nx,xmi,xma,nq,qmi,qma,iosp);

    int iq1, iq2, jxmin;
    if(jq == 0) {
        iq1   = 1;
        iq2   = nq;
        jxmin = 0;
    }
    else {
        iq1   = jq;
        iq2   = jq;
        jxmin = ixmin;
    }

    double dif = 0.0;
    for(int iq = iq1; iq <= iq2; iq++) {
        for (int ix = 1; ix <= nx; ix++) {
            double val1 = QCDNUM::bvalij(jset1,id,ix,iq,1);
            double val2 = QCDNUM::bvalij(jset2,id,ix,iq,1);
            double difi = val1-val2;
            dif  = max(dif,abs(difi));
            if(ix <= jxmin) {
                printf("%3i,%3i,%10.4F,%10.4F,%13.3E,%5i,%3i,%3i,\n",
                       id,ix,val1,val2,difi,jq,jset1,jset2);
            }
        }
    }
    if(jxmin == 0) {
            printf(" id = %2i",id);  printf(" dif = %15.5E \n",dif);
    }
}

//  ----------------------------------------------------------------------------

double func2(int *ipdf, double *x) {                                //for EVSGNS

//  This function takes the input pdfs from iset

    static int iset, iq0;
    double func = 0.0;

    if(*ipdf == -1) {
//      Read pdf set and input scale
        double val;
        QCDNUM::qstore("read",1,val);
        iset = int(val);
        QCDNUM::qstore("read",2,val);
        iq0  = int(val);
    }
    else {
//     Get start pdf from iset
        int ix = QCDNUM::ixfrmx(*x);
        func = QCDNUM::bvalij(iset,*ipdf,ix,iq0,1);
    }

    return func;
}
