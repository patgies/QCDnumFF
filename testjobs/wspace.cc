/*
 * ---------------------------------------------------------------------
 *     program wspace
 * ---------------------------------------------------------------------
 *     Examples described in the WSTORE write-up
 * ---------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;
using namespace WSTORE;

//- Routine protptypes -------------------------------------------------

int NW = 20000;
int NT = 10;
void mytab(double *w, int &ias, int &iax, int &iaq, int &ixq);
void K3(double *table, int (&kk)[5]);
int iP3(double *table, int i, int j, int k);
double val3(double *table, int i, int j, int k);
double val1(double *table, int i);
void fastloop(double *table, int i1, int i2, int j1,
                             int j2, int k1, int k2);
void wtable(double *w, int &ias, int &iaw);
void IniTab(double *table, double val);

//----------------------------------------------------------------------
int main() {
        
    cout << endl << " You are using WSTORE version " <<
                                  iws_Version() << " ---------" << endl;
    
    double w[NW];
    int iaS = iws_WsInit(w,NW,NT,"Please increase parameter NW");
    
    cout << endl <<
             " Example of a table-set ------------------------" << endl;
    int iax, iaq, ixq;
    mytab(w,iaS,iax,iaq,ixq);
    double *xtab = w+iax;
    double *qtab = w+iaq;
    double *tab3 = w+ixq;
    cout << " Limits of x-bin 10 :  " << val1(xtab,10) << "  "
                                      << val1(xtab,11) << endl;
    cout << " Limits of q-bin 15 :  " << val1(qtab,15) << "  "
                                      << val1(qtab,16) << endl;
    cout << " Value of T(10,5,4) :  " << val3(tab3,10,5,4) << endl;
    
    cout << endl <<
             " Fast 3-dim loop construct ---------------------" << endl;
    int i1 = 11;
    int i2 = 13;
    int j1 = 15;
    int j2 = 17;
    int k1 = 3;
    int k2 = 4;
    fastloop(tab3,i1,i2,j1,j2,k1,k2);

    cout << endl <<
             " Convolution integral --------------------------" << endl;
    int iaw;
    wtable(w,iaS,iaw);
    double *wtab = w+iaw;
    int ix  = 10;
    int iq  = 20;
    int nf  =  4;
    int jaF = iP3(tab3,1,iq,nf); //Address of F(1,iq,nf);
    int jaW = iP3(wtab,1,ix,nf); //Address of W(1,ix,nf);
    double FxW = MBUTIL::dmb_vdotv(tab3+jaF,wtab+jaW,ix);
    cout << " T(10,20,4) x W     :  " << FxW << endl;
    
    cout << endl <<
             " Initialise table ------------------------------" << endl;
    IniTab(tab3,-1);
    fastloop(tab3,i1,i2,j1,j2,k1,k2);
    
    return 0;
}

//----------------------------------------------------------------------
void mytab(double *w, int &ias, int &iax, int &iaq, int &ixq) {
    int imi[] = {1,1,3};
    int ima[] = {50,25,6};
    
    if(NT < 3) {
        cout << "MYTAB: increase parameter NT to at least 3";
        exit(1);
    }
//  Create table set
    ias      = iws_NewSet(w);
    int imax = ima[0]+1;
    iax      = iws_WTable(w,&imi[0],&imax,1);
    imax     = ima[1]+1;
    iaq      = iws_WTable(w,&imi[1],&imax,1);
    ixq      = iws_WTable(w,imi,ima,3);
    int iat  = iws_IaFirstTag(w,ias);
    w[iat]   = double(iax-ias);
    w[iat+1] = double(iaq-ias);
    w[iat+2] = double(ixq-ias);
//  Fill x-table
    int ia = iws_BeginTbody(w,iax)-1;
    for(int i=imi[0]; i<=ima[0]+1; i++){
        ia += 1;
        w[ia] = double(i);
    }
    int ja = iws_EndTbody(w,iax);
    if(ia != ja) {
        cout << "MYTAB: something wrong with x-indexing" << endl;
        exit(1);
    }
//  Fill q-table
    ia = iws_BeginTbody(w,iaq)-1;
    for(int i=imi[1]; i<=ima[1]+1; i++){
        ia += 1;
        w[ia] = double(i);
    }
    ja = iws_EndTbody(w,iaq);
    if(ia != ja) {
        cout << "MYTAB: something wrong with q-indexing" << endl;
        exit(1);
    }
//  Fill xq-table
    ia = iws_BeginTbody(w,ixq)-1;
    for(int k=imi[2]; k<=ima[2]; k++){
        for(int j=imi[1]; j<=ima[1]; j++){
            for(int i=imi[0]; i<=ima[0]; i++){
                ia += 1;
                w[ia] = double(10000*i+100*j+k);
            }
        }
    }
    ja = iws_EndTbody(w,ixq);
    if(ia != ja) {
        cout << "MYTAB: something wrong with xq-indexing" << endl;
        exit(1);
    }
}

//----------------------------------------------------------------------
void K3(double *table, int (&kk)[5]) {
    int ia = int(*(table+1));
    double *w = table-ia;
    int iak   = iws_IaKARRAY(w,ia);
    kk[0]     = int(*(w+iak));
    kk[1]     = int(*(w+iak+1));
    kk[2]     = int(*(w+iak+2));
    kk[3]     = int(*(w+iak+3));
    kk[4]     = iws_FingerPrint(w,ia);
}

//----------------------------------------------------------------------
int iP3(double *table, int i, int j, int k) {
    static int kk[5];
    int ia    = int(*(table+1));
    double *w = table-ia;
    int ifp   = iws_FingerPrint(w,ia);
    if(kk[4] != ifp) { K3(table, kk); }
    int ip = kk[0] + kk[1]*i + kk[2]*j + kk[3]*k;
    return ip;
}

//----------------------------------------------------------------------
double val3(double *table, int i, int j, int k) {
    int ip = iP3(table,i,j,k);
    return *(table+ip);
}

//----------------------------------------------------------------------
double val1(double *table, int i) {
    int ia    = int(*(table+1));
    double *w = table-ia;
    int iak   = iws_IaKARRAY(w,ia);
    int k0    = int(*(w+iak));
    return *(table+k0+i);
}

//----------------------------------------------------------------------
void fastloop(double *table, int i1, int i2, int j1,
                             int j2, int k1, int k2) {
    int kk[5];
    K3(table, kk);
    int inci = kk[1];
    int incj = kk[2];
    int inck = kk[3];
    int ia   = kk[0] + kk[1]*i1 + kk[2]*j1 + kk[3]*k1;
    for( int i = i1; i <= i2; i++) {
        int ja = ia;
        for(int j = j1; j <= j2; j++) {
            int ka = ja;
            for(int k = k1; k <= k2; k++) {
            cout << " T(" << i << "," << j << "," << k << ") = "
                 << *(table+ka) << endl;
                ka += inck;
            }
            ja += incj;
        }
        ia += inci;
    }
}

//----------------------------------------------------------------------
void wtable(double *w, int &ias, int &iaw) {
    int imi[] = {1,1,3};
    int ima[] = {50,50,6};
    if(NT < 1) {
        cout << "WTABLE: increase parameter NT to at least 1";
        exit(1);
    }
//  Create table-set
    ias     = iws_NewSet(w);
    iaw     = iws_WTable(w,imi,ima,3);
    int iat = iws_IaFirstTag(w,ias);
    w[iat]  = double(iaw-ias);
//  Fill weight table with delta-function
    for(int k = imi[2]; k <= ima[2]; k++) {
        for(int j = imi[1]; j <= ima[1]; j++) {
            w[iP3(w+iaw,j,j,k)+iaw] = 1;
        }
    }
}

//----------------------------------------------------------------------
void IniTab(double *table, double val) {
    int ia = int(*(table+1));
    double *w = table-ia;
    int i1 = iws_BeginTbody(w,ia);
    int i2 = iws_EndTbody(w,ia);
    for(int i = i1; i <= i2; i++) {
        *(w+i) = val;
    }
}

