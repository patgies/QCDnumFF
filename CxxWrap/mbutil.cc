// C++ definition

#include "QCDNUM/QCDNUM.h"
//#include "QCDNUM/CHARACTER.h"
#include "QCDNUM/QCDNUMfw.h"

#include <sstream>
#include <cstring>
#include <string>
#include <iostream>
using namespace std;

#define STR_EXPAND(top) #top
#define STR(tok) STR_EXPAND(tok)

namespace MBUTIL {

  /**************************************************************/
  /*  MBUTIL comparison routines from compa.f                   */
  /**************************************************************/

//----------------------------------------------------------------
  int lmb_eq(double a, double b, double epsi)
  {
    return flmb_eq(&a, &b, &epsi);
  }
//----------------------------------------------------------------
  int lmb_ne(double a, double b, double epsi)
  {
    return flmb_ne(&a, &b, &epsi);
  }
//----------------------------------------------------------------
  int lmb_ge(double a, double b, double epsi)
  {
    return flmb_ge(&a, &b, &epsi);
  }
//----------------------------------------------------------------
  int lmb_le(double a, double b, double epsi)
  {
    return flmb_le(&a, &b, &epsi);
  }
//----------------------------------------------------------------
  int lmb_gt(double a, double b, double epsi)
  {
    return flmb_gt(&a, &b, &epsi);
  }
//----------------------------------------------------------------
  int lmb_lt(double a, double b, double epsi)
  {
    return flmb_lt(&a, &b, &epsi);
  }

  /**************************************************************/
  /*  MBUTIL hash routines from hash.f                          */
  /**************************************************************/

  //--------------------------------------------------------------
    int imb_ihash(int iseed, int *imsg, int n)
    {
      return fimb_ihash(&iseed, imsg, &n);
    }
  //--------------------------------------------------------------
    int imb_jhash(int iseed, double *dmsg, int n)
    {
      return fimb_jhash(&iseed, dmsg, &n);
    }
  //--------------------------------------------------------------
    int imb_dhash(int iseed, double *dmsg, int n)
    {
      return fimb_dhash(&iseed, dmsg, &n);
    }

  /**************************************************************/
  /*  MBUTIL character routines from mchar.f                    */
  /**************************************************************/

//----------------------------------------------------------------
  void smb_itoch(int in, string &chout, int &leng)
  {
    int ls = 20;
    char *mychar = new char[ls+1];
    fsmb_itochcpp(&in,mychar,&ls,&leng);
    chout = "";
    for(int i = 0; i < leng; i++) {
      chout = chout + mychar[i];
    }
    delete[] mychar;
  }
//----------------------------------------------------------------
  void smb_dtoch(double dd, int n, string &chout, int &leng)
  {
    int ls = 20;
    char *mychar = new char[ls+1];
    fsmb_dtochcpp(&dd,&n,mychar,&ls,&leng);
    chout = "";
    for(int i = 0; i < leng; i++) {
      chout = chout + mychar[i];
    }
    delete[] mychar;
  }
//----------------------------------------------------------------
  void smb_hcode(int ihash, string &hcode)
  {
    int ls = 15;
    char *mychar = new char[ls+1];
    fsmb_hcodecpp(&ihash,mychar,&ls);
    hcode = "";
    for(int i = 0; i < ls; i++) {
      hcode = hcode + mychar[i];
    }
    delete[] mychar;
  }

  /**************************************************************/
  /*  MBUTIL interpolation routines from polint.f               */
  /**************************************************************/

  //--------------------------------------------------------------
    void smb_spline(int n, double *x, double *y, double *b, double *c, double *d)
    {
      fsmb_spline(&n, x, y, b, c, d);
    }
  //--------------------------------------------------------------
    double dmb_seval(int n, double u, double *x, double *y, double *b, double *c, double *d)
    {
      return fdmb_seval(&n, &u, x, y, b, c, d);
    }

  /**************************************************************/
  /*  MBUTIL utility routines from utils.f                      */
  /**************************************************************/

  //--------------------------------------------------------------
    int imb_version()
    {
      return fimb_version();
    }
  //--------------------------------------------------------------
    double dmb_gauss(double (*f)(double*), double a, double b,
                     double &eps)
    {
      return fdmb_gauss(f, &a, &b, &eps);
    }
  //--------------------------------------------------------------
    double dmb_gaus1(double (*f)(double*), double a, double b)
    {
      return fdmb_gaus1(f, &a, &b);
    }
  //--------------------------------------------------------------
    double dmb_gaus2(double (*f)(double*), double a, double b)
    {
      return fdmb_gaus2(f, &a, &b);
    }
  //--------------------------------------------------------------
    double dmb_gaus3(double (*f)(double*), double a, double b)
    {
      return fdmb_gaus3(f, &a, &b);
    }
  //--------------------------------------------------------------
    double dmb_gaus4(double (*f)(double*), double a, double b)
    {
      return fdmb_gaus4(f, &a, &b);
    }

  /**************************************************************/
  /*  MBUTIL vector routines from vectors.f                     */
  /**************************************************************/

//----------------------------------------------------------------
  void smb_ifill(int *ia, int n, int ival)
  {
    fsmb_ifill(ia, &n, &ival);
  }
//----------------------------------------------------------------
  void smb_vfill(double *a, int n, double val)
  {
    fsmb_vfill(a, &n, &val);
  }
//----------------------------------------------------------------
  void smb_vmult(double *a, int n, double val)
  {
    fsmb_vmult(a, &n, &val);
  }
//----------------------------------------------------------------
  void smb_vcopy(double *a, double* b, int n)
  {
    fsmb_vcopy(a, b, &n);
  }
//----------------------------------------------------------------
  void smb_icopy(int *ia, int* ib, int n)
  {
    fsmb_icopy(ia, ib, &n);
  }
//----------------------------------------------------------------
  void smb_vitod(int *ia, double *b, int n)
  {
    fsmb_vitod(ia, b, &n);
  }
//----------------------------------------------------------------
  void smb_vdtoi(double *a, int *ib, int n)
  {
    fsmb_vdtoi(a, ib, &n);
  }
//----------------------------------------------------------------
  void smb_vaddv(double *a, double *b, double *c, int n)
  {
    fsmb_vaddv(a, b, c, &n);
  }
//----------------------------------------------------------------
  void smb_vminv(double *a, double *b, double *c, int n)
  {
    fsmb_vminv(a, b, c, &n);
  }
//----------------------------------------------------------------
  double dmb_vdotv(double *a, double *b, int n)
  {
    return fdmb_vdotv(a, b, &n);
  }
//----------------------------------------------------------------
  double dmb_vnorm(int m, double *a, int n)
  {
    return fdmb_vnorm(&m, a, &n);
  }
//----------------------------------------------------------------
  int lmb_vcomp(double *a, double *b, int n, double epsi)
  {
    return flmb_vcomp(a, b, &n, &epsi);
  }
//----------------------------------------------------------------
  double dmb_vpsum(double *a, double *w, int n)
  {
    return fdmb_vpsum(a, w, &n);
  }

}
