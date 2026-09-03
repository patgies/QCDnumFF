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

namespace QCDNUM {

  /**************************************************************/
  /*  ZMSTF routines                                            */
  /**************************************************************/

  //--------------------------------------------------------------
    int izmvers()
    {
      return fizmvers();
    }
//----------------------------------------------------------------
  void zmwords(int &ntotal, int &nused)
  {
    fzmwords(&ntotal,&nused);
  }
//----------------------------------------------------------------
  void zmdefq2(double a, double b)
  {
    fzmdefq2(&a,&b);
  }
//----------------------------------------------------------------
  void zmabval(double &a, double &b)
  {
    fzmabval(&a,&b);
  }
//----------------------------------------------------------------
  double zmqfrmu(double qmu2)
  {
    return fzmqfrmu(&qmu2);
  }
//----------------------------------------------------------------
  double zmufrmq(double Q2)
  {
    return fzmufrmq(&Q2);
  }
//----------------------------------------------------------------
  void zswitch(int iset)
  {
    fzswitch(&iset);
  }
//----------------------------------------------------------------
  void zmstfun(int istf, double *def, double *x, double *Q2, double *f, int n, int ichk)
  {
    fzmstfun(&istf,def,x,Q2,f,&n,&ichk);
  }
//----------------------------------------------------------------
  void zmfillw(int &nused)
  {
    fzmfillw(&nused);
  }
//----------------------------------------------------------------
  void zmdumpw(int lun, string fname)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fzmdumpwcpp(&lun,cfname,&ls);
    delete[] cfname;
  }
//----------------------------------------------------------------
  void zmreadw(int lun, string fname, int &nused, int &ierr)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fzmreadwcpp(&lun,cfname,&ls,&nused,&ierr);
    delete[] cfname;
  }
//----------------------------------------------------------------
  void zmwfile(string fname)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fzmwfilecpp(cfname,&ls);
    delete[] cfname;
  }

}
