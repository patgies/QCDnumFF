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
  /*  HQSTF routines                                            */
  /**************************************************************/

  //--------------------------------------------------------------
    int ihqvers()
    {
      return fihqvers();
    }
//----------------------------------------------------------------
  void hqwords(int &ntotal, int &nused)
  {
    fhqwords(&ntotal,&nused);
  }
//----------------------------------------------------------------
  void hqparms(double *qmass, double &a, double &b)
  {
    fhqparms(qmass,&a,&b);
  }
//----------------------------------------------------------------
  double hqqfrmu(double qmu2)
  {
    return fhqqfrmu(&qmu2);
  }
//----------------------------------------------------------------
  double hqmufrq(double Q2)
  {
    return fhqmufrq(&Q2);
  }
//----------------------------------------------------------------
  void hswitch(int iset)
  {
    fhswitch(&iset);
  }
//----------------------------------------------------------------
  void hqstfun(int istf, int icbt, double *def, double *x, double *Q2, double *f, int n, int ichk)
  {
   fhqstfun(&istf,&icbt,def,x,Q2,f,&n,&ichk);
  }
//----------------------------------------------------------------
  void hqfillw(int istf, double *qmass, double aq, double bq, int &nused)
  {
    fhqfillw(&istf,qmass,&aq,&bq,&nused);
  }
//----------------------------------------------------------------
  void hqdumpw(int lun, string fname)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fhqdumpwcpp(&lun,cfname,&ls);
    delete[] cfname;
  }
//----------------------------------------------------------------
  void hqreadw(int lun, string fname, int &nused, int &ierr)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fhqreadwcpp(&lun,cfname,&ls,&nused,&ierr);
    delete[] cfname;
  }

}
