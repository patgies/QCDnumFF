#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <time.h>
#include "QCDNUM/QCDNUM.h"
using namespace std;

#include "../testjobs/inputpdfCxx.h"

//   coarse x - Q2 grid
 
double xxtab[] = { 1.0e-5, 2.0e-5, 5.0e-5, 1.0e-4, 2.0e-4, 5.0e-4,
           1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 1.0e-1,
           1.5e-1, 2.0e-1, 3.0e-1, 4.0e-1, 5.5e-1, 7.0e-1, 9.0e-1 };
int nxtab  = 20;
double xmi = xxtab[0];
double xma = xxtab[nxtab-1];

double qqtab[] = { 2.0e0, 2.7e0, 3.6e0, 5.0e0, 7.0e0, 1.0e1, 1.4e1,
		    2.0e1, 3.0e1, 5.0e1, 7.0e1, 1.0e2, 2.0e2, 5.0e2, 1.0e3,
		    3.0e3, 1.0e4, 4.0e4, 2.0e5, 1.0e6 };
int nqtab  = 20;
double qmi = qqtab[0];
double qma = qqtab[nqtab-1];

//==============================================
void mypdfs(double x, double qmu2, double *xpdf)
//==============================================
{
  QCDNUM::allfxq(1,x,qmu2,xpdf,0,1);
  return;
}

//========
int main()
//========
{
  /*
   *   Test the zmstf package
   */

  // qcdnum input
  double pdfin[] =
    // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
    //-6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
    {  0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,   //dval
       0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,   //uval
       0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,   //sval
       0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,   //dbar
       0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,   //ubar
       0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   //sbar
       0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,   //cval
       0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   //cbar
       0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,   //bval
       0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   //bbar
      -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,   //tval
       1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }; //tbar
  // Holds q-grid definition
  double qarr[] = {qmi, qma};
  double warr[] = {1., 1.};

  // Define output distributions
  //              tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
  //              -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
  double dnv[] = { 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0. };
  double upv[] = { 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0. };
  double del[] = { 0., 0., 0., 0.,-1., 1., 0., 0., 0., 0., 0., 0., 0. };
  double uds[] = { 0., 0., 0., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0. };
  double pro[] = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };
  double duv[] = { 0., 0., 0., 0.,-1.,-1., 0., 1., 1., 0., 0., 0., 0. };

  // cbt quark masses
  double hqmass[] = { 2., 5., 80. };

  // LO/NLO/NNLO
  int iord = 3;
  // VFNS
  int nfix = 0;
  // Define alpha_s
  double alf = 0.35e0;
  double q2a = 2.0e0;
  // More grid parameters
  int iosp  = 3;
  int n_x   = 100;
  int n_q   = 60;
  double q0 = 2.e0;

  int lun = 6;
  int lunout = abs(lun);
  string fname = " ";
  QCDNUM::qcinit(lun,fname);
  // QCDNUM calls to pass the values defined above
  QCDNUM::setord(iord);
  QCDNUM::setalf(alf,q2a);
  // x-mu2 grid
  int nxout;
  double xmin[] = {xmi};
  int iwt[] = {1};
  QCDNUM::gxmake(xmin,iwt,1,n_x,nxout,iosp);
  int nqout;
  QCDNUM::gqmake(qarr,warr,2,n_q,nqout);
  int iq0 = QCDNUM::iqfrmq(q0);
  // Thresholds
  int iqc = QCDNUM::iqfrmq(hqmass[0]*hqmass[0]);
  int iqb = QCDNUM::iqfrmq(hqmass[1]*hqmass[1]);
  int iqt = QCDNUM::iqfrmq(hqmass[2]*hqmass[2]);    
  QCDNUM::setcbt(nfix,iqc,iqb,iqt);

  // Try to read the weight file and create one if that fails
  string wname = "../weights/unpolarised.wgt";
  QCDNUM::wtfile(1,wname);

  int nztot, nzuse;
  QCDNUM::zmwords(nztot,nzuse);
  cout << " nztot, nzuse = " <<  nztot << " " << nzuse << endl;

  // Try to read the weight file and create one if that fails
  int lunw      = QCDNUM::nxtlun(10);
  string zmname = "../weights/zmstf.wgt";
  int nwords    = 0;
  int ierr      = 0;
  QCDNUM::zmreadw(lunw,zmname,nwords,ierr);
  if(ierr != 0) {
    QCDNUM::zmfillw(nwords);
    QCDNUM::zmdumpw(lunw,zmname);
  }

  // Scale variation
  //QCDNUM::setabr(2.,0.);
  //QCDNUM::zmdefq2(2.,1.); 

  // Evolve
  double epsi = 1e-3;
  int iter = 1;
  const clock_t tim1 = clock();
  for(int i=0; i<iter; i++) QCDNUM::evolfg(1,func,pdfin,iq0,epsi);
  const clock_t tim2 = clock();

  // Charge squared weighted for proton
  for(int i=0; i<13; i++) pro[i] /= 9;

  int ichk = 1;
  int iset = 1;

  cout << endl;
  cout << "     mu2       x           uv         dv        del        uds         gl" << endl;
  cout << setprecision(4) << scientific;

  for(int iq=0; iq<nqtab; iq=iq+4) {
    double q = qqtab[iq];
    cout << endl;
    for(int ix=0; ix<nxtab; ix=ix+4) {
      double x  = xxtab[ix];
      double uv = QCDNUM::sumfxq(iset,upv,1,x,q,ichk);
      double dv = QCDNUM::sumfxq(iset,dnv,1,x,q,ichk);
      double de = QCDNUM::sumfxq(iset,del,1,x,q,ichk);
      double ud = QCDNUM::sumfxq(iset,uds,1,x,q,ichk);
      double gl = QCDNUM::fvalxq(iset,  0,x,q,ichk);
      cout << q << " " << x << "  " << uv << " " << dv << " " << de << " " << ud << " " << gl << endl;
    }
  }

  cout << endl;
  cout << "     mu2       x           pr        F2         FL        FLp        xF3" << endl;
     
  QCDNUM::zswitch(iset);

  for(int iq=0; iq<nqtab; iq=iq+4) {
    double q[] = {qqtab[iq]};
    double qp  = qqtab[iq];
    cout << endl;
    for(int ix=0; ix<nxtab; ix=ix+4) {
      double x[] = {xxtab[ix]};
      double xp  = xxtab[ix];
      double pr = QCDNUM::sumfxq(iset,pro,1,xp,qp,ichk);
      // Ideally, zmstfun should not be called in a loop...
      double F2[1], FL[1], FLp[1], xF3[1];
      QCDNUM::zmstfun(2,pro,x,q,F2, 1,ichk);
      QCDNUM::zmstfun(1,pro,x,q,FL, 1,ichk);
      QCDNUM::zmstfun(4,pro,x,q,FLp,1,ichk);
      QCDNUM::zmstfun(3,duv,x,q,xF3,1,ichk);
      cout << q[0] << " " << x[0] << "  " << pr << " " << F2[0] << " " << FL[0] << " " << FLp[0] << " " << xF3[0] << endl;
    }
  }

  const clock_t tim3 = clock();

  double xlist[100], qlist[100], flist[100];
  // Make a list of interpolation points      
  int nlist = 0;
  for(int iq=0; iq<nqtab; iq=iq+4) {
    for(int ix=0; ix<nxtab; ix=ix+4) {
      xlist[nlist] = xxtab[ix];
      qlist[nlist] = qqtab[iq];
      nlist++;
    }
  }

  // ZmStfun with the list of points, instead of calling it in a loop
  // This is to see how much faster that is (factor 10!)  
  QCDNUM::zmstfun(2,pro,xlist,qlist,flist,nlist,ichk);
  QCDNUM::zmstfun(1,pro,xlist,qlist,flist,nlist,ichk);
  QCDNUM::zmstfun(4,pro,xlist,qlist,flist,nlist,ichk);
  QCDNUM::zmstfun(3,pro,xlist,qlist,flist,nlist,ichk);

  const clock_t tim4 = clock();

  // Print alfas
  cout << endl;
  cout << "    mu2       alfas" << endl;
  for(int iq=0; iq<nqtab; iq=iq+4) {
    double q = qqtab[iq];
    int nf;
    double a = QCDNUM::asfunc(q,nf,ierr);
    cout << q << " " << a << " " << ierr << endl;
  }

  cout << endl;
  cout << "time spent evol: " << float( tim2 - tim1 ) / iter / CLOCKS_PER_SEC << endl;
  cout << "time spent stfs: " << float( tim3 - tim2 ) / CLOCKS_PER_SEC << endl;
  cout << "time spent stfs: " << float( tim4 - tim3 ) / CLOCKS_PER_SEC << endl;

}
      

