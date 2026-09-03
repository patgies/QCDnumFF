/*
 * ---------------------------------------------------------------------
 * Time-like (fragmentation function) evolution of the Kniehl-Kramer
 * D0 fragmentation function, LO, mu0=mQ convention.
 * B.A. Kniehl and G. Kramer, "Charmed-Hadron Fragmentation Functions
 * from CERN LEP1 Revisited", DESY 06-102 (arXiv:hep-ph/0607306).
 *
 * The paper has TWO separate non-perturbative inputs:
 *   - charm-quark FF (Peterson form) at mu0 = mc
 *   - bottom-quark FF (power law)    at mu0 = mb
 * with everything else zero at that flavour's own threshold.
 *
 * Because DGLAP evolution is linear, the full result is obtained
 * by evolving each input separately on the same grid and adding:
 *   D_total(x,Q) = D_c-evolved(x,Q) + D_b-evolved(x,Q)   for Q >= mb
 *   D_total(x,Q) = D_c-evolved(x,Q)                      for mc <= Q < mb
 * ---------------------------------------------------------------------
 */

#include <cmath>
#include <cstdio>
#include "QCDNUM/QCDNUM.h"

using namespace std;

//----------------------------------------------------------------------
// Charm-quark Peterson input D_c(x,mc^2), D0 meson, LO
// (Table I of hep-ph/0607306: N=0.694, epsilon=0.101)
double funcC(int* ipdf, double* xx) {
  int    i = *ipdf;
  double x = *xx;
  double f = 0;
  if (i == 8) {                              // cbar column -> c=cbar=D_c
    double N = 0.694, eps = 0.101;
    double den = (1-x)*(1-x) + eps*x;
    f = N * x*x * (1-x)*(1-x) / (den*den);   // x * D_c(x,mc^2)
  }
  return f;
}
//----------------------------------------------------------------------
// Bottom-quark power-law input D_b(x,mb^2), D0 meson, LO
// (Table I of hep-ph/0607306: N=81.7, alpha=1.81, beta=4.95)
double funcB(int* ipdf, double* xx) {
  int    i = *ipdf;
  double x = *xx;
  double f = 0;
  if (i == 10) {                             // bbar column -> b=bbar=D_b
    double N = 81.7, alfa = 1.81, beta = 4.95;
    f = N * pow(x, alfa+1) * pow(1-x, beta); // x * D_b(x,mb^2)
  }
  return f;
}
//----------------------------------------------------------------------
int main() {

  // Standard QCDNum flavour basis (shared by both evolutions):
  // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
  double def[] =
    { 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,      // dval
      0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,      // uval
      0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,      // sval
      0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,      // dbar
      0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,      // ubar
      0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // sbar
      0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,      // cval
      0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // cbar
      0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,      // bval
      0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // bbar
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // tval (zero)
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};     // tbar (zero)

  double qmc2 = 2.25, qmb2 = 25.0;                // mc=1.5, mb=5 GeV
  double as0  = 0.118, r20 = 8315.25;             // alphas(Mz2)=0.118, LO ref
  double xmin[] = {1.e-5};
  int    iwt[] = {1}, ng = 1, nxin = 300, iosp = 3;   // x grid, splord
  int    nqin = 100;
  // mu2 grid starts below qmc2 so that setcbt has grid points below
  // the charm threshold (required by QCDNUM)
  double qq[] = {1.0, 2.25, 25.0, 1.e6}, wt[] = {1.0, 1.0, 2.0, 1.0};

  int lun = 6; string outfile = " ";
  QCDNUM::qcinit(lun, outfile);
  int nx, nq;
  QCDNUM::gxmake(xmin, iwt, ng, nxin, nx, iosp);
  QCDNUM::gqmake(qq, wt, 4, nqin, nq);
  QCDNUM::wtfile(3, "../weights/timelike.wgt");       // itype=3: time-like wgts
  QCDNUM::setord(1);                                  // LO evolution
  QCDNUM::setalf(as0, r20);
  int iqc = QCDNUM::iqfrmq(qmc2);
  int iqb = QCDNUM::iqfrmq(qmb2);
  QCDNUM::setcbt(0, iqc, iqb, 999);                   // dynamic VFNS thresholds

  // Evolution A: charm-quark input at mu0 = mc -------------------------
  double epsC;
  int iq0C = QCDNUM::iqfrmq(qmc2);
  QCDNUM::evolfg(13, funcC, def, iq0C, epsC);   // itype=3(t-like),iset=1 -> jset=13

  // Evolution B: bottom-quark input at mu0 = mb -----------------------
  double epsB;
  int iq0B = QCDNUM::iqfrmq(qmb2);
  QCDNUM::evolfg(23, funcB, def, iq0B, epsB);   // itype=3(t-like),iset=2 -> jset=23

  // Get results  --------------------------------------------------------
  double x  = 0.1;
  double q2 = 100.0;                                  // an example scale > qmb2

  double pdfC[13], pdfB[13];
  QCDNUM::allfxq(1, x, q2, pdfC, 0, 1);

  double xDtot;
  if (q2 >= qmb2) {
    QCDNUM::allfxq(2, x, q2, pdfB, 0, 1);
    xDtot = pdfC[10] + pdfB[10];      // c-quark: Fortran index +4 -> C++ 10
  } else {
    xDtot = pdfC[10];
  }

  printf(" x, Q2, x*D_c^D0(x,Q2) = %13.5e %13.5e %14.5e\n", x, q2, xDtot);
  printf(" -> D_c^D0(x,Q2) (divide by x) = %14.5e\n", xDtot/x);

  return 0;
}
