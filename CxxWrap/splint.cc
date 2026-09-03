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

namespace SPLINT {

  /**************************************************************/
  /*  SPLINT routines from usrsplint.f                          */
  /**************************************************************/

  //--------------------------------------------------------------
    int isp_spvers()
    {
      return fisp_spvers();
    }
  //--------------------------------------------------------------
    void ssp_spinit(int nuser)
    {
      fssp_spinit(&nuser);
    }
  //--------------------------------------------------------------
    void ssp_uwrite(int i, double val)
    {
      fssp_uwrite(&i, &val);
    }
  //--------------------------------------------------------------
    double dsp_uread(int i)
    {
      return fdsp_uread(&i);
    }
  //--------------------------------------------------------------
    int isp_sxmake(int istep)
    {
      return fisp_sxmake(&istep);
    }
  //--------------------------------------------------------------
    int isp_sxuser(double* xarr, int nx)
    {
      return fisp_sxuser(xarr, &nx);
    }
  //--------------------------------------------------------------
    void ssp_sxfill(int ias, double (*fun)(int*,int*,bool*),
                    int iq)
    {
      fssp_sxfill(&ias, fun, &iq);
    }
  //--------------------------------------------------------------
    void ssp_sxfpdf(int ias, int iset, double *def, int isel,
                    int iq)
    {
      fssp_sxfpdf(&ias, &iset, def, &isel, &iq);
    }
  //--------------------------------------------------------------
    void ssp_sxf123(int ias, int iset, double *def, int istf,
                    int iq)
    {
      fssp_sxf123(&ias, &iset, def, &istf, &iq);
    }
  //--------------------------------------------------------------
    int isp_sqmake(int istep)
    {
      return fisp_sqmake(&istep);
    }
  //--------------------------------------------------------------
    int isp_squser(double* qarr, int nq)
    {
      return fisp_squser(qarr, &nq);
    }
  //--------------------------------------------------------------
    void ssp_sqfill(int ias, double (*fun)(int*,int*,bool*),
                    int ix)
    {
      fssp_sqfill(&ias, fun, &ix);
    }
  //--------------------------------------------------------------
    void ssp_sqfpdf(int ias, int iset, double *def, int isel,
                    int ix)
    {
      fssp_sqfpdf(&ias, &iset, def, &isel, &ix);
    }
  //--------------------------------------------------------------
    void ssp_sqf123(int ias, int iset, double *def, int istf,
                    int ix)
    {
      fssp_sqf123(&ias, &iset, def, &istf, &ix);
    }
  //--------------------------------------------------------------
    double dsp_funs1(int ia, double z, int ichk)
    {
      return fdsp_funs1(&ia, &z, &ichk);
    }
  //--------------------------------------------------------------
    double dsp_ints1(int ia, double u1, double u2)
    {
      return fdsp_ints1(&ia, &u1, &u2);
    }
  //--------------------------------------------------------------
    int isp_s2make(int istepx, int istepq)
    {
      return fisp_s2make(&istepx, &istepq);
    }
  //--------------------------------------------------------------
    int isp_s2user(double *xarr, int nx, double *qarr, int nq)
    {
      return fisp_s2user(xarr, &nx, qarr, &nq);
    }
  //--------------------------------------------------------------
    void ssp_s2fill(int ias, double (*fun)(int*,int*,bool*),
                    double rs)
    {
      fssp_s2fill(&ias, fun, &rs);
    }
  //--------------------------------------------------------------
    void ssp_s2fpdf(int ias, int iset, double *def, int isel,
                    double rs)
    {
      fssp_s2fpdf(&ias, &iset, def, &isel, &rs);
    }
  //--------------------------------------------------------------
    void ssp_s2f123(int ias, int iset, double *def, int istf,
                    double rs)
    {
      fssp_s2f123(&ias, &iset, def, &istf, &rs);
    }
  //--------------------------------------------------------------
    double dsp_funs2(int ia, double x, double q, int ichk)
    {
      return fdsp_funs2(&ia, &x, &q, &ichk);
    }
  //--------------------------------------------------------------
    double dsp_ints2(int ia, double x1, double x2, double q1,
                     double q2, double rs, int np)
    {
      return fdsp_ints2(&ia, &x1, &x2, &q1, &q2, &rs, &np);
    }
  //--------------------------------------------------------------
    void ssp_extrapu(int ia, int n)
    {
      fssp_extrapu(&ia, &n);
    }
  //--------------------------------------------------------------
    void ssp_extrapv(int ia, int n)
    {
      fssp_extrapv(&ia, &n);
    }
  //--------------------------------------------------------------
    int isp_splinetype(int ia)
    {
      return fisp_splinetype(&ia);
    }
  //--------------------------------------------------------------
    void ssp_splims(int ia, int &nu, double &umi, double &uma,
                    int &nv, double &vmi, double &vma, int &nb)
    {
      fssp_splims(&ia, &nu, &umi, &uma, &nv, &vmi, &vma, &nb);
    }
  //--------------------------------------------------------------
    void ssp_unodes(int ia, double* array, int n, int &nus)
    {
      fssp_unodes(&ia, array, &n, &nus);
    }
  //--------------------------------------------------------------
    void ssp_vnodes(int ia, double* array, int n, int &nvs)
    {
      fssp_vnodes(&ia, array, &n, &nvs);
    }
  //--------------------------------------------------------------
    void ssp_nprint(int ia)
    {
      fssp_nprint(&ia);
    }
  //--------------------------------------------------------------
    double dsp_rscut(int ia)
    {
      return fdsp_rscut(&ia);
    }
  //--------------------------------------------------------------
    double dsp_rsmax(int ia, double rs)
    {
      return fdsp_rsmax(&ia, &rs);
    }
  //--------------------------------------------------------------
    void ssp_mprint()
    {
      fssp_mprint();
    }
  //--------------------------------------------------------------
    int isp_spsize(int ia)
    {
      return fisp_spsize(&ia);
    }
  //--------------------------------------------------------------
    void ssp_erase(int ia)
    {
      fssp_erase(&ia);
    }
  //--------------------------------------------------------------
    void ssp_spdump(int ia, string fname)
    {
      int ls = fname.size();
      char *cfname = new char[ls+1];
      strcpy(cfname,fname.c_str());
      fssp_spdumpcpp(&ia, cfname, &ls);
      delete[] cfname;
    }
  //--------------------------------------------------------------
    int isp_spread(string fname)
    {
      int ls = fname.size();
      char *cfname = new char[ls+1];
      strcpy(cfname,fname.c_str());
      int ia = fisp_spreadcpp(cfname, &ls);
      delete[] cfname;
      return ia;
    }
  //--------------------------------------------------------------
    void ssp_spsetval(int ia, int i, double val)
    {
      fssp_spsetval(&ia, &i, &val);
    }
  //--------------------------------------------------------------
    double dsp_spgetval(int ia, int i)
    {
      return fdsp_spgetval(&ia, &i);
    }

}
