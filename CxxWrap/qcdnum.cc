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
  /*  QCDNUM address function from usrGlobalId.f                */
  /**************************************************************/

//----------------------------------------------------------------
  int ipdftab(int iset, int id)
  {
    return fipdftab(&iset, &id);
  }

  /**************************************************************/
  /*  QCDNUM convolution engine from usrcvol.f                  */
  /**************************************************************/

//----------------------------------------------------------------
  void makewta(double *w, int jd,
               double (*afun)(double*,double*,int*),
               double (*achi)(double*))
  {
    fmakewta(w, &jd, afun, achi);
  }
//----------------------------------------------------------------
  void makewtb(double *w, int jd,
               double (*bfun)(double*,double*,int*),
               double (*achi)(double*),int ndel)
  {
    fmakewtb(w, &jd, bfun, achi, &ndel);
  }
//----------------------------------------------------------------
  void makewrs(double *w, int jd,
               double (*rfun)(double*,double*,int*),
               double (*sfun)(double*,double*,int*),
               double (*achi)(double*),int ndel)
  {
    fmakewrs(w, &jd, rfun, sfun, achi, &ndel);
  }
//----------------------------------------------------------------
  void makewtd(double *w, int jd,
               double (*dfun)(double*,double*,int*),
               double (*achi)(double*))
  {
    fmakewtd(w, &jd, dfun, achi);
  }
//----------------------------------------------------------------
  void makewtx(double *w, int jd)
  {
    fmakewtx(w, &jd);
  }
//----------------------------------------------------------------
  int idspfun(string pname, int iord, int jset)
  {
    int ls = pname.size();
    char *cpname = new char[ls+1];
    strcpy(cpname,pname.c_str());
    return fidspfuncpp(cpname, &ls, &iord, &jset);
    delete[] cpname;
  }
//----------------------------------------------------------------
  void scalewt(double *w, double c, int jd)
  {
    fscalewt(w, &c, &jd);
  }
//----------------------------------------------------------------
  void copywgt(double *w, int jd1, int jd2, int iadd)
  {
    fcopywgt(w, &jd1, &jd2, &iadd);
  }
//----------------------------------------------------------------
  void wcrossw(double *w, int jda, int jdb, int jdc, int iadd)
  {
    fwcrossw(w, &jda, &jdb, &jdc, &iadd);
  }
//----------------------------------------------------------------
  void wtimesf(double *w, double (*fun)(int*,int*),
               int jd1, int jd2, int iadd)
  {
    fwtimesf(w, fun, &jd1, &jd2, &iadd);
  }
//----------------------------------------------------------------
  double fcrossk(double *w, int idw, int jdum, int idf,
                 int ix, int iq)
  {
    return ffcrossk(w, &idw, &jdum, &idf, &ix, &iq);
  }
//----------------------------------------------------------------
  double fcrossf(double *w, int idw, int jdum, int ida, int idb,
                 int ix, int iq)
  {
    return ffcrossf(w, &idw, &jdum, &ida, &idb, &ix, &iq);
  }
//----------------------------------------------------------------
  void efromqq(double *qvec, double *evec, int nf)
  {
    fefromqq(qvec, evec, &nf);
  }
//----------------------------------------------------------------
  void qqfrome(double *evec, double *qvec, int nf)
  {
    fqqfrome(evec, qvec, &nf);
  }
//----------------------------------------------------------------
  void stfunxq(double (*fun)(int*,int*), double *x, double *q,
               double *f, int n, int jchk)
  {
    fstfunxq(fun, x, q, f, &n, &jchk);
  }

  /**************************************************************/
  /*  QCDNUM error message routines from usrerr.f               */
  /**************************************************************/

//----------------------------------------------------------------
  void setumsg(string name)
  {
    int ls = name.size();
    char *cname = new char[ls+1];
    strcpy(cname,name.c_str());
    fsetumsgcpp(cname, &ls);
    delete[] cname;
  }
//----------------------------------------------------------------
  void clrumsg()
  {
    fclrumsg();
  }

  /**************************************************************/
  /*  QCDNUM toolbox evolution from usrevnn.f                   */
  /**************************************************************/

//----------------------------------------------------------------
  void evfilla(double *w, int id, double (*func)(int*,int*,int*))
  {
    fevfilla(w, &id, func);
  }
//----------------------------------------------------------------
  double evgetaa(double *w, int id, int iq, int &nf, int &ithresh)
  {
    return fevgetaa(w, &id, &iq, &nf, &ithresh);
  }
//----------------------------------------------------------------
  void evdglap(double *w, int idw, int ida, int idf, double *start, int m, int n, int *iqlim, int &nf, double &eps)
  {
    fevdglap(w, &idw, &ida, &idf, start, &m, &n, iqlim, &nf, &eps);
  }
//----------------------------------------------------------------
  void cpyparw(double *w, double *array, int n, int kset)
  {
    fcpyparw(w, array, &n, &kset);
  }
//----------------------------------------------------------------
  void useparw(double *w, int kset)
  {
    fuseparw(w, &kset);
  }
//----------------------------------------------------------------
  int keyparw(double *w, int kset)
  {
    return fkeyparw(w, &kset);
  }
//----------------------------------------------------------------
  int keygrpw(double *w, int kset, int igroup)
  {
    return fkeygrpw(w, &kset, &igroup);
  }
//----------------------------------------------------------------
  double evpdfij(double *w, int id, int ix, int jq, int jchk)
  {
    return fevpdfij(w, &id, &ix, &jq, &jchk);
  }
//----------------------------------------------------------------
  void evplist(double *w, int id, double *x, double *q, double *f, int n, int jchk)
  {
    fevplist(w, &id, x, q, f, &n, &jchk);
  }
//----------------------------------------------------------------
  void evtable(double *w, int id, double *xx, int nx, double *qq, int nq, double *pdf, int jchk)
  {
    fevtable(w, &id, xx, &nx, qq, &nq, pdf, &jchk);
  }
//----------------------------------------------------------------
  void evpcopy(double *w, int id, double *def, int n, int jset)
  {
    fevpcopy(w, &id, def, &n, &jset);
  }

  /**************************************************************/
  /*  QCDNUM evolution routines from usrevol.f                  */
  /**************************************************************/

//----------------------------------------------------------------
  double rfromf(double fscale2)
  {
    return frfromf(&fscale2);
  }
//----------------------------------------------------------------
  double ffromr(double rscale2)
  {
    return fffromr(&rscale2);
  }
//----------------------------------------------------------------
  double asfunc(double r2, int &nf, int &ierr)
  {
    return fasfunc(&r2,&nf,&ierr);
  }
//----------------------------------------------------------------
  double altabn(int iset, int iq, int n, int &ierr)
  {
    return faltabn(&iset,&iq,&n,&ierr);
  }
//----------------------------------------------------------------
  void evolfg(int iset, double (*func) (int*, double*), double *def, int iq0, double &epsi)
  {
    fevolfg(&iset,func,def,&iq0,&epsi);
  }
//----------------------------------------------------------------
  void evsgns(int iset, double (*func)(int*, double*), int *isns, int n, int iq0, double &epsi)
  {
    fevsgns(&iset,func,isns,&n,&iq0,&epsi);
  }
//----------------------------------------------------------------
  void pdfcpy(int iset1, int iset2)
  {
    fpdfcpy(&iset1,&iset2);
  }
//----------------------------------------------------------------
  void extpdf(double (*func)(int*, double*, double*, bool*), int iset, int n, double offset, double &epsi)
  {
    fextpdf(func,&iset,&n,&offset,&epsi);
  }
//----------------------------------------------------------------
  void usrpdf(double (*func)(int*, double*, double*, bool*), int iset, int n, double offset, double &epsi)
  {
    fusrpdf(func,&iset,&n,&offset,&epsi);
  }
//----------------------------------------------------------------
  int nptabs(int iset)
  {
    return fnptabs(&iset);
  }
//----------------------------------------------------------------
  int ievtyp(int iset)
  {
    return fievtyp(&iset);
  }

  /**************************************************************/
  /*  QCDNUM fast convolution routines from usrfast.f           */
  /**************************************************************/

//----------------------------------------------------------------
  void idscope(double *w, int jset)
  {
    fidscope(w, &jset);
  }
//----------------------------------------------------------------
  void fastini(double *xlist, double *qlist, int n, int jchk)
  {
    ffastini(xlist, qlist, &n, &jchk);
  }
//----------------------------------------------------------------
  void fastclr(int ibuf)
  {
    ffastclr(&ibuf);
  }
//----------------------------------------------------------------
  void fastinp(double *w, int jdf, double *coef, int jbuf, int iadd)
  {
    ffastinp(w, &jdf, coef, &jbuf, &iadd);
  }
//----------------------------------------------------------------
  void fastepm(int jdum, int jdf, int jbuf)
  {
    ffastepm(&jdum, &jdf, &jbuf);
  }
//----------------------------------------------------------------
  void fastsns(int jset, double *qvec, int isel, int jbuf)
  {
    ffastsns(&jset, qvec, &isel, &jbuf);
  }
//----------------------------------------------------------------
  void fastsum(int jset, double *coef, int jbuf)
  {
    ffastsum(&jset, coef, &jbuf);
  }
//----------------------------------------------------------------
  void fastfxk(double *w, int *idwt, int ibuf1, int jbuf2)
  {
    ffastfxk(w, idwt, &ibuf1, &jbuf2);
  }
//----------------------------------------------------------------
  void fastfxf(double *w, int idx, int ibuf1, int ibuf2, int jbuf3)
  {
    ffastfxf(w, &idx, &ibuf1, &ibuf2, &jbuf3);
  }
//----------------------------------------------------------------
  void fastkin(int ibuf, double (*fun)(int*,int*,int*,int*))
  {
    ffastkin(&ibuf, fun);
  }
//----------------------------------------------------------------
  void fastcpy(int ibuf1, int ibuf2, int iadd)
  {
    ffastcpy(&ibuf1, &ibuf2, &iadd);
  }
//----------------------------------------------------------------
  void fastfxq(int ibuf, double *stf, int n)
  {
    ffastfxq(&ibuf, stf, &n);
  }

  /**************************************************************/
  /*  QCDNUM grid routines from usrgrd.f                        */
  /**************************************************************/

//----------------------------------------------------------------
  void gxmake(double *xmin, int *iwt, int n, int nxin, int &nxout, int iord)
  {
    fgxmake(xmin,iwt,&n,&nxin,&nxout,&iord);
  }
//----------------------------------------------------------------
  int ixfrmx(double x)
  {
    return fixfrmx(&x);
  }
//----------------------------------------------------------------
  int xxatix(double x, int ix)
  {
    return fxxatix(&x,&ix);
  }
//----------------------------------------------------------------
  double xfrmix(int ix)
  {
    return fxfrmix(&ix);
  }
//----------------------------------------------------------------
  void gxcopy(double *array, int n, int &nx)
  {
    fgxcopy(array,&n,&nx);
  }
//----------------------------------------------------------------
  void gqmake(double *qarr, double *wgt, int n, int nqin, int &nqout)
  {
    fgqmake(qarr,wgt,&n,&nqin,&nqout);
  }
//----------------------------------------------------------------
  int iqfrmq(double q2)
  {
    return fiqfrmq(&q2);
  }
//----------------------------------------------------------------
  int qqatiq(double q2, int iq)
  {
    return fqqatiq(&q2,&iq);
  }
//----------------------------------------------------------------
  double qfrmiq(int iq)
  {
    return fqfrmiq(&iq);
  }
//----------------------------------------------------------------
  void gqcopy(double *array, int n, int &nq)
  {
    fgqcopy(array,&n,&nq);
  }
//----------------------------------------------------------------
  void grpars(int &nx, double &xmi, double &xma, int &nq, double &qmi, double &qma, int &iord)
  {
    fgrpars(&nx,&xmi,&xma,&nq,&qmi,&qma,&iord);
  }
//---------------------------------------------------------------------
  void setlim( int ixmin, int iqmin, int iqmax, double dummy)
  {
    fsetlim(&ixmin,&iqmin,&iqmax,&dummy);
  }
//---------------------------------------------------------------------
  void getlim( int iset, double &xmin, double &qmin, double &qmax, double &dummy)
  {
    fgetlim(&iset,&xmin,&qmin,&qmax,&dummy);
  }

  /**************************************************************/
  /*  QCDNUM init routines from usrini.f                        */
  /**************************************************************/

//----------------------------------------------------------------
  void qcinit(int lun, string fname)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fqcinitcpp(&lun,cfname,&ls);
    delete[] cfname;
  }
//----------------------------------------------------------------
  void setlun(int lun, string fname)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fsetluncpp(&lun,cfname,&ls);
    delete[] cfname;
  }
//----------------------------------------------------------------
  int nxtlun(int lmin)
  {
    return fnxtlun(&lmin);
  }
//----------------------------------------------------------------
  void setval(string opt, double val)
  {
    int ls = opt.size();
    char *copt = new char[ls+1];
    strcpy(copt,opt.c_str());
    fsetvalcpp(copt,&ls,&val);
    delete[] copt;
  }
//----------------------------------------------------------------
  void getval(string opt, double &val)
  {
    int ls = opt.size();
    char *copt = new char[ls+1];
    strcpy(copt,opt.c_str());
    fgetvalcpp(copt,&ls,&val);
    delete[] copt;
  }
//----------------------------------------------------------------
  void setint(string opt, int ival)
  {
    int ls = opt.size();
    char *copt = new char[ls+1];
    strcpy(copt,opt.c_str());
    fsetintcpp(copt,&ls,&ival);
    delete[] copt;
  }
//----------------------------------------------------------------
  void getint(string opt, int &ival)
  {
    int ls = opt.size();
    char *copt = new char[ls+1];
    strcpy(copt,opt.c_str());
    fgetintcpp(copt,&ls,&ival);
    delete[] copt;
  }
//----------------------------------------------------------------
  void qstore(string opt, int ival, double &val)
  {
    int ls = opt.size();
    char *copt = new char[ls+1];
    strcpy(copt,opt.c_str());
    fqstorecpp(copt,&ls,&ival,&val);
    delete[] copt;
  }

  /**************************************************************/
  /*  QCDNUM parameter routines from usrparams.f                */
  /**************************************************************/

//----------------------------------------------------------------
  void setord(int iord)
  {
    fsetord(&iord);
  }
//----------------------------------------------------------------
  void getord(int &iord)
  {
    fgetord(&iord);
  }
//----------------------------------------------------------------
  void setalf(double as0, double r20)
  {
    fsetalf(&as0,&r20);
  }
//----------------------------------------------------------------
  void getalf(double &as0, double &r20)
  {
    fgetalf(&as0,&r20);
  }
//----------------------------------------------------------------
  void setcbt(int nfix, int iqc, int iqb, int iqt)
  {
    fsetcbt(&nfix,&iqc,&iqb,&iqt);
  }
//----------------------------------------------------------------
  void mixfns(int nfix, double r2c, double r2b, double r2t)
  {
    fmixfns(&nfix,&r2c,&r2b,&r2t);
  }
//----------------------------------------------------------------
  void getcbt(int &nfix, double &q2c, double &q2b, double &q2t)
  {
    fgetcbt(&nfix,&q2c,&q2b,&q2t);
  }
//----------------------------------------------------------------
  void setabr(double ar, double br)
  {
    fsetabr(&ar,&br);
  }
//----------------------------------------------------------------
  void getabr(double &ar, double &br)
  {
    fgetabr(&ar,&br);
  }
//----------------------------------------------------------------
  int nfrmiq(int iset, int iq, int &ithresh)
  {
    return fnfrmiq(&iset,&iq,&ithresh);
  }
//----------------------------------------------------------------
  int nflavs(int iq, int &ithresh)
  {
    return fnflavs(&iq,&ithresh);
  }
//----------------------------------------------------------------
  void cpypar(double *array, int n, int iset)
  {
    fcpypar(array,&n,&iset);
  }
//----------------------------------------------------------------
  void usepar(int iset)
  {
    fusepar(&iset);
  }
//----------------------------------------------------------------
  int keypar(int iset)
  {
    return fkeypar(&iset);
  }
//----------------------------------------------------------------
  int keygrp(int iset, int igrp)
  {
    return fkeygrp(&iset,&igrp);
  }
//----------------------------------------------------------------
  void pushcp()
  {
    fpushcp();
  }
//----------------------------------------------------------------
  void pullcp()
  {
    fpullcp();
  }

  /**************************************************************/
  /*  QCDNUM pdf output routines from usrpdf.f                  */
  /**************************************************************/

//----------------------------------------------------------------
  void allfxq(int iset, double x, double qmu2, double *pdfs, int n, int ichk)
  {
    fallfxq(&iset,&x,&qmu2,pdfs,&n,&ichk);
  }
//----------------------------------------------------------------
  void allfij(int iset, int ix, int iq, double *pdfs, int n, int ichk)
  {
    fallfij(&iset,&ix,&iq,pdfs,&n,&ichk);
  }
//----------------------------------------------------------------
  double bvalxq(int iset, int id, double x, double qmu2, int ichk)
  {
    return fbvalxq(&iset,&id,&x,&qmu2,&ichk);
  }
//----------------------------------------------------------------
  double bvalij(int iset, int id, int ix, int iq, int ichk)
  {
    return fbvalij(&iset,&id,&ix,&iq,&ichk);
  }
//----------------------------------------------------------------
  double fvalxq(int iset, int id, double x, double qmu2, int ichk)
  {
    return ffvalxq(&iset,&id,&x,&qmu2,&ichk);
  }
//----------------------------------------------------------------
  double fvalij(int iset, int id, int ix, int iq, int ichk)
  {
    return ffvalij(&iset,&id,&ix,&iq,&ichk);
  }
//----------------------------------------------------------------
  double sumfxq(int iset, double *c, int isel, double x, double qmu2, int ichk)
  {
    return fsumfxq(&iset,c,&isel,&x,&qmu2,&ichk);
  }
//----------------------------------------------------------------
  double sumfij(int iset, double *c, int isel, int ix, int iq, int ichk)
  {
    return fsumfij(&iset,c,&isel,&ix,&iq,&ichk);
  }
//----------------------------------------------------------------
  double fsplne(int iset, int id, double x, int iq)
  {
    return ffsplne(&iset,&id,&x,&iq);
  }
//----------------------------------------------------------------
  double splchk(int iset, int id, int iq)
  {
    return fsplchk(&iset,&id,&iq);
  }
//----------------------------------------------------------------
  void fflist(int iset, double *c, int m, double *x, double *q, double *f, int n, int ichk)
  {
    ffflist(&iset,c,&m,x,q,f,&n,&ichk);
  }
//----------------------------------------------------------------
  void fftabl(int iset, double *c, int isel, double *x, int nx, double *q, int nq, double *table, int m, int ichk)
  {
    ffftabl(&iset,c,&isel,x,&nx,q,&nq,table,&m,&ichk);
  }
//----------------------------------------------------------------
  void ftable(int iset, double *c, int m, double *x, int nx, double *q, int nq, double *table, int ichk)
  {
    fftable(&iset,c,&m,x,&nx,q,&nq,table,&ichk);
  }
//----------------------------------------------------------------
  void ffplot(string fname, double(* fun)(int*, double*, bool*), int m, double zmi, double zma, int n, string txt)
  {
    int ls1 = fname.size();
    char *cfname = new char[ls1+1];
    strcpy(cfname,fname.c_str());
    int ls2 = txt.size();
    char *ctxt = new char[ls2+1];
    strcpy(ctxt,txt.c_str());
    fffplotcpp(cfname,&ls1,fun,&m,&zmi,&zma,&n,ctxt,&ls2);
    delete[] cfname;
    delete[] ctxt;
  }
//----------------------------------------------------------------
  void fiplot(string fname, double(* fun)(int*, double*, bool*), int m, double* zval, int n, string txt)
  {
    int ls1 = fname.size();
    char *cfname = new char[ls1+1];
    strcpy(cfname,fname.c_str());
    int ls2 = txt.size();
    char *ctxt = new char[ls2+1];
    strcpy(ctxt,txt.c_str());
    ffiplotcpp(cfname,&ls1,fun,&m,zval,&n,ctxt,&ls2);
    delete[] cfname;
    delete[] ctxt;
  }

  /**************************************************************/
  /*  QCDNUM workspace routines from usrstore.f                 */
  /**************************************************************/

//----------------------------------------------------------------
  void maketab(double *w, int nw, int *jtypes, int npar,
               int newt, int &jset, int &nwords)
  {
    fmaketab(w, &nw, jtypes, &npar, &newt, &jset, &nwords);
  }
//----------------------------------------------------------------
  void setparw(double *w, int jset, double *par, int n)
  {
    fsetparw(w, &jset, par, &n);
  }
//----------------------------------------------------------------
  void getparw(double *w, int jset, double *par, int n)
  {
    fgetparw(w, &jset, par, &n);
  }
//----------------------------------------------------------------
  void dumptab(double *w, int jset, int lun,
               string fnam, string fkey)
  {
    int lf = fnam.size();
    int lk = fkey.size();
    char *cfnam = new char[lf+1];
    char *cfkey = new char[lk+1];
    strcpy(cfnam,fnam.c_str());
    strcpy(cfkey,fkey.c_str());
    fdumptabcpp(w, &jset, &lun, cfnam, &lf, cfkey, &lk);
    delete[] cfnam;
    delete[] cfkey;
  }
//----------------------------------------------------------------
  void readtab(double *w, int nw, int lun, string fnam,
               string fkey, int newt, int &jset,
               int &nwords, int &ierr)
  {
    int lf = fnam.size();
    int lk = fkey.size();
    char *cfnam = new char[lf+1];
    char *cfkey = new char[lk+1];
    strcpy(cfnam,fnam.c_str());
    strcpy(cfkey,fkey.c_str());
    freadtabcpp(w, &nw, &lun, cfnam, &lf, cfkey, &lk,
                &newt, &jset, &nwords, &ierr);
    delete[] cfnam;
    delete[] cfkey;
  }

  /**************************************************************/
  /*  QCDNUM weight routines from usrwgt.f                      */
  /**************************************************************/

//----------------------------------------------------------------
  void fillwt(int itype, int &idmin, int &idmax, int &nwds)
  {
    ffillwt(&itype,&idmin,&idmax,&nwds);
  }
//----------------------------------------------------------------
  void dmpwgt(int itype, int lun, string fname)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fdmpwgtcpp(&itype,&lun,cfname,&ls);
    delete[] cfname;
  }
//----------------------------------------------------------------
  void readwt(int lun, string fname, int &idmin, int &idmax, int &nwds, int &ierr)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    freadwtcpp(&lun,cfname,&ls,&idmin,&idmax,&nwds,&ierr);
    delete[] cfname;
  }
//----------------------------------------------------------------
  void wtfile(int itype, string fname)
  {
    int ls = fname.size();
    char *cfname = new char[ls+1];
    strcpy(cfname,fname.c_str());
    fwtfilecpp(&itype,cfname,&ls);
    delete[] cfname;
  }
//----------------------------------------------------------------
  void nwused(int &nwtot, int &nwuse, int &ndummy)
  {
    fnwused(&nwtot,&nwuse,&ndummy);
  }

}
