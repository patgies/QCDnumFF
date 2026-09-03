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

namespace WSTORE {

  /**************************************************************/
  /*  WSTORE istore routines from istore.f                      */
  /**************************************************************/

  //--------------------------------------------------------------
    void sws_IwInit(int *iw, int nw, string txt)
    {
      int ls = txt.size();
      char *ctxt = new char[ls+1];
      strcpy(ctxt,txt.c_str());
      fsws_iwinitcpp(iw, &nw, ctxt, &ls);
      delete[] ctxt;
    }
  //--------------------------------------------------------------
    void sws_SetIwN(int *iw, int nw)
    {
     fsws_setiwn(iw, &nw);
    }
  //--------------------------------------------------------------
    int iws_Iarray(int* iw, int imi, int ima)
    {
     int ia = fiws_iarray(iw, &imi, &ima);
     return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    int iws_IAread(int* iw, int* iarr, int n)
    {
     int ia = fiws_iaread(iw, iarr, &n);
     return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    int iws_DAread(int* iw, double* darr, int n)
    {
     int ia = fiws_daread(iw, darr, &n);
     return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    void sws_IwWipe(int* obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     fsws_iwwipe(iw, &ia);
    }
  //--------------------------------------------------------------
    int iws_IhSize()
    {
     return fiws_ihsize();
    }
  //--------------------------------------------------------------
    int iws_IsaIstore(int *obj)
    {
     int  ia = *(obj+1);
     int *iw = obj-ia;
     return fiws_isaistore(iw);
    }
  //--------------------------------------------------------------
    int iws_IwSize(int *obj)
    {
     int  ia = *(obj+1);
     int *iw = obj-ia;
     return fiws_iwsize(iw);
    }
  //--------------------------------------------------------------
    int iws_IwNused(int *obj)
    {
     int  ia = *(obj+1);
     int *iw = obj-ia;
     return fiws_iwnused(iw);
    }
  //--------------------------------------------------------------
    int iws_IwNarrays(int *obj)
    {
     int  ia = *(obj+1);
     int *iw = obj-ia;
     return fiws_iwnarrays(iw);
    }
  //--------------------------------------------------------------
    int iws_IaLastObj(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = fiws_ialastobj(iw);
     return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    int iws_IwNheader(int *obj)
    {
     int  ia = *(obj+1);
     int *iw = obj-ia;
     return fiws_iwnheader(iw);
    }
  //--------------------------------------------------------------
    int iws_IwObjectType(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     return fiws_iwobjecttype(iw,&ia);
    }
  //--------------------------------------------------------------
    int iws_IwAsize(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     return fiws_iwasize(iw,&ia);
    }
  //--------------------------------------------------------------
    int iws_IwAnumber(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     return fiws_iwanumber(iw,&ia);
    }
  //--------------------------------------------------------------
    int iws_IwAfprint(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     return fiws_iwafprint(iw,&ia);
    }
  //--------------------------------------------------------------
    int iws_IwAdim(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     return fiws_iwadim(iw,&ia);
    }
  //--------------------------------------------------------------
    int iws_IwAimin(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     return fiws_iwaimin(iw,&ia);
    }
  //--------------------------------------------------------------
    int iws_IwAimax(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     return fiws_iwaimax(iw,&ia);
    }
  //--------------------------------------------------------------
    int iws_IaAbegin(int *obj)
    {
     int  ja   = *(obj+1);
     int *iw   = obj-ja;
     int  ia   = iaCtoF(ja);
     int iak   = fiws_iaabegin(iw, &ia);
     return iaFtoC(iak);
    }
  //--------------------------------------------------------------
    int iws_IaAend(int *obj)
    {
     int  ja   = *(obj+1);
     int *iw   = obj-ja;
     int  ia   = iaCtoF(ja);
     int iak   = fiws_iaaend(iw, &ia);
     return iaFtoC(iak);
    }
  //--------------------------------------------------------------
    int iws_IwKnul(int *obj)
    {
     int  ja = *(obj+1);
     int *iw = obj-ja;
     int  ia = iaCtoF(ja);
     return fiws_iwknul(iw,&ia);
    }
  //--------------------------------------------------------------
    int iws_ArrayI(int *obj, int i)
    {
     int  ja   = *(obj+1);
     int *iw   = obj-ja;
     int  ia   = iaCtoF(ja);
     int iak   = fiws_arrayi(iw, &ia, &i);
     return iaFtoC(iak);
    }
  //--------------------------------------------------------------
    void sws_IwTree(int *iw)
    {
     int iroot    = 0;
     fswsiwtree(iw,&iroot);
    }
  //--------------------------------------------------------------
    void sws_IwHead(int *obj)
    {
     int ja    = *(obj+1);
     int *iw   = obj-ja;
     int ia    = iaCtoF(ja);
     fsws_iwhead(iw,&ia);
    }

  /**************************************************************/
  /*  WSTORE workspace routines from wstore.f                   */
  /**************************************************************/

  //--------------------------------------------------------------
    int iws_Version()
    {
      return fiws_version();
    }
  //--------------------------------------------------------------
    int iws_WsInit(double *w, int nw, int nt, string txt)
    {
      int ls = txt.size();
      char *ctxt = new char[ls+1];
      strcpy(ctxt,txt.c_str());
      int ia = fiwswsinitcpp(w, &nw, &nt, ctxt, &ls);
      delete[] ctxt;
      return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    void sws_SetWsN(double *w, int nw)
    {
     fsws_setwsn(w, &nw);
    }
  //--------------------------------------------------------------
    void sws_Stampit(double *w)
    {
     fsws_stampit(w);
    }
  //--------------------------------------------------------------
    int iws_NewSet(double *w)
    {
     int ia = fiws_newset(w);
     return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    int iws_WTable(double* w, int* imi, int* ima, int n)
    {
     int ia = fiws_wtable(w, imi, ima, &n);
     return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    int iws_WClone(double *w1, int ia1, double *w2)
    {
     int ja1    = iaCtoF(ia1);
     int ia     = fiws_wclone(w1, &ja1, w2);
     return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    void sws_TbCopy(double *w1, int ia1,
                    double *w2, int ia2, int itag)
    {
     int ja1    = iaCtoF(ia1);
     int ja2    = iaCtoF(ia2);
     fsws_tbcopy(w1, &ja1, w2, &ja2, &itag);
    }
  //--------------------------------------------------------------
    void sws_WsWipe(double *w, int ia)
    {
     int ja     = iaCtoF(ia);
     fsws_wswipe(w, &ja);
    }
  //--------------------------------------------------------------
    void sws_TsDump(string fname, int key,
                    double *w, int ia, int &ierr)
    {
     int ls = fname.size();
     char *cfname = new char[ls+1];
     strcpy(cfname,fname.c_str());
     int ja     = iaCtoF(ia);
     fswstsdumpcpp(cfname, &ls, &key, w, &ja, &ierr);
    }
  //--------------------------------------------------------------
    int iws_TsRead(string fname, int key, double *w, int &ierr)
    {
     int ls = fname.size();
     char *cfname = new char[ls+1];
     strcpy(cfname,fname.c_str());
     int ia = fiwstsreadcpp(cfname, &ls, &key, w, &ierr);
     return iaFtoC(ia);
    }
  //--------------------------------------------------------------
    void sws_WsMark(int &mws, int &mset, int &mtab)
    {
     fsws_wsmark(&mws, &mset, &mtab);
    }
  //--------------------------------------------------------------
    int iws_HdSize()
    {
     return fiws_hdsize();
    }
  //--------------------------------------------------------------
    int iws_TbSize(int *imi, int *ima, int ndim)
    {
     return fiws_tbsize(imi, ima, &ndim);
    }
  //--------------------------------------------------------------
    int iws_Tpoint(double *w, int ia, int *index, int n)
    {
     int ja    = iaCtoF(ia);
     int iadr  = fiws_tpoint(w, &ja, index, &n);
     return iaFtoC(iadr);
    }
  //--------------------------------------------------------------
    int iws_IsaWorkspace(double *w)
    {
     return fiws_isaworkspace(w);
    }
  //--------------------------------------------------------------
    int iws_SizeOfW(double *w)
    {
     return fiws_sizeofw(w);
    }
  //--------------------------------------------------------------
    int iws_WordsUsed(double *w)
    {
     return fiws_wordsused(w);
    }
  //--------------------------------------------------------------
    int iws_Nheader(double *w)
    {
     return fiws_nheader(w);
    }
  //--------------------------------------------------------------
    int iws_Ntags(double *w)
    {
     return fiws_ntags(w);
    }
  //--------------------------------------------------------------
    int iws_HeadSkip(double *w)
    {
     return fiws_headskip(w);
    }
  //--------------------------------------------------------------
    int iws_IaRoot()
    {
     int ja    = fiws_iaroot();
     return iaFtoC(ja);
    }
  //--------------------------------------------------------------
    int iws_IaDrain(double *w)
    {
     int ja    = fiws_iadrain(w);
     return iaFtoC(ja);
    }
  //--------------------------------------------------------------
    int iws_IaNull(double *w)
    {
     int ja    = fiws_ianull(w);
     return iaFtoC(ja);
    }
  //--------------------------------------------------------------
    int iws_ObjectType(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_objecttype(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_ObjectSize(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_objectsize(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_Nobjects(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_nobjects(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_ObjectNumber(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_objectnumber(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_FingerPrint(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_fingerprint(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_IaFirstTag(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     int ift   = fiws_iafirsttag(w, &ja);
     return iaFtoC(ift);
    }
  //--------------------------------------------------------------
    int iws_TableDim(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_tabledim(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_IaKARRAY(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     int iak   = fiws_iakarray(w, &ja);
     return iaFtoC(iak);
    }
  //--------------------------------------------------------------
    int iws_IaIMIN(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     int iak   = fiws_iaimin(w, &ja);
     return iaFtoC(iak);
    }
  //--------------------------------------------------------------
    int iws_IaIMAX(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     int iak   = fiws_iaimax(w, &ja);
     return iaFtoC(iak);
    }
  //--------------------------------------------------------------
    int iws_BeginTbody(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     int iak   = fiws_begintbody(w, &ja);
     return iaFtoC(iak);
    }
  //--------------------------------------------------------------
    int iws_EndTbody(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     int iak   = fiws_endtbody(w, &ja);
     return iaFtoC(iak);
    }
  //--------------------------------------------------------------
    int iws_TFskip(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_tfskip(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_TBskip(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_tbskip(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_SFskip(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_sfskip(w, &ja);
    }
  //--------------------------------------------------------------
    int iws_SBskip(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     return fiws_sbskip(w, &ja);
    }
  //--------------------------------------------------------------
    void sws_WsTree(double *w)
    {
     int iroot    = 0;
     fswswstree(w,&iroot);
    }
  //--------------------------------------------------------------
    void sws_WsHead(double *w, int ia)
    {
     int ja    = iaCtoF(ia);
     fsws_wshead(w,&ja);
    }

}
