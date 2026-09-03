#ifndef QCDNUMfw_H
#define QCDNUMfw_H

// Declarations for the Fortran/C interface
#include "QCDNUM/FortranWrappers.h"

extern "C" {


  /**************************************************************/
  /*  MBUTIL comparison routines from compa.f                   */
  /**************************************************************/

  #define flmb_eq FC_FUNC(lmb_eq,LMB_EQ)
    int flmb_eq(double*,double*,double*);
  #define flmb_ne FC_FUNC(lmb_ne,LMB_NE)
    int flmb_ne(double*,double*,double*);
  #define flmb_ge FC_FUNC(lmb_ge,LMB_GE)
    int flmb_ge(double*,double*,double*);
  #define flmb_le FC_FUNC(lmb_le,LMB_LE)
    int flmb_le(double*,double*,double*);
  #define flmb_gt FC_FUNC(lmb_gt,LMB_GT)
    int flmb_gt(double*,double*,double*);
  #define flmb_lt FC_FUNC(lmb_lt,LMB_LT)
    int flmb_lt(double*,double*,double*);

  /**************************************************************/
  /*  MBUTIL hash routines from hash.f                          */
  /**************************************************************/

  #define fimb_ihash FC_FUNC(imb_ihash,IMB_IHASH)
    int fimb_ihash(int*, int*, int*);
  #define fimb_jhash FC_FUNC(imb_jhash,IMB_JHASH)
    int fimb_jhash(int*, double*, int*);
  #define fimb_dhash FC_FUNC(imb_dhash,IMB_DHASH)
    int fimb_dhash(int*, double*, int*);

  /**************************************************************/
  /*  MBUTIL character routines from mchar.f                    */
  /**************************************************************/

  #define fsmb_itochcpp FC_FUNC(smb_itochcpp,SMB_ITOCHCPP)
    void fsmb_itochcpp(int*, char*, int*, int*);
  #define fsmb_dtochcpp FC_FUNC(smb_dtochcpp,SMB_DTOCHCPP)
    void fsmb_dtochcpp(double*, int*, char*, int*, int*);
  #define fsmb_hcodecpp FC_FUNC(smb_hcodecpp,SMB_HCODECPP)
    void fsmb_hcodecpp(int*, char*, int*);

  /**************************************************************/
  /*  MBUTIL interpolation routines from polint.f               */
  /**************************************************************/

  #define fsmb_spline FC_FUNC(smb_spline,SMB_SPLINE)
    void fsmb_spline(int*, double*, double*, double*, double*, double*);
  #define fdmb_seval FC_FUNC(dmb_seval,DMB_SEVAL)
    double fdmb_seval(int*, double*, double*, double*, double*, double*, double*);

  /**************************************************************/
  /*  MBUTIL utility routines from utils.f                      */
  /**************************************************************/

  #define fimb_version FC_FUNC(imb_version,IMB_VERSION)
    int fimb_version();
  #define fdmb_gauss FC_FUNC(dmb_gauss,DMB_GAUSS)
    int fdmb_gauss(double(*)(double*),double*,double*,double*);
  #define fdmb_gaus1 FC_FUNC(dmb_gaus1,DMB_GAUS1)
    int fdmb_gaus1(double(*)(double*),double*,double*);
  #define fdmb_gaus2 FC_FUNC(dmb_gaus2,DMB_GAUS2)
    int fdmb_gaus2(double(*)(double*),double*,double*);
  #define fdmb_gaus3 FC_FUNC(dmb_gaus3,DMB_GAUS3)
    int fdmb_gaus3(double(*)(double*),double*,double*);
  #define fdmb_gaus4 FC_FUNC(dmb_gaus4,DMB_GAUS4)
    int fdmb_gaus4(double(*)(double*),double*,double*);

  /**************************************************************/
  /*  MBUTIL vector routines from vectors.f                     */
  /**************************************************************/

  #define fsmb_ifill FC_FUNC(smb_ifill,SMB_IFILL)
    void fsmb_ifill(int*,int*,int*);
  #define fsmb_vfill FC_FUNC(smb_vfill,SMB_VFILL)
    void fsmb_vfill(double*,int*,double*);
  #define fsmb_vmult FC_FUNC(smb_vmult,SMB_VMULT)
    void fsmb_vmult(double*,int*,double*);
  #define fsmb_vcopy FC_FUNC(smb_vcopy,SMB_VCOPY)
    void fsmb_vcopy(double*,double*,int*);
  #define fsmb_icopy FC_FUNC(smb_icopy,SMB_ICOPY)
    void fsmb_icopy(int*,int*,int*);
  #define fsmb_vitod FC_FUNC(smb_vitod,SMB_VITOD)
    void fsmb_vitod(int*,double*,int*);
  #define fsmb_vdtoi FC_FUNC(smb_vdtoi,SMB_VDTOI)
    void fsmb_vdtoi(double*,int*,int*);
  #define fsmb_vaddv FC_FUNC(smb_vaddv,SMB_VADDV)
    void fsmb_vaddv(double*,double*,double*,int*);
  #define fsmb_vminv FC_FUNC(smb_vminv,SMB_VMINV)
    void fsmb_vminv(double*,double*,double*,int*);
  #define fdmb_vdotv FC_FUNC(dmb_vdotv,DMB_VDOTV)
    double fdmb_vdotv(double*,double*,int*);
  #define fdmb_vnorm FC_FUNC(dmb_vnorm,DMB_VNORM)
    double fdmb_vnorm(int*,double*,int*);
  #define flmb_vcomp FC_FUNC(lmb_vcomp,LMB_VCOMP)
    int flmb_vcomp(double*,double*,int*,double*);
  #define fdmb_vpsum FC_FUNC(dmb_vpsum,DMB_VPSUM)
    double fdmb_vpsum(double*,double*,int*);

  /**************************************************************/
  /*  WSTORE istore routines from istore.f                      */
  /**************************************************************/

  #define fsws_iwinitcpp FC_FUNC(sws_iwinitcpp,SWS_IWINITCPP)
    void fsws_iwinitcpp(int*, int*, char*, int*);
  #define fsws_setiwn FC_FUNC(sws_setiwn,SWS_SETIWN)
    void fsws_setiwn(int*, int*);
  #define fiws_iarray FC_FUNC(iws_iarray,IWS_IARRAY)
    int fiws_iarray(int*, int*, int*);
  #define fiws_iaread FC_FUNC(iws_iaread,IWS_IAREAD)
    int fiws_iaread(int*, int*, int*);
  #define fiws_daread FC_FUNC(iws_daread,IWS_DAREAD)
    int fiws_daread(int*, double*, int*);
  #define fsws_iwwipe FC_FUNC(sws_iwwipe,SWS_IWWIPE)
    void fsws_iwwipe(int*, int*);
  #define fiws_ihsize FC_FUNC(iws_ihsize,IWS_IHSIZE)
    int fiws_ihsize();

  #define fiws_isaistore FC_FUNC(iws_isaistore,IWS_ISAISTORE)
    int fiws_isaistore(int*);
  #define fiws_iwsize FC_FUNC(iws_iwsize,IWS_IWSIZE)
    int fiws_iwsize(int*);
  #define fiws_iwnused FC_FUNC(iws_iwnused,IWS_IWNUSED)
    int fiws_iwnused(int*);
  #define fiws_iwnarrays FC_FUNC(iws_iwnarrays,IWS_IWNARRAYS)
    int fiws_iwnarrays(int*);
  #define fiws_ialastobj FC_FUNC(iws_ialastobj,IWS_IALASTOBJ)
    int fiws_ialastobj(int*);
  #define fiws_iwnheader FC_FUNC(iws_iwnheader,IWS_IWNHEADER)
    int fiws_iwnheader(int*);
  #define fiws_iwobjecttype FC_FUNC(iws_iwobjecttype,IWS_IWOBJECTTYPE)
    int fiws_iwobjecttype(int*,int*);
  #define fiws_iwasize FC_FUNC(iws_iwasize,IWS_IWASIZE)
    int fiws_iwasize(int*,int*);
  #define fiws_iwanumber FC_FUNC(iws_iwanumber,IWS_IWANUMBER)
    int fiws_iwanumber(int*,int*);
  #define fiws_iwafprint FC_FUNC(iws_iwafprint,IWS_IWAFPRINT)
    int fiws_iwafprint(int*,int*);
  #define fiws_iwadim FC_FUNC(iws_iwadim,IWS_IWADIM)
    int fiws_iwadim(int*,int*);
  #define fiws_iwaimin FC_FUNC(iws_iwaimin,IWS_IWAIMIN)
    int fiws_iwaimin(int*,int*);
  #define fiws_iwaimax FC_FUNC(iws_iwaimax,IWS_IWAIMAX)
    int fiws_iwaimax(int*,int*);
  #define fiws_iaabegin FC_FUNC(iws_iaabegin ,IWS_IAABEGIN)
    int fiws_iaabegin(int*, int*);
  #define fiws_iaaend FC_FUNC(iws_iaaend ,IWS_IAAEND)
    int fiws_iaaend(int*, int*);
  #define fiws_iwknul FC_FUNC(iws_iwknul,IWS_IWKNUL)
    int fiws_iwknul(int*,int*);
  #define fiws_arrayi FC_FUNC(iws_arrayi ,IWS_ARRAYI)
    int fiws_arrayi(int*, int*, int*);
  #define fswsiwtree FC_FUNC(swsiwtree ,SWSIWTREE)
    int fswsiwtree(int*, int*);
  #define fsws_iwhead FC_FUNC(sws_iwhead ,SWS_IWHEAD)
    int fsws_iwhead(int*, int*);

  /**************************************************************/
  /*  WSTORE workspace routines from wstore.f                   */
  /**************************************************************/

  #define fiws_version FC_FUNC(iws_version,IWS_VERSION)
    int fiws_version();
  #define fiwswsinitcpp FC_FUNC(iwswsinitcpp,IWS_WSINITCPP)
    int fiwswsinitcpp(double*, int*, int*, char*, int*);
  #define fsws_setwsn FC_FUNC(sws_setwsn,SWS_SETWSN)
    void fsws_setwsn(double*, int*);
  #define fsws_stampit FC_FUNC(sws_stampit,SWS_STAMPIT)
    void fsws_stampit(double*);
  #define fiws_newset FC_FUNC(iws_newset,IWS_NEWSET)
    int fiws_newset(double*);
  #define fiws_wtable FC_FUNC(iws_wtable,IWS_WTABLE)
    int fiws_wtable(double*, int*, int*, int*);
  #define fiws_wclone FC_FUNC(iws_wclone,IWS_WCLONE)
    int fiws_wclone(double*, int*, double*);
  #define fsws_tbcopy FC_FUNC(sws_tbcopy,SWS_TBCOPY)
    void fsws_tbcopy(double*, int*, double*, int*, int*);
  #define fsws_wswipe FC_FUNC(sws_wswipe,SWS_WSWIPE)
    void fsws_wswipe(double*, int*);
  #define fswstsdumpcpp FC_FUNC(swstsdumpcpp,SWSTSDUMPCPP)
    void fswstsdumpcpp(char*, int*, int*, double*, int*, int*);
  #define fiwstsreadcpp FC_FUNC(iwstsreadcpp,IWSTSREADCPP)
    int fiwstsreadcpp(char*, int*, int*, double*, int*);
  #define fsws_wsmark FC_FUNC(sws_wsmark,SWS_WSMARK)
    void fsws_wsmark(int*, int*, int*);
  #define fiws_hdsize FC_FUNC(iws_hdsize,IWS_HDSIZE)
    int fiws_hdsize();
  #define fiws_tbsize FC_FUNC(iws_tbsize,IWS_TBSIZE)
    int fiws_tbsize(int*, int*, int*);
  #define fiws_tpoint FC_FUNC(iws_tpoint,IWS_TPOINT)
    int fiws_tpoint(double*, int*, int*, int*);

  #define fiws_isaworkspace FC_FUNC(iws_isaworkspace,IWS_ISAWORKSPACE)
    int fiws_isaworkspace(double*);
  #define fiws_sizeofw FC_FUNC(iws_sizeofw,IWS_SIZEOFW)
    int fiws_sizeofw(double*);
  #define fiws_wordsused FC_FUNC(iws_wordsused,IWS_WORDSUSED)
    int fiws_wordsused(double*);
  #define fiws_nheader FC_FUNC(iws_nheader,IWS_NHEADER)
    int fiws_nheader(double*);
  #define fiws_ntags FC_FUNC(iws_ntags,IWS_NTAGS)
    int fiws_ntags(double*);
  #define fiws_headskip FC_FUNC(iws_headskip,IWS_HEADSKIP)
    int fiws_headskip(double*);
  #define fiws_iaroot FC_FUNC(iws_iaroot,IWS_IAROOT)
    int fiws_iaroot();
  #define fiws_iadrain FC_FUNC(iws_iadrain,IWS_IADRAIN)
    int fiws_iadrain(double*);
  #define fiws_ianull FC_FUNC(iws_ianull,IWS_IANULL)
    int fiws_ianull(double*);
  #define fiws_objecttype FC_FUNC(iws_objecttype,IWS_OBJECTTYPE)
    int fiws_objecttype(double*, int*);
  #define fiws_objectsize FC_FUNC(iws_objectsize,IWS_OBJECTSIZE)
    int fiws_objectsize(double*, int*);
  #define fiws_nobjects FC_FUNC(iws_nobjects,IWS_NOBJECTS)
    int fiws_nobjects(double*, int*);
  #define fiws_objectnumber FC_FUNC(iws_objectnumber,IWS_OBJECTNUMBER)
    int fiws_objectnumber(double*, int*);
  #define fiws_fingerprint FC_FUNC(iws_fingerprint,IWS_FINGERPRINT)
    int fiws_fingerprint(double*, int*);
  #define fiws_iafirsttag FC_FUNC(iws_iafirsttag,IWS_IAFIRSTTAG)
    int fiws_iafirsttag(double*, int*);
  #define fiws_tabledim FC_FUNC(iws_tabledim,IWS_TABLEDIM)
    int fiws_tabledim(double*, int*);
  #define fiws_iakarray FC_FUNC(iws_iakarray,IWS_IAKARRAY)
    int fiws_iakarray(double*, int*);
  #define fiws_iaimin FC_FUNC(iws_iaimin ,IWS_IAIMIN)
    int fiws_iaimin(double*, int*);
  #define fiws_iaimax FC_FUNC(iws_iaimax ,IWS_IAIMAX)
    int fiws_iaimax(double*, int*);
  #define fiws_begintbody FC_FUNC(iws_begintbody ,IWS_BEGINTBODY)
    int fiws_begintbody(double*, int*);
  #define fiws_endtbody FC_FUNC(iws_endtbody ,IWS_ENDTBODY)
    int fiws_endtbody(double*, int*);
  #define fiws_tfskip FC_FUNC(iws_tfskip ,IWS_TFSKIP)
    int fiws_tfskip(double*, int*);
  #define fiws_tbskip FC_FUNC(iws_tbskip ,IWS_TBSKIP)
    int fiws_tbskip(double*, int*);
  #define fiws_sfskip FC_FUNC(iws_sfskip ,IWS_SFSKIP)
    int fiws_sfskip(double*, int*);
  #define fiws_sbskip FC_FUNC(iws_sbskip ,IWS_SBSKIP)
    int fiws_sbskip(double*, int*);
  #define fswswstree FC_FUNC(swswstree ,SMBWSTREE)
    int fswswstree(double*, int*);
  #define fsws_wshead FC_FUNC(sws_wshead ,SWS_WSHEAD)
    int fsws_wshead(double*, int*);

  /**************************************************************/
  /*  SPLINT routines from usrsplint.f                          */
  /**************************************************************/

  #define fisp_spvers FC_FUNC(isp_spvers,ISP_SPVERS)
    int fisp_spvers();
  #define fssp_spinit FC_FUNC(ssp_spinit,SSP_SPINIT)
    void fssp_spinit(int*);
  #define fssp_uwrite FC_FUNC(ssp_uwrite,SSP_UWRITE)
    void fssp_uwrite(int*, double*);
  #define fdsp_uread FC_FUNC(dsp_uread,DSP_UREAD)
    double fdsp_uread(int*);
  #define fisp_sxmake FC_FUNC(isp_sxmake,ISP_SXMAKE)
    int fisp_sxmake(int*);
  #define fisp_sxuser FC_FUNC(isp_sxuser,ISP_SXUSER)
    int fisp_sxuser(double *, int*);
  #define fssp_sxfill FC_FUNC(ssp_sxfill,SSP_SXFILL)
    void fssp_sxfill(int*, double (*)(int*,int*,bool*), int*);
  #define fssp_sxfpdf FC_FUNC(ssp_sxfpdf,SSP_SXFPDF)
    void fssp_sxfpdf(int*, int*, double*, int*, int*);
  #define fssp_sxf123 FC_FUNC(ssp_sxf123,SSP_SXF123)
    void fssp_sxf123(int*, int*, double*, int*, int*);
  #define fisp_sqmake FC_FUNC(isp_sqmake,ISP_SQMAKE)
    int fisp_sqmake(int*);
  #define fisp_squser FC_FUNC(isp_squser,ISP_SQUSER)
    int fisp_squser(double *, int*);
  #define fssp_sqfill FC_FUNC(ssp_sqfill,SSP_SQFILL)
    void fssp_sqfill(int*, double (*)(int*,int*,bool*), int*);
  #define fssp_sqfpdf FC_FUNC(ssp_sqfpdf,SSP_SQFPDF)
    void fssp_sqfpdf(int*, int*, double*, int*, int*);
  #define fssp_sqf123 FC_FUNC(ssp_sqf123,SSP_SQF123)
    void fssp_sqf123(int*, int*, double*, int*, int*);
  #define fdsp_funs1 FC_FUNC(dsp_funs1,DSP_FUNS1)
    double fdsp_funs1(int*,double*,int*);
  #define fdsp_ints1 FC_FUNC(dsp_ints1,DSP_INTS1)
    double fdsp_ints1(int*,double*,double*);
  #define fisp_s2make FC_FUNC(isp_s2make,ISP_S2MAKE)
    int fisp_s2make(int*, int*);
  #define fisp_s2user FC_FUNC(isp_s2user,ISP_s2user)
    int fisp_s2user(double*, int*, double*, int*);
  #define fssp_s2fill FC_FUNC(ssp_s2fill,SSP_S2FILL)
    void fssp_s2fill(int*, double (*)(int*,int*,bool*), double*);
  #define fssp_s2fpdf FC_FUNC(ssp_s2fpdf,SSP_S2FPDF)
    void fssp_s2fpdf(int*, int*, double*, int*, double*);
  #define fssp_s2f123 FC_FUNC(ssp_s2f123,SSP_S2F123)
    void fssp_s2f123(int*, int*, double*, int*, double*);
  #define fdsp_funs2 FC_FUNC(dsp_funs2,DSP_FUNS2)
    double fdsp_funs2(int*,double*,double*,int*);
  #define fdsp_ints2 FC_FUNC(dsp_ints2,DSP_INTS2)
    double fdsp_ints2(int*,double*,double*,double*,double*,double*,int*);
  #define fssp_extrapu FC_FUNC(ssp_extrapu,SSP_EXTRAPU)
    void fssp_extrapu(int*, int*);
  #define fssp_extrapv FC_FUNC(ssp_extrapv,SSP_EXTRAPV)
    void fssp_extrapv(int*, int*);
  #define fisp_splinetype FC_FUNC(isp_splinetype,ISP_SPLINETYPE)
    int fisp_splinetype(int*);
  #define fssp_splims FC_FUNC(ssp_splims,SSP_SPLIMS)
    void fssp_splims(int*, int*, double*, double*, int*, double*,
                     double*, int*);
  #define fssp_unodes FC_FUNC(ssp_unodes,SSP_UNODES)
    void fssp_unodes(int*, double*, int*, int*);
  #define fssp_vnodes FC_FUNC(ssp_vnodes,SSP_VNODES)
    void fssp_vnodes(int*, double*, int*, int*);
  #define fssp_nprint FC_FUNC(ssp_nprint,SSP_NPRINT)
    void fssp_nprint(int*);
  #define fdsp_rscut FC_FUNC(dsp_rscut,DSP_RSCUT)
    double fdsp_rscut(int*);
  #define fdsp_rsmax FC_FUNC(dsp_rsmax,DSP_RSMAX)
    double fdsp_rsmax(int*,double*);
  #define fssp_mprint FC_FUNC(ssp_mprint,SSP_MPRINT)
    void fssp_mprint();
  #define fisp_spsize FC_FUNC(isp_spsize,ISP_SPSIZE)
    int fisp_spsize(int*);
  #define fssp_erase FC_FUNC(ssp_erase,SSP_ERASE)
    void fssp_erase(int*);
  #define fssp_spdumpcpp FC_FUNC(ssp_spdumpcpp,SSP_SPDUMPCPP)
    void fssp_spdumpcpp(int*,char*,int*);
  #define fisp_spreadcpp FC_FUNC(isp_spreadcpp,ISP_SPREADCPP)
    int fisp_spreadcpp(char*,int*);
  #define fssp_spsetval FC_FUNC(ssp_spsetval,SSP_SPSETVAL)
    void fssp_spsetval(int*,int*,double*);
  #define fdsp_spgetval FC_FUNC(dsp_spgetval,DSP_SPGETVAL)
    double fdsp_spgetval(int*,int*);

  /**************************************************************/
  /*  QCDNUM address function from usrGlobalId.f                */
  /**************************************************************/

  #define fipdftab FC_FUNC(ipdftab,IPDFTAB)
    int fipdftab(int*,int*);

  /**************************************************************/
  /*  QCDNUM convolution engine from usrcvol.f                  */
  /**************************************************************/

  #define fmakewta FC_FUNC(makewta,MAKEWTA)
    void fmakewta(double*,int*,double(*)(double*,double*,int*),
                  double(*)(double*));
  #define fmakewtb FC_FUNC(makewtb,MAKEWTB)
    void fmakewtb(double*,int*,double(*)(double*,double*,int*),
                  double(*)(double*),int*);
  #define fmakewrs FC_FUNC(makewrs,MAKEWRS)
    void fmakewrs(double*,int*,double(*)(double*,double*,int*),
                  double(*)(double*,double*,int*),
                  double(*)(double*),int*);
  #define fmakewtd FC_FUNC(makewtd,MAKEWTD)
    void fmakewtd(double*,int*,double(*)(double*,double*,int*),
                  double(*)(double*));
  #define fmakewtx FC_FUNC(makewtx,MAKEWTX)
    void fmakewtx(double*,int*);
  #define fidspfuncpp FC_FUNC(idspfuncpp,IDSPFUNCPP)
    int fidspfuncpp(char*,int*,int*,int*);
  #define fscalewt FC_FUNC(scalewt,SCALEWT)
    void fscalewt(double*,double*,int*);
  #define fcopywgt FC_FUNC(copywgt,COPYWGT)
    void fcopywgt(double*,int*,int*,int*);
  #define fwcrossw FC_FUNC(wcrossw,WCROSSW)
    void fwcrossw(double*,int*,int*,int*,int*);
  #define fwtimesf FC_FUNC(wtimesf,WTIMESF)
    void fwtimesf(double*,double(*)(int*,int*),int*,int*,int*);
  #define ffcrossk FC_FUNC(fcrossk,FCROSSK)
    double ffcrossk(double*,int*,int*,int*,int*,int*);
  #define ffcrossf FC_FUNC(fcrossf,FCROSSF)
    double ffcrossf(double*,int*,int*,int*,int*,int*,int*);
  #define fefromqq FC_FUNC(efromqq,EFROMQQ)
    void fefromqq(double*,double*,int*);
  #define fqqfrome FC_FUNC(qqfrome,QQFROME)
    void fqqfrome(double*,double*,int*);
  #define fstfunxq FC_FUNC(stfunxq,STFUNXQ)
    void fstfunxq(double(*)(int*,int*),double*,double*,double*,
                  int*,int*);

  /**************************************************************/
  /*  QCDNUM error message routines from usrerr.f               */
  /**************************************************************/

  #define fsetumsgcpp FC_FUNC(setumsgcpp,SETUMSGCPP)
    void fsetumsgcpp(char*,int*);
  #define fclrumsg FC_FUNC(clrumsg,CLRUMSG)
    void fclrumsg();

  /**************************************************************/
  /*  QCDNUM toolbox evolution from usrevnn.f                   */
  /**************************************************************/

  #define fevfilla FC_FUNC(evfilla,EVFILLA)
    void fevfilla(double*,int*,double(*)(int*,int*,int*));
  #define fevgetaa FC_FUNC(evgetaa,EVGETAA)
    double fevgetaa(double*,int*,int*,int*,int*);
  #define fevdglap FC_FUNC(evdglap,EVDGLAP)
    void fevdglap(double*,int*,int*,int*,double*,int*,int*,int*,int*,double*);
  #define fcpyparw FC_FUNC(cpyparw,CPYPARW)
    void fcpyparw(double*,double*,int*,int*);
  #define fuseparw FC_FUNC(useparw,USEPARW)
    void fuseparw(double*,int*);
  #define fkeyparw FC_FUNC(keyparw,KEYPARW)
    int fkeyparw(double*,int*);
  #define fkeygrpw FC_FUNC(keygrpw,KEYGRPW)
    int fkeygrpw(double*,int*,int*);
  #define fevpdfij FC_FUNC(evpdfij,EVPDFIJ)
    double fevpdfij(double*,int*,int*,int*,int*);
  #define fevplist FC_FUNC(evplist,EVPLIST)
    void fevplist(double*,int*,double*,double*,double*,int*,int*);
  #define fevtable FC_FUNC(evtable,EVTABLE)
    void fevtable(double*,int*,double*,int*,double*,int*,double*,int*);
  #define fevpcopy FC_FUNC(evpcopy,EVPCOPY)
    void fevpcopy(double*,int*,double*,int*,int*);

  /**************************************************************/
  /*  QCDNUM evolution routines from usrevol.f                  */
  /**************************************************************/

  #define frfromf FC_FUNC(rfromf,RFROMF)
    double frfromf(double*);
  #define fffromr FC_FUNC(ffromr,FFROMR)
    double fffromr(double*);
  #define fasfunc FC_FUNC(asfunc,ASFUNC)
    double fasfunc(double*, int*, int*);
  #define faltabn FC_FUNC(altabn,ALTABN)
    double faltabn(int*, int*, int*, int*);
  #define fevolfg FC_FUNC(evolfg,EVOLFG)
    void fevolfg(int*, double (*)(int* ,double*), double*, int*, double*);
  #define fevsgns FC_FUNC(evsgns,EVSGNS)
    void fevsgns(int*, double (*)(int* ,double*), int*, int*, int*, double*);
  #define fpdfcpy FC_FUNC(pdfcpy,PDFCPY)
    void fpdfcpy(int*, int*);
  #define fextpdf FC_FUNC(extpdf,EXTPDF)
    void fextpdf(double (*)(int*, double* ,double* ,bool*), int*, int*, double*, double*);
  #define fusrpdf FC_FUNC(usrpdf,USRPDF)
    void fusrpdf(double (*)(int*, double* ,double* ,bool*), int*, int*, double*, double*);
  #define fnptabs FC_FUNC(nptabs,NPTABS)
    int fnptabs(int*);
  #define fievtyp FC_FUNC(ievtyp,IEVTYP)
    int fievtyp(int*);

  /**************************************************************/
  /*  QCDNUM fast convolution routines from usrfast.f           */
  /**************************************************************/

  #define fidscope FC_FUNC(idscope,IDSCOPE)
    void fidscope(double*,int*);
  #define ffastini FC_FUNC(fastini,FASTINI)
    void ffastini(double*,double*,int*,int*);
  #define ffastclr FC_FUNC(fastclr,FASTCLR)
    void ffastclr(int*);
  #define ffastinp FC_FUNC(fastinp,FASTINP)
    void ffastinp(double*,int*,double*,int*,int*);
  #define ffastepm FC_FUNC(fastepm,FASTEPM)
    void ffastepm(int*,int*,int*);
  #define ffastsns FC_FUNC(fastsns,FASTSNS)
    void ffastsns(int*,double*,int*,int*);
  #define ffastsum FC_FUNC(fastsum,FASTSUM)
    void ffastsum(int*,double*,int*);
  #define ffastfxk FC_FUNC(fastfxk,FASTFXK)
    void ffastfxk(double*,int*,int*,int*);
  #define ffastfxf FC_FUNC(fastfxf,FASTFXF)
    void ffastfxf(double*,int*,int*,int*,int*);
  #define ffastkin FC_FUNC(fastkin,FASTKIN)
    void ffastkin(int*,double(*)(int*,int*,int*,int*));
  #define ffastcpy FC_FUNC(fastcpy,FASTCPY)
    void ffastcpy(int*,int*,int*);
  #define ffastfxq FC_FUNC(fastfxq,FASTFXQ)
    void ffastfxq(int*,double*,int*);

  /**************************************************************/
  /*  QCDNUM grid routines from usrgrd.f                        */
  /**************************************************************/

  #define fgxmake FC_FUNC(gxmake,GXMAKE)
    void fgxmake(double*, int*, int*, int*, int*, int*);
  #define fixfrmx FC_FUNC(ixfrmx,IXFRMX)
    int fixfrmx(double*);
  #define fxxatix FC_FUNC(xxatix,XXATIX)
    int fxxatix(double*, int*);
  #define fxfrmix FC_FUNC(xfrmix,XFRMIX)
    double fxfrmix(int*);
  #define fgxcopy FC_FUNC(gxcopy,GXCOPY)
    void fgxcopy(double*, int*, int*);
  #define fgqmake FC_FUNC(gqmake,GQMAKE)
    void fgqmake(double*, double*, int*, int*, int*);
  #define fiqfrmq FC_FUNC(iqfrmq,IQFRMQ)
    int fiqfrmq(double*);
  #define fqqatiq FC_FUNC(qqatiq,QQATIQ)
    int fqqatiq(double*, int*);
  #define fqfrmiq FC_FUNC(qfrmiq,QFRMIQ)
    double fqfrmiq(int*);
  #define fgqcopy FC_FUNC(gqcopy,GQCOPY)
    void fgqcopy(double*, int*, int*);
  #define fgrpars FC_FUNC(grpars,GRPARS)
    void fgrpars(int*, double*, double*, int*, double*, double*, int*);
  #define fsetlim FC_FUNC(setlim,SETLIM)
    void fsetlim(int*, int*, int*, double*);
  #define fgetlim FC_FUNC(getlim,GETLIM)
    void fgetlim(int*, double*, double*, double*, double*);

  /**************************************************************/
  /*  QCDNUM init routines from usrini.f                        */
  /**************************************************************/

  #define fqcinitcpp FC_FUNC(qcinitcpp,QCINITCPP)
    void fqcinitcpp(int*, char*, int*);
  #define fsetluncpp FC_FUNC(setluncpp,SETLUNCPP)
    void fsetluncpp(int*, char*, int*);
  #define fnxtlun FC_FUNC(nxtlun,NXTLUN)
    int fnxtlun(int*);
  #define fsetvalcpp FC_FUNC(setvalcpp,SETVALCPP)
    void fsetvalcpp(char*, int*, double*);
  #define fgetvalcpp FC_FUNC(getvalcpp,GETVALCPP)
    void fgetvalcpp(char*, int*, double*);
  #define fsetintcpp FC_FUNC(setintcpp,SETINTCPP)
    void fsetintcpp(char*, int*, int*);
  #define fgetintcpp FC_FUNC(getintcpp,GETINTCPP)
    void fgetintcpp(char*, int*, int*);
  #define fqstorecpp FC_FUNC(qstorecpp,QSTORECPP)
    void fqstorecpp(char*, int*, int*, double*);

  /**************************************************************/
  /*  QCDNUM parameter routines from usrparams.f                */
  /**************************************************************/

  #define fsetord FC_FUNC(setord,SETORD)
    void fsetord(int*);
  #define fgetord FC_FUNC(getord,GETORD)
    void fgetord(int*);
  #define fsetalf FC_FUNC(setalf,SETALF)
    void fsetalf(double*, double*);
  #define fgetalf FC_FUNC(getalf,GETALF)
    void fgetalf(double*, double*);
  #define fsetcbt FC_FUNC(setcbt,SETCBT)
    void fsetcbt(int*, int*, int*, int*);
  #define fmixfns FC_FUNC(mixfns,MIXFNS)
    void fmixfns(int*, double*, double*, double*);
  #define fgetcbt FC_FUNC(getcbt,GETCBT)
    void fgetcbt(int*, double*, double*, double*);
  #define fsetabr FC_FUNC(setabr,SETABR)
    void fsetabr(double*, double*);
  #define fgetabr FC_FUNC(getabr,GETABR)
    void fgetabr(double*, double*);
  #define fnfrmiq FC_FUNC(nfrmiq,NFRMIQ)
    int fnfrmiq(int*, int*, int*);
  #define fnflavs FC_FUNC(nflavs,NFLAVS)
    int fnflavs(int*, int*);
  #define fcpypar FC_FUNC(cpypar,CPYPAR)
    void fcpypar(double*, int*, int*);
  #define fusepar FC_FUNC(usepar,USEPAR)
    void fusepar(int*);
  #define fkeypar FC_FUNC(keypar,KEYPAR)
    int fkeypar(int*);
  #define fkeygrp FC_FUNC(keygrp,KEYGRP)
    int fkeygrp(int*, int*);
  #define fpushcp FC_FUNC(pushcp,PUSHCP)
    void fpushcp();
  #define fpullcp FC_FUNC(pullcp,PULLCP)
    void fpullcp();

  /**************************************************************/
  /*  QCDNUM pdf output routines from usrpdf.f                  */
  /**************************************************************/

  #define fallfxq FC_FUNC(allfxq,ALLFXQ)
    void fallfxq(int*, double*, double*, double*, int*, int*);
  #define fallfij FC_FUNC(allfij,ALLFIJ)
    void fallfij(int*, int*, int*, double*, int*, int*);
  #define fbvalxq FC_FUNC(bvalxq,BVALXQ)
    double fbvalxq(int*, int*, double*, double*, int*);
  #define fbvalij FC_FUNC(bvalij,BVALIJ)
    double fbvalij(int*, int*, int*, int*, int*);
  #define ffvalxq FC_FUNC(fvalxq,FVALXQ)
    double ffvalxq(int*, int*, double*, double*, int*);
  #define ffvalij FC_FUNC(fvalij,FVALIJ)
    double ffvalij(int*, int*, int*, int*, int*);
  #define fsumfxq FC_FUNC(sumfxq,SUMFXQ)
    double fsumfxq(int*, double*, int*, double*, double*, int*);
  #define fsumfij FC_FUNC(sumfij,SUMFIJ)
    double fsumfij(int*, double*, int*, int*, int*, int*);
  #define ffsplne FC_FUNC(fsplne,FSPLNE)
    double ffsplne(int*, int*, double*, int*);
  #define fsplchk FC_FUNC(splchk,SPLCHK)
    double fsplchk(int*, int*, int*);
  #define ffflist FC_FUNC(fflist,FFLIST)
   void ffflist(int*, double*, int*, double*, double*, double*, int*, int*);
  #define ffftabl FC_FUNC(fftabl,FFTABL)
   void ffftabl(int*, double*, int*, double*, int*, double*, int*, double*, int*, int*);
  #define fftable FC_FUNC(ftable,FTABLE)
   void fftable(int*, double*, int*, double*, int*, double*, int*, double*, int*);
  #define fffplotcpp FC_FUNC(ffplotcpp,FFPLOTCPP)
    void fffplotcpp(char*, int*, double (*)(int*, double*, bool*), int*, double*,  double*, int*, char*, int*);
  #define ffiplotcpp FC_FUNC(fiplotcpp,FIPLOTCPP)
    void ffiplotcpp(char*, int*, double (*)(int*, double*, bool*), int*, double*, int*, char*, int*);

  /**************************************************************/
  /*  QCDNUM workspace routines from usrstore.f                 */
  /**************************************************************/

  #define fmaketab FC_FUNC(maketab,MAKETAB)
    void fmaketab(double*,int*,int*,int*,int*,int*,int*);
  #define fsetparw FC_FUNC(setparw,SETPARW)
    void fsetparw(double*,int*,double*,int*);
  #define fgetparw FC_FUNC(getparw,GETPARW)
    void fgetparw(double*,int*,double*,int*);
  #define fdumptabcpp FC_FUNC(dumptabcpp,DUMPTABCPP)
    void fdumptabcpp(double*,int*,int*,char*,int*,char*,int*);
  #define freadtabcpp FC_FUNC(readtabcpp,READTABCPP)
    void freadtabcpp(double*,int*,int*,char*,int*,char*,int*,
                     int*,int*,int*,int*);

  /**************************************************************/
  /*  QCDNUM weight routines from usrwgt.f                      */
  /**************************************************************/

  #define ffillwt FC_FUNC(fillwt,FILLWT)
    void ffillwt(int*, int*, int*, int*);
  #define fdmpwgtcpp FC_FUNC(dmpwgtcpp,DMPWGTCPP)
    void fdmpwgtcpp(int*, int*, char*, int*);
  #define freadwtcpp FC_FUNC(readwtcpp,READWTCPP)
    void freadwtcpp(int*, char*, int*, int*, int*, int*, int*);
  #define fwtfilecpp FC_FUNC(wtfilecpp,WTFILECPP)
    void fwtfilecpp(int*, char*, int*);
  #define fnwused FC_FUNC(nwused,NWUSED)
    void fnwused(int*, int*, int*);

  /**************************************************************/
  /*  HQSTF routines                                            */
  /**************************************************************/

  #define fihqvers FC_FUNC(ihqvers,IHQVERS)
    int fihqvers();
  #define fhqwords FC_FUNC(hqwords,HQWORDS)
    void fhqwords(int*, int*);
  #define fhqparms FC_FUNC(hqparms,HQPARMS)
    void fhqparms(double*, double*, double*);
  #define fhqqfrmu FC_FUNC(hqqfrmu,HQQFRMU)
    double fhqqfrmu(double*);
  #define fhqmufrq FC_FUNC(hqmufrq,HQMUFRQ)
    double fhqmufrq(double*);
  #define fhswitch FC_FUNC(hswitch,HSWITCH)
    void fhswitch(int*);
  #define fhqstfun FC_FUNC(hqstfun,HQSTFUN)
    void fhqstfun(int*, int*, double*, double*, double*, double*, int*, int*);
  #define fhqfillw FC_FUNC(hqfillw,HQFILLW)
    void fhqfillw(int*, double*, double*, double*, int*);
  #define fhqdumpwcpp FC_FUNC(hqdumpwcpp,HQDUMPWCPP)
    void fhqdumpwcpp(int*, char*, int*);
  #define fhqreadwcpp FC_FUNC(hqreadwcpp,HQREADWCPP)
    void fhqreadwcpp(int*, char*, int*, int*, int*);

  /**************************************************************/
  /*  ZMSTF routines                                            */
  /**************************************************************/

  #define fizmvers FC_FUNC(izmvers,IZMVERS)
    int fizmvers();
  #define fzmwords FC_FUNC(zmwords,ZMWORDS)
    void fzmwords(int*, int*);
  #define fzmdefq2 FC_FUNC(zmdefq2,ZMDEFQ2)
    void fzmdefq2(double*, double*);
  #define fzmabval FC_FUNC(zmabval,ZMABVAL)
    void fzmabval(double*, double*);
  #define fzmqfrmu FC_FUNC(zmqfrmu,ZMQFRMU)
    double fzmqfrmu(double*);
  #define fzmufrmq FC_FUNC(zmufrmq,ZMUFRMQ)
    double fzmufrmq(double*);
  #define fzswitch FC_FUNC(zswitch,ZSWITCH)
    void fzswitch(int*);
  #define fzmstfun FC_FUNC(zmstfun,ZMSTFUN)
    void fzmstfun(int*, double*, double*, double*, double*, int*, int*);
  #define fzmfillw FC_FUNC(zmfillw,ZMFILLW)
    void fzmfillw(int*);
  #define fzmdumpwcpp FC_FUNC(zmdumpwcpp,ZMDUMPWCPP)
    void fzmdumpwcpp(int*, char*, int*);
  #define fzmreadwcpp FC_FUNC(zmreadwcpp,ZMREADWCPP)
    void fzmreadwcpp(int*, char*, int*, int*, int*);
  #define fzmwfilecpp FC_FUNC(zmwfilecpp,ZMWFILECPP)
    void fzmwfilecpp(char*,int*);
}

#endif

