
C--   file usrcvol.f containing convolution engine user interface

C--   subroutine MakeWtA(w,id,afun,achi)
C--   subroutine MakeWtB(w,id,bfun,achi,ndel)
C--   subroutine MakeWRS(w,id,rfun,sfun,achi,ndel)
C--   subroutine MakeWtD(w,id,dfun,achi)
C--   subroutine MakeWtX(w,id)

C--   integer function idSpfunCPP(pname,ls,iord,jset)
C--   integer function idSpfun(pname,iord,jset)
C--   integer function iPdftab(jset,id)

C--   subroutine ScaleWt(w,c,id)
C--   subroutine CopyWgt(w,jd1,id2,iadd)
C--   subroutine WcrossW(w,jda,jdb,idc,iadd)
C--   subroutine WtimesF(w,fun,jd1,id2,iadd)

C--   double precision function FcrossK(w,idw,jset,idf,ix,iq)
C--   double precision function FcrossF(w,idw,jset,ida,idb,ix,iq)

C--   subroutine EfromQQ(qvec,evec,nf)
C--   subroutine QQfromE(evec,qvec,nf)

C--   subroutine StfunXq(fun,x,q,f,n,jchk)
C--   subroutine sqcStfLstMpt(subnam,fun,x,q,f,n,jchk)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM convolution engine from usrcvol.f                */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM convolution engine from usrcvol.f                  */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM convolution engine from usrcvol.f                  */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    void makewta(double *w, int jd,
CXXHDR                 double (*afun)(double*,double*,int*),
CXXHDR                 double (*achi)(double*));
C-----------------------------------------------------------------------
CXXHFW  #define fmakewta FC_FUNC(makewta,MAKEWTA)
CXXHFW    void fmakewta(double*,int*,double(*)(double*,double*,int*),
CXXHFW                  double(*)(double*));
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void makewta(double *w, int jd,
CXXWRP               double (*afun)(double*,double*,int*),
CXXWRP               double (*achi)(double*))
CXXWRP  {
CXXWRP    fmakewta(w, &jd, afun, achi);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================
      subroutine MakeWtA(w,jd,afun,achi)
C     ==================================

C--   Make weight table for regular contribution A
C--
C--   w     (in) :   workspace
C--   jd    (in) :   table identifier in global format
C--   afun  (in) :   function declared external in the calling routine 
C--   achi  (in) :   function declared external in the calling routine 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      external afun, achi

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     1   / , iflg /   0,    0  /
      data icma  /    1,     4   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam /'MAKEWTA ( W, ID, AFUN, ACHI )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check id
      igl = iqcSjekId(subnam,'ID',w,jd,icmi,icma,iflg,lint)

C--   Do the work
      call sqcUweitA(w,igl,afun,achi,ierr)
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,'Function achi(qmu2) < 1 encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void makewtb(double *w, int jd,
CXXHDR                 double (*bfun)(double*,double*,int*),
CXXHDR                 double (*achi)(double*),int ndel);
C-----------------------------------------------------------------------
CXXHFW  #define fmakewtb FC_FUNC(makewtb,MAKEWTB)
CXXHFW    void fmakewtb(double*,int*,double(*)(double*,double*,int*),
CXXHFW                  double(*)(double*),int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void makewtb(double *w, int jd,
CXXWRP               double (*bfun)(double*,double*,int*),
CXXWRP               double (*achi)(double*),int ndel)
CXXWRP  {
CXXWRP    fmakewtb(w, &jd, bfun, achi, &ndel);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================================
      subroutine MakeWtB(w,jd,bfun,achi,ndel)
C     =======================================

C--   Make weight table for singular contribution B
C--
C--   w      (in) :  workspace
C--   jd     (in) :  table identifier in global format
C--   bfun   (in) :  function declared external in the calling routine
C--   achi   (in) :  function declared external in the calling routine
C--   ndel   (in) :  1 = no delta(1-x) contribution; otherwise yes

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      external bfun, achi

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign   iset  itype              pzi  cfil
      data icmi  /    1,     1   / , iflg /   0,    0  /
      data icma  /    1,     4   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam  /'MAKEWTB ( W, ID, BFUN, ACHI, NODELTA )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check id
      igl = iqcSjekId(subnam,'ID',w,jd,icmi,icma,iflg,lint)

C--   Do the work (note flip ndel 0 <--> 1)
      call sqcUweitB(w,igl,bfun,achi,1-ndel,ierr)
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,'Function achi(qmu2) < 1 encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void makewrs(double *w, int jd,
CXXHDR                 double (*rfun)(double*,double*,int*),
CXXHDR                 double (*sfun)(double*,double*,int*),
CXXHDR                 double (*achi)(double*),int ndel);
C-----------------------------------------------------------------------
CXXHFW  #define fmakewrs FC_FUNC(makewrs,MAKEWRS)
CXXHFW    void fmakewrs(double*,int*,double(*)(double*,double*,int*),
CXXHFW                  double(*)(double*,double*,int*),
CXXHFW                  double(*)(double*),int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void makewrs(double *w, int jd,
CXXWRP               double (*rfun)(double*,double*,int*),
CXXWRP               double (*sfun)(double*,double*,int*),
CXXWRP               double (*achi)(double*),int ndel)
CXXWRP  {
CXXWRP    fmakewrs(w, &jd, rfun, sfun, achi, &ndel);
CXXWRP  }
C-----------------------------------------------------------------------

C     ============================================
      subroutine MakeWRS(w,jd,rfun,sfun,achi,ndel)
C     ============================================

C--   Make weight table for contribution RS
C--
C--   w      (in) :  workspace
C--   jd     (in) :  table identifier in global format
C--   rfun   (in) :  function declared external in the calling routine
C--   sfun   (in) :  function declared external in the calling routine 
C--   achi   (in) :  function declared external in the calling routine
C--   ndel   (in) :  1 = no delta(1-x) contribution; otherwise yes

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      external rfun, sfun, achi

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     1   / , iflg /   0,    0  /
      data icma  /    1,     4   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam /'MAKEWRS ( W, ID, RFUN, SFUN, ACHI, NODELTA )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check id
      igl = iqcSjekId(subnam,'ID',w,jd,icmi,icma,iflg,lint)

C--   Do the work (note flip ndel 0 <--> 1)
      call sqcUwgtRS(w,igl,rfun,sfun,achi,1-ndel,ierr)
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,'Function achi(qmu2) < 1 encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void makewtd(double *w, int jd,
CXXHDR                 double (*dfun)(double*,double*,int*),
CXXHDR                 double (*achi)(double*));
C-----------------------------------------------------------------------
CXXHFW  #define fmakewtd FC_FUNC(makewtd,MAKEWTD)
CXXHFW    void fmakewtd(double*,int*,double(*)(double*,double*,int*),
CXXHFW                  double(*)(double*));
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void makewtd(double *w, int jd,
CXXWRP               double (*dfun)(double*,double*,int*),
CXXWRP               double (*achi)(double*))
CXXWRP  {
CXXWRP    fmakewtd(w, &jd, dfun, achi);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================
      subroutine MakeWtD(w,jd,dfun,achi)
C     ==================================

C--   Make weight table for factor*delta(1-x) contribution
C--
C--   w      (in) :  workspace
C--   jd     (in) :  table identifier in global format
C--   dfun   (in) :  function declared external in the calling routine
C--   achi   (in) :  function declared external in the calling routine


      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      external dfun,achi

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     1   / , iflg /   0,    0  /
      data icma  /    1,     4   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam  /'MAKEWTD ( W, ID, DFUN, ACHI )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check id
      igl = iqcSjekId(subnam,'ID',w,jd,icmi,icma,iflg,lint)

C--   Do the work
      call sqcUweitD(w,igl,dfun,achi,ierr)
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,'Function achi(qmu2) < 1 encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void makewtx(double *w, int jd);
C-----------------------------------------------------------------------
CXXHFW  #define fmakewtx FC_FUNC(makewtx,MAKEWTX)
CXXHFW    void fmakewtx(double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void makewtx(double *w, int jd)
CXXWRP  {
CXXWRP    fmakewtx(w, &jd);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ========================
      subroutine MakeWtX(w,jd)
C     ========================

C--   Make weight table for FxF convolution
C--
C--   w    (in) :    workspace
C--   jd   (in) :    table identifier in global format

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     1   / , iflg /   0,    0  /
      data icma  /    1,     4   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam /'MAKEWTX ( W, ID )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check id
      igl = iqcSjekId(subnam,'ID',w,jd,icmi,icma,iflg,lint)

C--   Do the work
      call sqcUweitX(w,igl,ierr)
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,'Error condition encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    int idspfun(string pname, int iord, int jset);
C-----------------------------------------------------------------------
CXXHFW  #define fidspfuncpp FC_FUNC(idspfuncpp,IDSPFUNCPP)
CXXHFW    int fidspfuncpp(char*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int idspfun(string pname, int iord, int jset)
CXXWRP  {
CXXWRP    int ls = pname.size();
CXXWRP    char *cpname = new char[ls+1];
CXXWRP    strcpy(cpname,pname.c_str());
CXXWRP    return fidspfuncpp(cpname, &ls, &iord, &jset);
CXXWRP    delete[] cpname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ===============================================
      integer function idSpfunCPP(pname,ls,iord,jset)
C     ===============================================

      implicit double precision (a-h,o-z)

      character*(100) pname

      if(ls.gt.100) stop 'idSpfunCPP: input PIJ size > 100 characters'

      idSpfunCPP = idSpfun(pname(1:ls),iord,jset)

      return
      end


C     =========================================
      integer function idSpfun(pname,iord,jset)
C     =========================================

C--   Returns splitting function table index encoded as -(1000*jset+idspl)
C--   Returns -99999 if splitting function table not available
C--   NB: jset 1=unpol, 2=pol, 3=timelike

      implicit double precision (a-h,o-z)      
      
      include 'qcdnum.inc'
      include 'qstor7.inc'
      
      character*(*) pname
      character*3   input
      character*3   ptab(12)
C--                 1     2     3     4     5     6     7
      data ptab / 'PQQ','PQG','PGQ','PGG','PPL','PMI','PVA',
C--                 8     9    10    11    12
     +            'AGQ','AGG','AQQ','AHQ','AHG' /
      
      idspfun = -1
      if(jset.lt.1. .or. jset.gt.3) return
      if(mxord7(jset).eq.0)         return

      input        = '   '
      len          = min(imb_lenoc(pname),3)
      input(1:len) = pname(1:len)
      call smb_cltou(input)
      id = 0
      do i = 1,7
        if(input.eq.ptab(i))   id = idPij7(i,iord,jset)
      enddo
      if(input.eq.ptab( 8)) id = idAijk7(1,2,iord,jset)
      if(input.eq.ptab( 9)) id = idAijk7(1,1,iord,jset)
      if(input.eq.ptab(10)) id = idAijk7(2,2,iord,jset)
      if(input.eq.ptab(11)) id = idAijk7(3,2,iord,jset)
      if(input.eq.ptab(12)) id = idAijk7(3,1,iord,jset)

      if(id.eq.0) then
        idspfun = -99999
      else
        idspfun = -id
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void scalewt(double *w, double c, int jd);
C-----------------------------------------------------------------------
CXXHFW  #define fscalewt FC_FUNC(scalewt,SCALEWT)
CXXHFW    void fscalewt(double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void scalewt(double *w, double c, int jd)
CXXWRP  {
CXXWRP    fscalewt(w, &c, &jd);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================
      subroutine ScaleWt(w,c,jd)
C     ==========================

C--   Multiply weight table by a constant
C--
C--   w   (in)     workspace
C--   c   (in)     constant
C--   jd  (in)     weight table identifier in w (global format)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     1   / , iflg /   0,    0  /
      data icma  /    1,     4   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam /'SCALEWT ( W, C, ID )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check id
      igl = iqcSjekId(subnam,'ID',w,jd,icmi,icma,iflg,lint)

C--   Do the work
      call sqcScaleWt(w,c,igl)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void copywgt(double *w, int jd1, int jd2, int iadd);
C-----------------------------------------------------------------------
CXXHFW  #define fcopywgt FC_FUNC(copywgt,COPYWGT)
CXXHFW    void fcopywgt(double*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void copywgt(double *w, int jd1, int jd2, int iadd)
CXXWRP  {
CXXWRP    fcopywgt(w, &jd1, &jd2, &iadd);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================
      subroutine CopyWgt(w,jd1,jd2,iadd)
C     ==================================

C--   Add content of jd1 to jd2
C--
C--   w    (in) : workspace
C--   jd1  (in) : input id (global format), < 0 table in internal memory
C--   jd2  (in) : output id (global format) in workspace w
C--   iadd (in) : -1,0,1 = subtract/copy/add id1 to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi1(2), icma1(2), iflg1(2)
C--                isign  itype               pzi  cfil
      data icmi1  /   -1,     1   / , iflg1 /   0,    1  /
      data icma1  /    1,     6   /

      dimension icmi2(2), icma2(2), iflg2(2)
C--                isign  itype               pzi  cfil
      data icmi2  /    1,     1   / , iflg2 /   0,    0  /
      data icma2  /    1,     6   /

      logical lint1, lint2

      dimension w(*)

C--   iotyp(ityp_in,ityp_out) = 0 --> combination not allowed      
      dimension iotyp(6,6)
C--   ityp_in      1  2  3  4  5  6
      data iotyp / 1, 0, 0, 0, 0, 0,   !ityp_out 1
     +             1, 1, 0, 0, 0, 0,   !ityp_out 2
     +             1, 0, 1, 0, 0, 0,   !ityp_out 3
     +             1, 1, 1, 1, 0, 0,   !ityp_out 4
     +             0, 0, 0, 0, 1, 0,   !ityp_out 5
     +             0, 0, 0, 0, 0, 1 /  !ityp_out 6

      character*80 subnam
      data subnam /'COPYWGT ( W, ID1, ID2, IADD )'/

C--   Avoid compiler warning
      ifirst = 0
C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Overwrite, no thank you...
      if(jd2.eq.jd1) then
        call sqcErrMsg(subnam,'ID2 cannot be equal to ID1')
      endif      

C--   Check iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')

C--   Check id1
      igl1 = iqcSjekId(subnam,'ID1',w,jd1,icmi1,icma1,iflg1,lint1)

C--   Check id2
      igl2 = iqcSjekId(subnam,'ID2',w,jd2,icmi2,icma2,iflg2,lint2)

C--   Check input/output table types
      call sqcChkTyp12(subnam,'ID1','ID2',igl1,igl2,iotyp)

C--   Do the work
      if(lint1) then
        call sqcChkIoy12(subnam,'ID1','ID2',stor7,igl1,w,igl2)
        call sqcCopyWt(stor7,igl1,w,igl2,iadd)
      else
        call sqcChkIoy12(subnam,'ID1','ID2',w,igl1,w,igl2)
        call sqcCopyWt(w,igl1,w,igl2,iadd)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void wcrossw(double *w, int jda, int jdb, int jdc, int iadd);
C-----------------------------------------------------------------------
CXXHFW  #define fwcrossw FC_FUNC(wcrossw,WCROSSW)
CXXHFW    void fwcrossw(double*,int*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void wcrossw(double *w, int jda, int jdb, int jdc, int iadd)
CXXWRP  {
CXXWRP    fwcrossw(w, &jda, &jdb, &jdc, &iadd);
CXXWRP  }
C-----------------------------------------------------------------------

C     ======================================
      subroutine WcrossW(w,jda,jdb,jdc,iadd)
C     ======================================

C--   Make weight table for convolution C = A cross B

C--   w    (in) : workspace
C--   jda  (in) : input id (global format), < 0 table in internal memory
C--   jdb  (in) : input id (global format), < 0 table in internal memory
C--   jdc  (in) : output id (global format) in workspace w
C--   iadd (in) : -1,0,1 = subtract/copy/add id1 to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmia(2), icmaa(2), iflga(2)
C--                isign  itype               pzi  cfil
      data icmia  /   -1,     1   / , iflga /   0,    1  /
      data icmaa  /    1,     4   /
      dimension icmib(2), icmab(2), iflgb(2)
C--                isign  itype               pzi  cfil
      data icmib  /   -1,     1   / , iflgb /   0,    1  /
      data icmab  /    1,     4   /
      dimension icmic(2), icmac(2), iflgc(2)
C--                isign  itype               pzi  cfil
      data icmic  /    1,     1   / , iflgc /   0,    0  /
      data icmac  /    1,     4   /

      logical linta, lintb, lintc

      dimension w(*)
      
C--   iotyp(ityp_in,ityp_out) = 0 --> combination not allowed      
      dimension iotyp(6,6)
C--   ityp_in      1  2  3  4  5  6
      data iotyp / 1, 0, 0, 0, 0, 0,   !ityp_out 1
     +             1, 1, 0, 0, 0, 0,   !ityp_out 2
     +             1, 0, 1, 0, 0, 0,   !ityp_out 3
     +             1, 1, 1, 1, 0, 0,   !ityp_out 4
     +             0, 0, 0, 0, 0, 0,   !ityp_out 5
     +             0, 0, 0, 0, 0, 0 /  !ityp_out 6

      character*80 subnam
      data subnam /'WCROSSW ( W, IDA, IDB, IDC, IADD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Overwrite, no thank you...
      if(jdc.eq.jda .or. jdc.eq.jdb) then
        call sqcErrMsg(subnam,'IDC cannot be equal to IDA or IDB')
      endif
      
C--   Check iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')
      
C--   Check ida
      igla = iqcSjekId(subnam,'IDA',w,jda,icmia,icmaa,iflga,linta)

C--   Check idb
      iglb = iqcSjekId(subnam,'IDB',w,jdb,icmib,icmab,iflgb,lintb)

C--   Check idc
      iglc = iqcSjekId(subnam,'IDC',w,jdc,icmic,icmac,iflgc,lintc)

C--   Check input/output table types
      call sqcChkTyp12(subnam,'IDA','IDC',igla,iglc,iotyp)
      call sqcChkTyp12(subnam,'IDA','IDC',iglb,iglc,iotyp)

C--   Global scratch table identifiers
      idt1 = iqcIdPdfLtoG(0,1)
      idt2 = iqcIdPdfLtoG(0,2)
C--   Do the work
      if(linta .and. lintb)     then
        call sqcChkIoy12(subnam,'IDA','IDC',stor7,igla,w,iglc)
        call sqcChkIoy12(subnam,'IDB','IDC',stor7,iglb,w,iglc)
        call sqcWcrossW(stor7,igla,stor7,iglb,w,iglc,idt1,idt2,iadd)
      elseif(linta .and. .not.lintb) then
        call sqcChkIoy12(subnam,'IDA','IDC',stor7,igla,w,iglc)
        call sqcChkIoy12(subnam,'IDB','IDC',w    ,iglb,w,iglc)
        call sqcWcrossW(stor7,igla,w,iglb,w,iglc,idt1,idt2,iadd)
      elseif(.not.linta .and. lintb) then
        call sqcChkIoy12(subnam,'IDA','IDC',w    ,igla,w,iglc)
        call sqcChkIoy12(subnam,'IDB','IDC',stor7,iglb,w,iglc)
        call sqcWcrossW(w,igla,stor7,iglb,w,iglc,idt1,idt2,iadd)
      else
        call sqcChkIoy12(subnam,'IDA','IDC',w    ,igla,w,iglc)
        call sqcChkIoy12(subnam,'IDB','IDC',w    ,iglb,w,iglc)
        call sqcWcrossW(w,igla,w,iglb,w,iglc,idt1,idt2,iadd)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void wtimesf(double *w, double (*fun)(int*,int*),
CXXHDR                 int jd1, int jd2, int iadd);
C-----------------------------------------------------------------------
CXXHFW  #define fwtimesf FC_FUNC(wtimesf,WTIMESF)
CXXHFW    void fwtimesf(double*,double(*)(int*,int*),int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void wtimesf(double *w, double (*fun)(int*,int*),
CXXWRP               int jd1, int jd2, int iadd)
CXXWRP  {
CXXWRP    fwtimesf(w, fun, &jd1, &jd2, &iadd);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ======================================
      subroutine WtimesF(w,fun,jd1,jd2,iadd)
C     ======================================

C--   Multiply id1 by fun(iq,nf) and store result in id2
C--
C--   w    (in) : workspace
C--   fun  (in) : user defined function
C--   jd1  (in) : input id (global format), < 0 table in internal memory
C--   jd2  (in) : output id (global format) in workspace w
C--   iadd (in) : -1,0,1  subtract,store,add result to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi1(2), icma1(2), iflg1(2)
C--                isign  itype               pzi  cfil
      data icmi1  /   -1,     1   / , iflg1 /   0,    1  /
      data icma1  /    1,     4   /
      dimension icmi2(2), icma2(2), iflg2(2)
C--                isign  itype               pzi  cfil
      data icmi2  /    1,     1   / , iflg2 /   0,    0  /
      data icma2  /    1,     4   /

      logical lint1, lint2

      dimension w(*)
      
C--   iotyp(ityp_in,ityp_out) = 0 --> combination not allowed
      dimension iotyp0(6,6)
C--   ityp_in      1  2  3  4  5  6    !fun depends on nothing
      data iotyp0/ 1, 0, 0, 0, 0, 0,   !ityp_out 1
     +             1, 1, 0, 0, 0, 0,   !ityp_out 2
     +             1, 0, 1, 0, 0, 0,   !ityp_out 3
     +             1, 1, 1, 1, 0, 0,   !ityp_out 4
     +             0, 0, 0, 0, 0, 0,   !ityp_out 5
     +             0, 0, 0, 0, 0, 0 /  !ityp_out 6
      dimension iotyp1(6,6)
C--   ityp_in      1  2  3  4  5  6    !fun depends on iq only
      data iotyp1/ 0, 0, 0, 0, 0, 0,   !ityp_out 1
     +             0, 1, 0, 0, 0, 0,   !ityp_out 2
     +             1, 0, 1, 0, 0, 0,   !ityp_out 3
     +             1, 1, 1, 1, 0, 0,   !ityp_out 4
     +             0, 0, 0, 0, 0, 0,   !ityp_out 5
     +             0, 0, 0, 0, 0, 0 /  !ityp_out 6
      dimension iotyp2(6,6)
C--   ityp_in      1  2  3  4  5  6    !fun depends on nf only
      data iotyp2/ 0, 0, 0, 0, 0, 0,   !ityp_out 1
     +             1, 1, 0, 0, 0, 0,   !ityp_out 2
     +             0, 0, 0, 0, 0, 0,   !ityp_out 3
     +             1, 1, 1, 1, 0, 0,   !ityp_out 4
     +             0, 0, 0, 0, 0, 0,   !ityp_out 5
     +             0, 0, 0, 0, 0, 0 /  !ityp_out 6
      dimension iotyp3(6,6)
C--   ityp_in      1  2  3  4  5  6    !fun depends on both iq and nf
      data iotyp3/ 0, 0, 0, 0, 0, 0,   !ityp_out 1
     +             0, 0, 0, 0, 0, 0,   !ityp_out 2
     +             0, 0, 0, 0, 0, 0,   !ityp_out 3
     +             1, 1, 1, 1, 0, 0,   !ityp_out 4
     +             0, 0, 0, 0, 0, 0,   !ityp_out 5
     +             0, 0, 0, 0, 0, 0 /  !ityp_out 6

      external fun

      character*80 subnam
      data subnam /'WTIMESF ( W, FUN, ID1, ID2, IADD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')

C--   Check id1
      igl1 = iqcSjekId(subnam,'ID1',w,jd1,icmi1,icma1,iflg1,lint1)

C--   Check id2
      igl2 = iqcSjekId(subnam,'ID2',w,jd2,icmi2,icma2,iflg2,lint2)

C--   Check fun dependence
      idep = 0                                    !no dpendence
      if(fun(1,3).ne.fun(2,3)) idep = idep+1      !iq dependent
      if(fun(1,3).ne.fun(1,4)) idep = idep+2      !nf dependent

C--   Check input/output table types
      if(idep.eq.0) then
        call sqcChkTyp12(subnam,'ID1','ID2',igl1,igl2,iotyp0)
      elseif(idep.eq.1) then
        call sqcChkTyp12(subnam,'ID1','ID2',igl1,igl2,iotyp1)
      elseif(idep.eq.2) then
        call sqcChkTyp12(subnam,'ID1','ID2',igl1,igl2,iotyp2)
      elseif(idep.eq.3) then
        call sqcChkTyp12(subnam,'ID1','ID2',igl1,igl2,iotyp3)
      endif

C--   Do the work
      if(lint1) then
        call sqcChkIoy12(subnam,'ID1','ID2',stor7,igl1,w,igl2)
        call sqcWtimesF(fun,stor7,igl1,w,igl2,iadd)
      else
        call sqcChkIoy12(subnam,'ID1','ID2',w    ,igl1,w,igl2)
        call sqcWtimesF(fun,w,igl1,w,igl2,iadd)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C==   ===============================================================
C==   Convolution engine ============================================
C==   ===============================================================

C-----------------------------------------------------------------------
CXXHDR    double fcrossk(double *w, int idw, int jdum, int idf,
CXXHDR                   int ix, int iq);
C-----------------------------------------------------------------------
CXXHFW  #define ffcrossk FC_FUNC(fcrossk,FCROSSK)
CXXHFW    double ffcrossk(double*,int*,int*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double fcrossk(double *w, int idw, int jdum, int idf,
CXXWRP                 int ix, int iq)
CXXWRP  {
CXXWRP    return ffcrossk(w, &idw, &jdum, &idf, &ix, &iq);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================================================
      double precision function FcrossK(w,idw,jdum,idf,ix,iq)
C     =======================================================

C--   Convolution F cross K
C--
C--   w     (in) :   workspace
C--   idw   (in) :   weight table identifier (global format)
C--   jdum  (in) :   dummy variable
C--   idf   (in) :   identifier of pdf F (global format)
C--   ix    (in) :   x-grid index
C--   iq    (in) :   mu2 grid index can be < 0

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmiw(2), icmaw(2), iflgw(2)
C--                isign  itype               pzi  cfil
      data icmiw  /   -1,     1   / , iflgw /   0,    1  /
      data icmaw  /    1,     4   /
      dimension icmif(2), icmaf(2), iflgf(2)
C--                isign  itype               pzi  cfil
      data icmif  /   -1,     5   / , iflgf /   0,    1  /
      data icmaf  /    1,     5   /

      logical lintw, lintf

      dimension w(*)

      character*80 subnam
      data subnam /'FCROSSK ( W, IDW, IDUM, IDF, IX, IQ )'/

      jset = jdum !avoid compoler warning

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check idw
      iglw = iqcSjekId(subnam,'IDW',w,idw,icmiw,icmaw,iflgw,lintw)

C--   Check idf
      iglf = iqcSjekId(subnam,'IDF',w,idf,icmif,icmaf,iflgf,lintf)

C--   Catch x = 1
      if(ix.eq.nyy2(0)+1) then
        FcrossK = 0.D0
        return
      endif

C--   Check cuts on pdf and point to correct parameters
      if(lintf) then
        ifail = iqcChkLmij(subnam,stor7,iglf,ix,abs(iq),1)
        ksetw = iglf/1000
        ipver = int(dparGetPar(stor7,ksetw,idipver8))
      else
        ifail = iqcChkLmij(subnam,  w  ,iglf,ix,abs(iq),1)
        ksetw = iglf/1000
        ipver = int(dparGetPar(w,ksetw,idipver8))
      endif

C--   Point to the correct parameters
      call sparParTo5(ipver)

C--   Do the work
      iy = nyy2(0) + 1 - ix
      if(lintw .and. lintf)     then
        FcrossK = dqcFcrossK(stor7,iglw,stor7,iglf,iy,iq)
      elseif(lintw .and. .not.lintf) then
        FcrossK = dqcFcrossK(stor7,iglw,w    ,iglf,iy,iq)
      elseif(.not.lintw .and. lintf) then
        FcrossK = dqcFcrossK(w    ,iglw,stor7,iglf,iy,iq)
      else
        FcrossK = dqcFcrossK(w    ,iglw,w    ,iglf,iy,iq)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    double fcrossf(double *w, int idw, int jdum, int ida, int idb,
CXXHDR                   int ix, int iq);
C-----------------------------------------------------------------------
CXXHFW  #define ffcrossf FC_FUNC(fcrossf,FCROSSF)
CXXHFW    double ffcrossf(double*,int*,int*,int*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double fcrossf(double *w, int idw, int jdum, int ida, int idb,
CXXWRP                 int ix, int iq)
CXXWRP  {
CXXWRP    return ffcrossf(w, &idw, &jdum, &ida, &idb, &ix, &iq);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ===========================================================
      double precision function FcrossF(w,idw,jdum,ida,idb,ix,iq)
C     ===========================================================

C--   Convolution Fa cross Fb
C--
C--   w     (in) :   store dimensioned in the calling routine
C--   idw   (in) :   weight table identifier (global format)
C--   jdum  (in) :   dummy variable
C--   ida   (in) :   identifier Fa (global format)
C--   idb   (in) :   identifier Fb (global format)
C--   ix    (in) :   x-grid index
C--   iq    (in) :   mu2 grid index can be < 0

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmiw(2), icmaw(2), iflgw(2)
C--                isign  itype               pzi  cfil  iset
      data icmiw  /    1,     1   / , iflgw /   0,    1  /
      data icmaw  /    1,     4   /
      dimension icmia(2), icmaa(2), iflga(2)
C--                isign  itype               pzi  cfil  iset
      data icmia  /   -1,     5   / , iflga /   0,    1  /
      data icmaa  /    1,     5   /
      dimension icmib(2), icmab(2), iflgb(2)
C--                isign  itype               pzi  cfil  iset
      data icmib  /   -1,     5   / , iflgb /   0,    1  /
      data icmab  /    1,     5   /

      logical lintw, linta, lintb

      dimension w(*)

      character*80 subnam
      data subnam /'FCROSSF ( W, IDW, IDUM, IDA, IDB, IX, IQ )'/

      jset = jdum !avoid compoler warning

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check idw
      iglw = iqcSjekId(subnam,'IDW',w,idw,icmiw,icmaw,iflgw,lintw)

C--   Check ida
      igla = iqcSjekId(subnam,'IDA',w,ida,icmia,icmaa,iflga,linta)

C--   Check idb
      iglb = iqcSjekId(subnam,'IDB',w,idb,icmib,icmab,iflgb,lintb)

C--   Catch x = 1
      if(ix.eq.nyy2(0)+1) then
        FcrossF = 0.D0
        return
      endif

C--   Check cuts on pdf a
      if(linta) then
        ifail = iqcChkLmij(subnam,stor7,igla,ix,abs(iq),1)
        kseta = igla/1000
        keya  = int(dparGetPar(stor7,kseta,idipver8))
      else
        ifail = iqcChkLmij(subnam,  w  ,igla,ix,abs(iq),1)
        kseta = igla/1000
        keya  = int(dparGetPar(w,kseta,idipver8))
      endif

C--   Check cuts on pdf b
      if(lintb) then
        ifail = iqcChkLmij(subnam,stor7,iglb,ix,abs(iq),1)
        ksetb = iglb/1000
        keyb  = int(dparGetPar(stor7,ksetb,idipver8))
      else
        ifail = iqcChkLmij(subnam,  w  ,iglb,ix,abs(iq),1)
        keyb  = iglb/1000
        ipver = int(dparGetPar(w,ksetb,idipver8))
      endif

C--   Now see if pdfs can be combined
      if(keya.ne.keyb) call sqcErrMsg(subnam,
     +    'Cannot combine pdfs with different evolution parameters')

C--   Point to correct parameters
      call sparParTo5(keya)

C--   Do the work
      iy = nyy2(0) + 1 - ix
      if    (     lintw .and.      linta .and.      lintb)     then
        FcrossF = dqcFcrossF(stor7,iglw,stor7,igla,stor7,iglb,iy,iq)
      elseif(     lintw .and.      linta .and. .not.lintb)     then
        FcrossF = dqcFcrossF(stor7,iglw,stor7,igla,w    ,iglb,iy,iq)
      elseif(     lintw .and. .not.linta .and.      lintb)     then
        FcrossF = dqcFcrossF(stor7,iglw,w    ,igla,stor7,iglb,iy,iq)
      elseif(     lintw .and. .not.linta .and. .not.lintb)     then
        FcrossF = dqcFcrossF(stor7,iglw,w    ,igla,w    ,iglb,iy,iq)
      elseif(.not.lintw .and.      linta .and.      lintb)     then
        FcrossF = dqcFcrossF(w    ,iglw,stor7,igla,stor7,iglb,iy,iq)
      elseif(.not.lintw .and.      linta .and. .not.lintb)     then
        FcrossF = dqcFcrossF(w    ,iglw,stor7,igla,w    ,iglb,iy,iq)
      elseif(.not.lintw .and. .not.linta .and.      lintb)     then
        FcrossF = dqcFcrossF(w    ,iglw,w    ,igla,stor7,iglb,iy,iq)
      else
        FcrossF = dqcFcrossF(w    ,iglw,w    ,igla,w    ,iglb,iy,iq)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C==   ==================================================================
C==   Utility routines =================================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void efromqq(double *qvec, double *evec, int nf);
C-----------------------------------------------------------------------
CXXHFW  #define fefromqq FC_FUNC(efromqq,EFROMQQ)
CXXHFW    void fefromqq(double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void efromqq(double *qvec, double *evec, int nf)
CXXWRP  {
CXXWRP    fefromqq(qvec, evec, &nf);
CXXWRP  }
C-----------------------------------------------------------------------

C     ================================
      subroutine EfromQQ(qvec,evec,nf)
C     ================================

C--   Transform coefficients from flavor basis to si/ns basis
C--
C--         -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   qvec  tb bb cb sb ub db  g  d  u  s  c  b  t
C--
C--   |qpm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--   |epm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  si 2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpdfs7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension qvec(-6:6), evec(12)

      character*80 subnam
      data subnam /'EFROMQQ ( QVEC, EVEC, NF )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Do the work      
      call sqcEfromQQ(qvec,evec,nf,nf)
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void qqfrome(double *evec, double *qvec, int nf);
C-----------------------------------------------------------------------
CXXHFW  #define fqqfrome FC_FUNC(qqfrome,QQFROME)
CXXHFW    void fqqfrome(double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void qqfrome(double *evec, double *qvec, int nf)
CXXWRP  {
CXXWRP    fqqfrome(evec, qvec, &nf);
CXXWRP  }
C-----------------------------------------------------------------------

C     ================================
      subroutine QQfromE(evec,qvec,nf)
C     ================================

C--   Transform si/ns basis to flavor basis
C--
C--   |epm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  si 2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-
C--
C--   |qpm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--         -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   qvec  tb bb cb sb ub db  g  d  u  s  c  b  t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpdfs7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension evec(12), qpm(12), qvec(-6:6)

      character*80 subnam
      data subnam /'QQFROME ( EVEC, QVEC, NF )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Transform  (eq 2.29)
      do i = 1,nf
        biplu = 0.D0
        bimin = 0.D0
        do j = 1,nf
          biplu = biplu + evec(j  )*umateq7(j,i,nf)
          bimin = bimin + evec(6+j)*umateq7(j,i,nf)
        enddo
        qpm(i  ) = biplu
        qpm(i+6) = bimin
      enddo
C--   Get q,qbar coeff from qpm coeff  (eq 2.20)
      do i = -6,6
        qvec(i) = 0.D0
      enddo
      do i = 1,nf
        qvec( i) = qpm(i) + qpm(6+i)
        qvec(-i) = qpm(i) - qpm(6+i)
      enddo
 
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C==   ===============================================================
C==   Interpolation =================================================
C==   ===============================================================

C-----------------------------------------------------------------------
CXXHDR    void stfunxq(double (*fun)(int*,int*), double *x, double *q,
CXXHDR                 double *f, int n, int jchk);
C-----------------------------------------------------------------------
CXXHFW  #define fstfunxq FC_FUNC(stfunxq,STFUNXQ)
CXXHFW    void fstfunxq(double(*)(int*,int*),double*,double*,double*,
CXXHFW                  int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void stfunxq(double (*fun)(int*,int*), double *x, double *q,
CXXWRP               double *f, int n, int jchk)
CXXWRP  {
CXXWRP    fstfunxq(fun, x, q, f, &n, &jchk);
CXXWRP  }
C-----------------------------------------------------------------------

C     ====================================
      subroutine StfunXq(fun,x,q,f,n,jchk)
C     ====================================

C--   Interpolate fun(ix,iq)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      external fun
      dimension x(*),q(*),f(*)

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'STFUNXQ ( STFUN, X, QMU2, STF, N, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check number of interpolation points
      if(n.le.0) call sqcErrMsg(subnam,'N should be larger than zero')

C--   Start scoping
      Lscopechek6 = .true.
      iscopekey6  = int(dparGetPar(pars8,iscopeslot6,idipver8))

C--   Point to correct parameters
      call sparParTo5(iscopekey6)

C--   Dummy call to fun
*      dummy = fun(-1,0)

      nmax  = mpt0
*      nmax  = 10
      nlast = 0
      ntodo = min(n,nmax)
      do while( ntodo.gt.0 )
        i1    = nlast+1
        call sqcStfLstMpt(subnam,fun,x(i1),q(i1),f(i1),ntodo,jchk)
        nlast = nlast+ntodo
        ntodo = min(n-nlast,nmax)
      enddo

C--   Terminate scoping
      Lscopechek6 = .false.

      return
      end

C     ================================================
      subroutine sqcStfLstMpt(subnam,fun,x,q,f,n,jchk)
C     ================================================

C--   Fast interpolation of max mpt0 points
C--
C--   subnam     (in) : subroutine name for error message
C--   fun(ix,iq) (in) : function to interpolate
C--   x,q        (in) : list of interpolation points
C--   f         (out) : list of interpolated results
C--   n          (in) : number of items in x,q,f <= mpt0
C--   jchk       (in) : error checking flag

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical lmb_eq

      external fun

      dimension x(*), q(*), f(*)
      dimension yy(mpt0),tt(mpt0),ipoint(mpt0),ff(mpt0)

      logical lqcInside

      character*80 subnam

C--   Workspace
      dimension ww(11+44*mpt0)

C--   Weed points outside grid and convert to y,t
      npt = 0
      do i = 1,n
        if(lmb_eq(x(i),1.D0,-aepsi6)) then
          f(i)        =  0.D0
        elseif(lqcInside(x(i),q(i))) then
          f(i)        =  0.D0
          npt         =  npt+1
          yy(npt)     = -log(x(i))
          tt(npt)     =  log(q(i))
          ipoint(npt) =  i
        elseif(jchk.ne.0) then
          call sqcDlele(subnam,'X(i)',xmic5,x(i),1.D0,' ')
          call sqcDlele(subnam,'Q(i)',qmic5,q(i),qmac5,' ')
        else
          f(i)        = qnull6
        endif
      enddo

      if(npt.eq.0) return

C--   Initialise list
      call sqcLstIni(yy,tt,npt,ww,11+45*mpt0,nused,ierr)
      if(ierr.eq.1) stop 'STFUNXQ Init: not enough space in ww'
      if(ierr.eq.2) stop 'STFUNXQ Init: no scratch buffer available'

C--   Fill buffer with function values
      call sqcFillBuffij(fun,ww,ierr)
      if(ierr.eq.1) stop 'STFUNXQ Fill: ww not initialised'
      if(ierr.eq.2) stop 'STFUNXQ Fill: evolution parameter change'
      if(ierr.eq.3) stop 'STFUNXQ Fill: no scratch buffer available'
      if(ierr.eq.4) stop 'STFUNXQ Fill: error from STFUN'

C--   Interpolate
      call sqcLstFun(ww,ff,mpt0,nout,ierr)
      if(ierr.eq.1) stop 'STFUNXQ LstF: ww not initialised'
      if(ierr.eq.2) stop 'STFUNXQ LstF: evolution parameter change'
      if(ierr.eq.3) stop 'STFUNXQ LstF: no buffer to interpolate'

C--   Update f
      do i = 1,nout
        f(ipoint(i)) = ff(i)
      enddo

      return
      end

