
C--   file usrfast.f containing fast convolution

C--   subroutine IdScope(w,idg)
C--   subroutine sqcUFBook(subnam)
C--   subroutine sqcUFIni(subnam,xlist,qlist,n,jchk)
C--   subroutine FastIni(xlist,qlist,n,jchk)
C--   subroutine FastClr(id)
C--   subroutine FastInp(w,jdf,jbuf)
C--   subroutine FastEpm(jset,id1,id2)
C--   subroutine FastSns(jset,qvec,isel,id)
C--   subroutine FastSum(iset,coef,id)
C--   subroutine FastFxK(w,idwt,id1,id2)
C--   subroutine FastFxF(w,idx,ida,idb,jbufo)
C--   subroutine FastKin(id,fun)
C--   subroutine FastCpy(id1,id2,iadd)
C--   subroutine FastFxq(id,stf,n)
C--
C--               -----------------
C--               Scoping in QCDNUM
C--               -----------------
C--
C--   The four scoping parameters in /qpars6/ are:
C--
C--   1. Lscopechek6  :   yes/no check the pdf scope (default no)
C--   2. iscopeslot6  :   slot in pars8 repository (default 1 = base)
C--   3. iscopeuse6   :   value of iscopeslot6 at start of scoping
C--   4. iscopekey6   :   scope key
C--
C--   Set scope by user call to IDSCOPE(w,iset)
C--   - Sets iscopeslot6 to key of iset, or to 1 (= base set = default)
C--
C--   Start scoping (by FASTINI, FASTCLR and STFUNXQ)
C--   - set Lscopechek6  = .true.
C--   - set iscopeuse6   = iscopeslot6
C--   - set iscopekey6   = key(iscopeslot6)
C--
C--   Check scope (all routines that access a pdf)
C--   - if Lscopechek6 is true then
C--       a. Check if pdf_key = iscopekey6 --> error if not
C--       b. Optionally check if iscopeslot6 = iscopeuse6 to detect an
C--          in-between change of slot by a call to IDSCOPE
C--
C--   Terminate scoping by FASTINP, FASTFXQ or call to IDSCOPE(w,-iset)
C--   - set Lscopechek6 = .false.

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM fast convolution routines from usrfast.          */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM fast convolution routines from usrfast.f           */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM fast convolution routines from usrfast.f           */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------


C==   ==================================================================
C==   Scoping ==========================================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void idscope(double *w, int jset);
C-----------------------------------------------------------------------
CXXHFW  #define fidscope FC_FUNC(idscope,IDSCOPE)
CXXHFW    void fidscope(double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void idscope(double *w, int jset)
CXXWRP  {
CXXWRP    fidscope(w, &jset);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================
      subroutine IdScope(w,jset)
C     ==========================

C--   Set the pdf scope

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      dimension w(*)

      logical lqcIsetExists

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'IDSCOPE ( W, ISET )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Release scoping
      if(jset.le.0) then
        Lscopechek6 = .false.
        return
      endif

      iscopeslot6 = 0
      kset        = abs(jset)

      if(int(w(1)).ne.654321) then
*mb        write(6,*) 'IDSCOPE internal memory iset = ',kset,isetf7(kset)
        call sqcIlele(subnam,'ISET',0,kset,mset0,' ')
        if(kset.eq.0) then
          iscopeslot6 = 1                                      !base set
        elseif(isetf7(kset).ne.0) then
          iscopeslot6 = int(dparGetPar(stor7,isetf7(kset),idipver8))
        else
          call sqcSetMsg(subnam,'ISET',kset)
        endif

      elseif(lqcIsetExists(w,kset)) then
*mb        write(6,*) 'IDSCOPE workspace kset = ',kset
        iscopeslot6 = int(dparGetPar(w,kset,idipver8))
      else
*mb        write(6,*) 'IDSCOPE workspace kset = ',kset
        call sqcSetMsg(subnam,'ISET',kset)
      endif

C--   Pdf set not filled
      if(iscopeslot6.eq.0) call sqcSetMsg(subnam,'ISET',kset)

      return
      end

C==   ==================================================================
C==   Fast convolutions        =========================================
C==   ==================================================================

C--   Remark:  The number of buffers nbuf is stored in mbf0
C--            The buffer index runs from 1 - nbuf

C     ============================
      subroutine sqcUFBook(subnam)
C     ============================

C--   Book scratch buffers
C--   Stop with error message if not enough space
C--   Called by FASTINI and FFLIST
C--
C--   subnam   (in) : subroutine name for error message

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      
      character*80 subnam 
      
C--   Book scratch buffers (or flag as empty if already booked)
      call sqcFastBook(nwords,ierr)
C--   Handle error code
      call sqcMemMsg(subnam,nwords,ierr)

      return
      
      end
      
C     ==============================================      
      subroutine sqcUFIni(subnam,xlist,qlist,n,jchk)
C     ==============================================

C--   Wrapper for sqcSetMark (setup interpolation) and clear buffers
C--   If jchk.ne.0 stop with error message when outside grid
C--   Called by FASTINI and FFLIST
C--
C--   subnam            (in) : s/r name for error message
C--   xlist(i),qlist(i) (in) : list of interpolation points
C--   n                 (in) : number of interpolation  points
C--   jchk              (in) : yes/no check points outside grid
C--
C--   The interpolation is setup by a call to sqcSetMark which fills
C--   the common block qfast9.inc

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'
      
      character*80 subnam
      
      dimension xlist(*),qlist(*)
      
      margin = 0
      call sqcSetMark(xlist,qlist,n,margin,ierr)
C--   At least one x,qmu2 outside grid
      if(jchk.eq.1 .and. ierr.eq.1) then
        call sqcErrMsg(subnam,'At least one x, mu2 outside grid')
      endif
C--   Clear buffers
      do i = 1,mbf0
         idg = iqcIdPdfLtoG(-1,i)
         call sqcPreset(idg,0.D0)
         isparse9(i) = 0          
      enddo

      return
      
      end

C-----------------------------------------------------------------------
CXXHDR    void fastini(double *xlist, double *qlist, int n, int jchk);
C-----------------------------------------------------------------------
CXXHFW  #define ffastini FC_FUNC(fastini,FASTINI)
CXXHFW    void ffastini(double*,double*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastini(double *xlist, double *qlist, int n, int jchk)
CXXWRP  {
CXXWRP    ffastini(xlist, qlist, &n, &jchk);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ======================================
      subroutine FastIni(xlist,qlist,n,jchk)
C     ======================================

C--   Pass list of x, mu2 values
C--   Also book scratch tables, if not already done
C--
C--   xlist = list of x values
C--   qlist = list of mu2 values
C--   n     = number of items in the list
C--   jchk  = 0 do not check xlist,qlist
C--         = 1 insist that all x,mu2 are inside grid   

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension xlist(*),qlist(*)
      
      character*80 subnam
      data subnam /'FASTINI ( X, QMU2, N, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check user input
      call sqcIlele(subnam,'N',1,n,mpt0,
     +'Please see the example program longlist.f to handle more points')
C--   Book scratch buffers (or flag as empty if already booked)
      call sqcUFBook(subnam)
C--   Setup list of interpolation points and clear buffers
      call sqcUFIni(subnam,xlist,qlist,n,jchk)
C--   Remember checking flag
      jchk9 = jchk
C--   Get key
      key = int(dparGetPar(pars8,iscopeslot6,idipver8))
C--   Store key
      iscopekey6  = key
      iscopeuse6  = iscopeslot6
      Lscopechek6 = .true.
C--   Point to correct parameters
      call sparParTo5(iscopekey6)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastclr(int ibuf);
C-----------------------------------------------------------------------
CXXHFW  #define ffastclr FC_FUNC(fastclr,FASTCLR)
CXXHFW    void ffastclr(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastclr(int ibuf)
CXXWRP  {
CXXWRP    ffastclr(&ibuf);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ========================
      subroutine FastClr(ibuf)
C     ========================

C--   Clear scratch buffer. If id = 0, clear all buffers  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'FASTCLR ( IBUF )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check identifier
      call sqcIlele(subnam,'IBUF',0,ibuf,mbf0,' ')
      if(ibuf.eq.0) then
        idmin = 1
        idmax = mbf0
      else
        idmin = ibuf
        idmax = ibuf
      endif
C--   Clear buffers
      jset = -1     !pdf set number assigned to fast buffers
      do i = idmin,idmax
         idg = iqcIdPdfLtoG(jset,i)
         call sqcPreset(idg,0.D0)
         isparse9(i) = 0       
      enddo
C--   For a clear-all, store current evolution parameters
      if(ibuf.eq.0) then
C--     Get key
        key = int(dparGetPar(pars8,iscopeslot6,idipver8))
C--     Store key
        iscopekey6  = key
        iscopeuse6  = iscopeslot6
        Lscopechek6 = .true.
       endif
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastinp(double *w, int jdf, double *coef, int jbuf, int iadd);
C-----------------------------------------------------------------------
CXXHFW  #define ffastinp FC_FUNC(fastinp,FASTINP)
CXXHFW    void ffastinp(double*,int*,double*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastinp(double *w, int jdf, double *coef, int jbuf, int iadd)
CXXWRP  {
CXXWRP    ffastinp(w, &jdf, coef, &jbuf, &iadd);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ========================================
      subroutine FastInp(w,jdf,coef,jbuf,iadd)
C     ========================================

C--   Copy weighted gluon or quark basis pdf to buffer jbuf
C--
C--   w         (in) : local workspace
C--   jdf       (in) : pdf indentifier in global format
C--   coef(3:6) (in) : nf-dependent weight
C--   jbuf      (in) : id of output scratch buffer [1,mbf0]
C--                    > 0 dense  buffer
C--                    < 0 sparse buffer
C--   iadd      (in) : store (0), add(1), substract(-1) to buffer

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /   -1,     5   / , iflg /   0,    1  /
      data icma  /    1,     5   /

      logical lint

      dimension w(*), coef(3:6)

      character*80 subnam
      data subnam /'FASTINP ( W, IDF, COEF, IBUF, IADD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check evolution parameters did not change
      call sqcFstMsg(subnam)

C--   Check jdf
      isi = 0
      igf = iqcSjekId(subnam,'IDF',w,jdf,icmi,icma,iflg,lint)

C--   Output identifier
      ibuf = abs(jbuf)      
      call sqcIlele(subnam,'IBUF',1,ibuf,mbf0,' ')

C--   Iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')

C--   Sparse or dense, thats the question
      if(jbuf.lt.0) then
        isparse9(ibuf) = 1     !sparse table
        idense         = 0
      elseif(isparse9(ibuf).eq.0 .or. iadd.eq.0) then
        isparse9(ibuf) = 2     !dense  table
        idense         = 1
      else
        idense         = isparse9(ibuf)-1
      endif

C--   Global buffer identifier
      ibg = iqcIdPdfLtoG(-1,ibuf)

C--   Do the work
      if(lint) then
        ksetw = igf/1000
        ipver = int(dparGetPar(stor7,ksetw,idipver8))
        call sparParTo5(ipver)
        call sqcFastInp(stor7,igf,ibg,iadd,coef,idense)
      else
        ksetw = igf/1000
        ipver = int(dparGetPar(w,ksetw,idipver8))
        call sparParTo5(ipver)
        call sqcFastInp(w    ,igf,ibg,iadd,coef,idense)
      endif

C--   Release scope check
      Lscopechek6 = .false.

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastepm(int jdum, int jdf, int jbuf);
C-----------------------------------------------------------------------
CXXHFW  #define ffastepm FC_FUNC(fastepm,FASTEPM)
CXXHFW    void ffastepm(int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastepm(int jdum, int jdf, int jbuf)
CXXWRP  {
CXXWRP    ffastepm(&jdum, &jdf, &jbuf);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      subroutine FastEpm(jdum,jdf,jbuf)
C     =================================

C--   Copy gluon or quark basis pdf to buffer jbuf
C--
C--   jdum  (in) : dummy variable
C--   jdf   (in) : pdf indentifier in global format
C--   jbuf  (in) : id of output scratch table [1,mbf0]
C--                > 0 dense  buffer
C--                < 0 sparse buffer

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /   -1,     5   / , iflg /   0,    1  /
      data icma  /   -1,     5   /

      logical lint

      dimension coef(3:6,0:12)

      character*80 subnam
      data subnam /'FASTEPM ( ISET, IDF, IBUF )'/

      jset = jdum !avoid compoler warning

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check jdf
      igf = iqcSjekId(subnam,'IDF',stor7,jdf,icmi,icma,iflg,lint)
C--   Check evolution parameters did not change
      call sqcFstMsg(subnam)
C--   Check jbuf identifier
      ibuf = abs(jbuf)      
      call sqcIlele(subnam,'IBUF',1,ibuf,mbf0,' ')
C--   Initialize
      isparse9(ibuf) = 0                !empty table
C--   Do the work
      call sqcIdPdfGtoL(igf,kset,idf)   !kset [1,mset0]  idf [0,12]
      do j = 3,6
        do i = 0,12
          coef(j,i) = 0.D0
        enddo
        coef(j,idf) = 1.D0
      enddo
C--   Get global id of gluon
      idg = iqcIdPdfLtoG(kset,0)
C--   Point to the right set
      ksetw = idg/1000
      ipver = int(dparGetPar(stor7,ksetw,idipver8))
      call sparParTo5(ipver)
C--   Get global id of ibuf
      ibg = iqcIdPdfLtoG(-1,ibuf)

C--   Do the linear combination
      if(jbuf.gt.0) then
        isparse9(ibuf) = 2                           !dense table
        call sqcFastPdf(idg,coef,ibg,1)
      else
        isparse9(ibuf) = 1                           !sparse table
        call sqcFastPdf(idg,coef,ibg,0)
      endif   
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastsns(int jset, double *qvec, int isel, int jbuf);
C-----------------------------------------------------------------------
CXXHFW  #define ffastsns FC_FUNC(fastsns,FASTSNS)
CXXHFW    void ffastsns(int*,double*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastsns(int jset, double *qvec, int isel, int jbuf)
CXXWRP  {
CXXWRP    ffastsns(&jset, qvec, &isel, &jbuf);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================================
      subroutine FastSns(jset,qvec,isel,jbuf)
C     =======================================

C--   Copy gluon or selected si/ns component to buffer jbuf
C--
C--   jset        (in) : pdf set identifier [1,mset0]
C--   qvec(-6:6)  (in) : lin combination of q, qbar
C--   isel        (in) : selection flag [0-7]
C--   jbuf        (in) : id of output buffer [1,mbf0]
C--                      > 0 dense  buffer
C--                      < 0 sparse buffer

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension qvec(-6:6),evec(12),coef(3:6,0:12)
      dimension mask(0:12,0:7)
C--              g s           v
C--              0 1 2 3 4 5 6 7 8 9 0 1 2      
      data mask /1,0,0,0,0,0,0,0,0,0,0,0,0,     !0=gluon
     +           0,1,0,0,0,0,0,0,0,0,0,0,0,     !1=singlet
     +           0,0,1,1,1,1,1,0,0,0,0,0,0,     !2=ns+
     +           0,0,0,0,0,0,0,1,0,0,0,0,0,     !3=valence
     +           0,0,0,0,0,0,0,0,1,1,1,1,1,     !4=ns-
     +           0,0,0,0,0,0,0,1,1,1,1,1,1,     !5=v and ns-
     +           0,0,1,1,1,1,1,1,1,1,1,1,1,     !6=all ns
     +           0,1,1,1,1,1,1,1,1,1,1,1,1   /  !7=all quarks

      character*80 subnam
      data subnam /'FASTSNS( ISET, DEF, ISEL, IBUF )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      ibuf = abs(jbuf)
C--   Check jset      
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Check status bits (also if jset exists)
      call sqcChkflg(jset,ichk,subnam)
C--   Check evolution parameters did not change
      call sqcFstMsg(subnam)
C--   Check jset filled with current parameters
      call sqcParMsg(subnam,'ISET',jset)
C--   Check selection flag
      call sqcIlele(subnam,'ISEL',0,isel,7,' ')
C--   Check ibuf identifier
      call sqcIlele(subnam,'IBUF',1,ibuf,mbf0,' ')
C--   Initialize
      isparse9(ibuf) = 0                    !empty table
C--   Do the work
      do nf = 3,6
        call sqcEfromQQ(qvec, evec, nf, nf) !quark coefficients
        coef(nf,0) = evec(1)*mask(0,isel)   !gluon = singlet * mask
        do i = 1,12
          coef(nf,i) = evec(i)*mask(i,isel) !apply selection mask quarks
        enddo
      enddo
C--   Get global id of gluon
      idg = iqcIdPdfLtoG(jset,0)
C--   Point to the right set
      ksetw = idg/1000
      ipver = int(dparGetPar(stor7,ksetw,idipver8))
      call sparParTo5(ipver)
C--   Get global identifier of ibuf
      ibg = iqcIdPdfLtoG(-1,ibuf)
C--   Do the linear combination
      if(jbuf.gt.0) then
        isparse9(ibuf) = 2                            !dense table
        call sqcFastPdf(idg,coef,ibg,1)
      else
        isparse9(ibuf) = 1                            !sparse table
        call sqcFastPdf(idg,coef,ibg,0)
      endif    

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastsum(int jset, double *coef, int jbuf);
C-----------------------------------------------------------------------
CXXHFW  #define ffastsum FC_FUNC(fastsum,FASTSUM)
CXXHFW    void ffastsum(int*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastsum(int jset, double *coef, int jbuf)
CXXWRP  {
CXXWRP    ffastsum(&jset, coef, &jbuf);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ==================================
      subroutine FastSum(jset,coef,jbuf)
C     ==================================

C--   Copy gluon or linear combination of quarks to buffer jbuf
C--
C--   jset            (in) : pdf set identifier [1,mset0]
C--   coef(0:12,3:6)  (in) : table of coefficients
C--   jbuf            (in) : id of output buffer [1,mbf0]
C--                          > 0 dense  buffer
C--                          < 0 sparse buffer

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension coef(0:12,3:6),coefn(3:6,0:12)

      character*80 subnam
      data subnam /'FASTSUM ( ISET, COEF, IBUF )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      ibuf = abs(jbuf)
C--   Check jset      
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ') 
C--   Check status bits (also if jset exists)
      call sqcChkflg(jset,ichk,subnam)
C--   Check evolution parameters did not change
      call sqcFstMsg(subnam)
C--   Check jset filled with current parameters
      call sqcParMsg(subnam,'ISET',jset)
C--   Check ibuf identifier
      call sqcIlele(subnam,'IBUF',1,ibuf,mbf0,' ')
C--   Initialize
      isparse9(ibuf) = 0                                !empty table
C--   Get global id of gluon
      idg = iqcIdPdfLtoG(jset,0)
C--   Point to the right set
      ksetw = idg/1000
      ipver = int(dparGetPar(stor7,ksetw,idipver8))
      call sparParTo5(ipver)
C--   Get global identifier of ibuf
      ibg = iqcIdPdfLtoG(-1,ibuf)
C--   Re-arrange coef
      do i = 0,12
        do j = 3,6
          coefn(j,i) = coef(i,j)
        enddo
      enddo
C--   Do the linear combination
      if(jbuf.gt.0) then
        isparse9(ibuf) = 2                              !dense table
        call sqcFastPdf(idg,coefn,ibg,1)
      else
        isparse9(ibuf) = 1                              !sparse table
        call sqcFastPdf(idg,coefn,ibg,0)
      endif 

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastfxk(double *w, int *idwt, int ibuf1, int jbuf2);
C-----------------------------------------------------------------------
CXXHFW  #define ffastfxk FC_FUNC(fastfxk,FASTFXK)
CXXHFW    void ffastfxk(double*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastfxk(double *w, int *idwt, int ibuf1, int jbuf2)
CXXWRP  {
CXXWRP    ffastfxk(w, idwt, &ibuf1, &jbuf2);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ======================================
      subroutine FastFxK(w,idwt,ibuf1,jbuf2)
C     ======================================

C--   Convolution F cross K
C--
C--   w       (in) : store filled with weight tables
C--   idwt(5) (in) : id(1)-(3) table ids LO, NLO, NNLO (0=no table)
C--                  id(4)     leading power of alfas  (0 or 1)
C--                  id(5)     not used
C--   ibuf1   (in) : scratch buffer w/pdfs previously filled by fastinp
C--   jbuf2   (in) : output scratch buffer
C--                  > 0 sparse buffer
C--                  < 0 sparse buffer

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     1   / , iflg /   1,    1  /
      data icma  /    1,     4   /

      logical lint

      dimension w(*),idwt(*),jdwt(5)

      character*80 subnam
      data subnam /'FASTFXK ( W, IDW, IBUF1, IBUF2 )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      ibuf2 = abs(jbuf2)
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check evolution parameters did not change
      call sqcFstMsg(subnam)
C--   Check weight table identifiers
      do i = 1,3
        jdwt(i) = iqcSjekId(subnam,'IDW',
     +                      w,idwt(i),icmi,icma,iflg,lint)
      enddo
      jdwt(4) = idwt(4)
C--   Check power of alphas
      call sqcIlele(subnam,'IDW(4)',0,idwt(4),1,' ')
C--   Where to get the alphas table from (aways base set)
      jdwt(5) = 1                                              !base set
C--   No overwrite, thank you
      if(ibuf1.eq.ibuf2) then
        call sqcErrMsg(subnam,'IBUF1 cannot be equal to IBUF2')
      endif
C--   Check identifiers           
      call sqcIlele(subnam,'IBUF1',1,ibuf1,mbf0,' ')
      call sqcIlele(subnam,'IBUF2',1,ibuf2,mbf0,' ')
C--   Check ibuf1 not empty or sparse
      if(isparse9(ibuf1).eq.0) then
        call sqcErrMsg(subnam,'IBUF1 empty buffer')
      endif  
C--   Check ibuf1 not sparse
      if(isparse9(ibuf1).eq.1) then
        call sqcErrMsg(subnam,'IBUF1 sparse buffer')
      endif
C--   Point to the correct parameters
      call sparParTo5(iscopekey6)
C--   Initialize
      isparse9(ibuf2) = 0                                   !empty table
C--   Global identifiers
      ibg1 = iqcIdPdfLtoG(-1,ibuf1)
      ibg2 = iqcIdPdfLtoG(-1,ibuf2)
C--   Do the work
      if(jbuf2.gt.0) then
        isparse9(ibuf2) = 1                                !sparse table
        call sqcFastFxK(w,jdwt,ibg1,ibg2,0,ierr)
      else
        isparse9(ibuf2) = 2                                 !dense table
        call sqcFastFxK(w,jdwt,ibg1,ibg2,1,ierr)
      endif 

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastfxf(double *w, int idx, int ibuf1, int ibuf2, int jbuf3);
C-----------------------------------------------------------------------
CXXHFW  #define ffastfxf FC_FUNC(fastfxf,FASTFXF)
CXXHFW    void ffastfxf(double*,int*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastfxf(double *w, int idx, int ibuf1, int ibuf2, int jbuf3)
CXXWRP  {
CXXWRP    ffastfxf(w, &idx, &ibuf1, &ibuf2, &jbuf3);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ===========================================
      subroutine FastFxF(w,idx,ibuf1,ibuf2,jbuf3)
C     ===========================================

C--   Convolution F cross F
C--
C--   w         (in) : store filled with weight tables
C--   idx       (in) : weight table filled by makewtx
C--   ibuf1,2   (in) : scratch buffers with pdfs
C--   jbuf3     (in) : output scratch buffer
C--                  > 0 sparse buffer
C--                  < 0 sparse buffer

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     1   / , iflg /   0,    1  /
      data icma  /    1,     4   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam /'FASTFXF ( W, IDX, IBUF1, IBUF2, IBUF3 )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      ibuf3 = abs(jbuf3)
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check evolution parameters did not change
      call sqcFstMsg(subnam)
C--   Check weight table identifier
      igx = iqcSjekId(subnam,'IDX',w,idx,icmi,icma,iflg,lint)
C--   No overwrite, thank you
      if(ibuf3.eq.ibuf1 .or. ibuf3.eq.ibuf2) then
        call sqcErrMsg(subnam,'IBUF3 cannot be equal to IBUF1 or IBUF2')
      endif
C--   Check identifiers           
      call sqcIlele(subnam,'IBUF1',1,ibuf1,mbf0,' ')
      call sqcIlele(subnam,'IBUF2',1,ibuf2,mbf0,' ')
      call sqcIlele(subnam,'IBUF3',1,ibuf3,mbf0,' ')
C--   Check ibuf1,2 not empty
      if(isparse9(ibuf1).eq.0) then
        call sqcErrMsg(subnam,'IBUF1 empty buffer')
      endif
      if(isparse9(ibuf2).eq.0) then
        call sqcErrMsg(subnam,'IBUF2 empty buffer')
      endif   
C--   Check ibuf1,2 not sparse
      if(isparse9(ibuf1).eq.1) then
        call sqcErrMsg(subnam,'IBUF1 sparse buffer')
      endif
      if(isparse9(ibuf2).eq.1) then
        call sqcErrMsg(subnam,'IBUF2 sparse buffer')
      endif
C--   Point to the correct parameters
      call sparParTo5(iscopekey6)
C--   Initialize
      isparse9(ibuf3) = 0                                !empty table
C--   Global identifiers
      ibg1 = iqcIdPdfLtoG(-1,ibuf1)
      ibg2 = iqcIdPdfLtoG(-1,ibuf2)
      ibg3 = iqcIdPdfLtoG(-1,ibuf3)
C--   Do the work
      if(jbuf3.gt.0) then
        isparse9(ibuf3) = 1                              !sparse table
        call sqcFastFxF(w,igx,ibg1,ibg2,ibg3,0)
      else
        isparse9(ibuf3) = 2                              !dense table
        call sqcFastFxF(w,igx,ibg1,ibg2,ibg3,1)
      endif 

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastkin(int ibuf, double (*fun)(int*,int*,int*,int*));
C-----------------------------------------------------------------------
CXXHFW  #define ffastkin FC_FUNC(fastkin,FASTKIN)
CXXHFW    void ffastkin(int*,double(*)(int*,int*,int*,int*));
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastkin(int ibuf, double (*fun)(int*,int*,int*,int*))
CXXWRP  {
CXXWRP    ffastkin(&ibuf, fun);
CXXWRP  }
C-----------------------------------------------------------------------

C     ============================
      subroutine FastKin(ibuf,fun)
C     ============================

C--   Multiply contents of pdf table id by fun(ix,iq)
C--
C--   Input:   ibuf = buffer filled with convolutions
C--            fun  = user defined function fun(ix,iq,nf,ithresh)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      external fun

      character*80 subnam
      data subnam /'FASTKIN ( IBUF, FUN )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check evolution parameters did not change
      call sqcFstMsg(subnam)
C--   Check identifier
      call sqcIlele(subnam,'IBUF',1,ibuf,mbf0,' ')
C--   Check id not empty
      if(isparse9(ibuf).eq.0) then
        call sqcErrMsg(subnam,'IBUF empty buffer')
      endif       
C--   Point to the correct parameters
      call sparParTo5(iscopekey6)
C--   Global identifier
      ibg = iqcIdPdfLtoG(-1,ibuf)
C--   Do the work
      if(isparse9(ibuf).eq.1) then
        call sqcFastKin(ibg,fun,0)
      else
        call sqcFastKin(ibg,fun,1)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastcpy(int ibuf1, int ibuf2, int iadd);
C-----------------------------------------------------------------------
CXXHFW  #define ffastcpy FC_FUNC(fastcpy,FASTCPY)
CXXHFW    void ffastcpy(int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastcpy(int ibuf1, int ibuf2, int iadd)
CXXWRP  {
CXXWRP    ffastcpy(&ibuf1, &ibuf2, &iadd);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ====================================
      subroutine FastCpy(ibuf1,ibuf2,iadd)
C     ====================================

C--   Add contents of id1 to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'FASTCPY ( IBUF1, IBUF2, IADD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check evolution parameters did not change
      call sqcFstMsg(subnam)
C--   Check identifiers
      if(ibuf1.eq.ibuf2) then
        call sqcErrMsg(subnam,'IBUF1 cannot be equal to IBUF2')
      endif
      call sqcIlele(subnam,'IBUF1',1,ibuf1,mbf0,' ')
      call sqcIlele(subnam,'IBUF2',1,ibuf2,mbf0,' ')
C--   Check id1 not empty
      if(isparse9(ibuf1).eq.0) then
        call sqcErrMsg(subnam,'IBUF1 empty buffer')
      endif       
C--   Check iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')
C--   Sparse or dense output thats the question
      if(isparse9(ibuf2).eq.0 .or. iadd.eq.0) then
        isparse9(ibuf2) = isparse9(ibuf1)
      else  
        isparse9(ibuf2) = min(isparse9(ibuf1),isparse9(ibuf2))
      endif

C--   Global buffer identifiers
      ibg1 = iqcIdPdfLtoG(-1,ibuf1)
      ibg2 = iqcIdPdfLtoG(-1,ibuf2)

C--   Point to the correct parameters
      call sparParTo5(iscopekey6)

C--   Do the work
      call sqcFastCpy(ibg1,ibg2,iadd,isparse9(ibuf2)-1)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fastfxq(int ibuf, double *stf, int n);
C-----------------------------------------------------------------------
CXXHFW  #define ffastfxq FC_FUNC(fastfxq,FASTFXQ)
CXXHFW    void ffastfxq(int*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fastfxq(int ibuf, double *stf, int n)
CXXWRP  {
CXXWRP    ffastfxq(&ibuf, stf, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================
      subroutine FastFxq(ibuf,stf,n)
C     ==============================

C--   Interpolation
C--
C--   Input:   ibuf  = scratch table filled with structure function
C--            xlst9 = list of x values in /qfast9/
C--            qlst9 = list of qmu2 values in /qfast9/ 
C--            n     = number of interpolations requested = min(n,nlst9)
C--
C--   Output:  stf   = array of interpolated stf values

C--   NB: the list of x,qmu2 values is entered via s/r fastini

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension stf(*)

      character*80 subnam
      data subnam /'FASTFXQ ( IBUF, F, N )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check evolution parameters do not change
      call sqcFstMsg(subnam)
C--   Check identifier
      call sqcIlele(subnam,'IBUF',1,ibuf,mbf0,' ')
C--   Check ibuf not empty
      if(isparse9(ibuf).eq.0) then
        call sqcErrMsg(subnam,'IBUF empty buffer')
      endif       
      call sqcIlele(subnam,'N',1,n,mpt0,
     +'Please see the example program longlist.f to handle more points')
C--   Point to the correct parameters
      call sparParTo5(iscopekey6)

C--   Do the work: pdf identifier is buffer identifier
      idg = iqcIbufGlobal(ibuf)
      call sqcFastFxq(stor7,idg,stf,n)

C--   Release scope check
      Lscopechek6 = .false.

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
