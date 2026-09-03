
C--   This is the file usrpdf.f with user pdf output routines

C--   subroutine AllFxq(jset,xx,qq,pdf,n,jchk)
C--   subroutine AllFij(jset,ix,iq,pdf,n,jchk)
C--   function   Bvalxq(jset,id,xx,qq,jchk)
C--   function   Bvalij(jset,id,ix,iq,jchk)
C--   function   Fvalxq(jset,xx,qq,jchk)
C--   function   Fvalij(jset,id,ix,iq,jchk)
C--   function   Sumfxq(jset,def,isel,xx,qq,jchk)
C--   function   Sumfij(jset,def,isel,ix,iq,jchk)
C--   function   fsplne(jset,id,xx,iq)
C--   function   splchk(jset,id,iq)
C--
C--   subroutine FFList(jset,w,isel,x,q,f,n,jchk)
C--   subroutine FFTabl(jset,w,isel,x,nx,q,nq,table,m,jchk)
C--   subroutine FTable(jset,w,isel,x,nx,q,nq,table,jchk)
C--   subroutine ffplotCPP(file,ls1,fun,m,zmi,zma,n,txt,ls2)
C--   subroutine FFPlot(file,fun,m,zmi,zma,n,txt)
C--   subroutine fiplotCPP(file,ls1,fun,m,zval,n,txt,ls2)
C--   subroutine FIPlot(file,fun,m,zval,n,txt)
C--
C--   subroutine sqcPdfLims(idg,iy1,iy2,it1,it2,jmax)
C--   subroutine sqcParForSumPdf(jset,coef,isel,par,n,nout,ierr)
C--   integer function iqcGimmeScratch()
C--   subroutine sqcReleaseScratch(idg)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM pdf output routines from usrpdf.f                */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM pdf output routines from usrpdf.f                  */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM pdf output routines from usrpdf.f                  */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------


C==   ==================================================================
C==   Pdf interpolation  ===============================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void allfxq(int iset, double x, double qmu2, double *pdfs, int n, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fallfxq FC_FUNC(allfxq,ALLFXQ)
CXXHFW    void fallfxq(int*, double*, double*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void allfxq(int iset, double x, double qmu2, double *pdfs, int n, int ichk)
CXXWRP  {
CXXWRP    fallfxq(&iset,&x,&qmu2,pdfs,&n,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------

C     ========================================
      subroutine AllFxq(jset,xx,qq,pdf,n,jchk)
C     ========================================

C--   Get all flavour pdfs of iset
C--
C--   jset        (in)  : pdf set [1,mset0]
C--   pdf(-6:6+n) (out) : the whole lot...
C--   n           (in)  : number of extra pdfs to store
C--   jchk        (in)  : [-1,1] check input and/or grid boundary
C--
C--                        input    limit
C--           jchk =  1    yes      yes
C--                =  0    yes      no
C--                = -1    no       yes

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension pdf(-6:6+n)

      character*80 subnam
      data subnam /'ALLFXQ ( ISET, X, Q2, PDF, N, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      if(jchk.ne.-1) then
C--     Check pdf set
        call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--     Check status bits (also that pdf set is filled)
        call sqcChkflg(jset,ichk,subnam)
C--     Check pdf set is filled with current parameters
        call sqcParMsg(subnam,'ISET',jset)
C--     Check number of extra pdfs
        call sqcIlele(subnam,'N',0,n,ilast7(jset)-12,
     +               'Attempt to access nonexisting extra pdfs in ISET')
C--     Check pdf type
        if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALXQ instead')
      endif

C--   Null by default
      do i = -6,6+n
        pdf(i) = qnull6
      enddo

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      yy = dqcXInside(subnam,xx,jchk)
      if(yy.eq.-1.D0) return
C--   Catch x = 1
      if(yy.eq.0.D0) then
        do i = -6,6+n
          pdf(i) = 0.D0
        enddo
        return
      endif
C--   Check if inside q-grid or cuts
      tt = dqcQInside(subnam,qq,jchk)
      if(tt.eq.0.D0) return
C--   Go ...
      idg   = iqcIdPdfLtoG(jset,0)
      call sqcAllFyt(idg,yy,tt,pdf,n)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void allfij(int iset, int ix, int iq, double *pdfs, int n, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fallfij FC_FUNC(allfij,ALLFIJ)
CXXHFW    void fallfij(int*, int*, int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void allfij(int iset, int ix, int iq, double *pdfs, int n, int ichk)
CXXWRP  {
CXXWRP    fallfij(&iset,&ix,&iq,pdfs,&n,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------
        
C     ========================================
      subroutine AllFij(jset,ix,jq,pdf,n,jchk)
C     ========================================

C--   Get all flavour pdfs of iset
C--
C--   jset        (in)  : pdf set [1,mset0]
C--   pdf(-6:6+n) (out) : the whole lot...
C--   n           (in)  : number of extra pdfs to store
C--   jchk        (in)  : [-1,1] check input and/or grid boundary
C--
C--                        input    grid
C--           jchk =  1      y        y
C--                =  0      y        n
C--                = -1      n        y
C--
C--   Remark: jq > 0 then follow the 4,5,6 convention at iqc,b,t
C--              < 0 then follow the 3,4,5 convention at iqc,b,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension pdf(-6:6+n)

      character*80 subnam
      data subnam /'ALLFIJ ( ISET, IX, IQ, PDF, N, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      if(jchk.ne.-1) then
C--     Check pdf set
        call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--     Check status bits (also if pdf set is filled)
        call sqcChkflg(jset,ichk,subnam)
C--     Check pdf set is filled with current parameters
        call sqcParMsg(subnam,'ISET',jset)
C--     Check number of extra pdfs
        call sqcIlele(subnam,'N',0,n,ilast7(jset)-12,
     +               'Attempt to access nonexisting extra pdfs in ISET')
C--     Check pdf type
        if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALIJ instead')
      endif

C--   Null by default
      do i = -6,6+n
        pdf(i) = qnull6
      enddo

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      iy = iqcIxInside(subnam,ix,jchk)
      if(iy.eq.-1) return
C--   Catch x = 1
      if(iy.eq.0) then
        do i = -6,6+n
          pdf(i) = 0.D0
        enddo
        return
      endif
C--   Check if inside q-grid or cuts
      it = iqcIqInside(subnam,jq,jchk)
      if(it.eq.0) return
C--   Go ...
      idg = iqcIdPdfLtoG(jset,0)
      call sqcAllFij(idg,iy,it,pdf,n)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double bvalxq(int iset, int id, double x, double qmu2, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fbvalxq FC_FUNC(bvalxq,BVALXQ)
CXXHFW    double fbvalxq(int*, int*, double*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double bvalxq(int iset, int id, double x, double qmu2, int ichk)
CXXWRP  {
CXXWRP    return fbvalxq(&iset,&id,&x,&qmu2,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ====================================================
      double precision function Bvalxq(jset,id,xx,qq,jchk)
C     ====================================================

C--   Get basis pdf stored in jset versus x,q2
C--
C--   jset  (in) : pdf set [1,mset0]
C--   id    (in) : g, si, ns identifier [0,12+n]
C--   jchk  (in) : [-1,1] check input and/or grid boundary
C--
C--                        input    grid
C--           jchk =  1      y        y
C--                =  0      y        n
C--                = -1      n        y

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

*      logical lmb_eq

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'BVALXQ ( ISET, ID, X, Q2, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      if(jchk .ne. -1) then
C--     Check jset is in range
        call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--     Check status bits (also check if pdf set is filled)
        call sqcChkflg(jset,ichk,subnam)
C--     Check id is in range
        call sqcIlele(subnam,'ID',0,id,ilast7(jset),' ')
C--     Check pdf set is filled with current parameters
        call sqcParMsg(subnam,'ISET',jset)
      endif

      Bvalxq = qnull6   !null by default

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      yy = dqcXInside(subnam,xx,jchk)
      if(yy.eq.-1.D0) return
C--   Catch x = 1
      if(yy.eq.0.D0) then
        Bvalxq = 0.D0
        return
      endif
C--   Check if inside q-grid or cuts
      tt = dqcQInside(subnam,qq,jchk)
      if(tt.eq.0.D0) return
C--   Go ...
      idg    = iqcIdPdfLtoG(jset,id)
      Bvalxq = dqcBvalyt(idg,yy,tt)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double bvalij(int iset, int id, int ix, int iq, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fbvalij FC_FUNC(bvalij,BVALIJ)
CXXHFW    double fbvalij(int*, int*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double bvalij(int iset, int id, int ix, int iq, int ichk)
CXXWRP  {
CXXWRP    return fbvalij(&iset,&id,&ix,&iq,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ====================================================
      double precision function Bvalij(jset,id,ix,jq,jchk)
C     ====================================================

C--   Get basis pdf stored in jset versus ix,iq
C--
C--   jset  (in) :  pdf set [1,mset0]
C--   id    (in) :  g, si, ns identifier [0,12+n]
C--   jchk  (in) :  [-1,1] check input and/or grid boundary
C--
C--                        input    grid
C--           jchk =  1      y        y
C--                =  0      y        n
C--                = -1      n        y
C--
C--   Remark: jq > 0 then follow the 4,5,6 convention at iqc,b,t
C--              < 0 then follow the 3,4,5 convention at iqc,b,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'BVALIJ ( ISET, ID, IX, IQ, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      if(jchk .ne. -1) then
C--     Check pdf set
        call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--     Check status bits  (also check if pdf set is filled)
        call sqcChkflg(jset,ichk,subnam)
C--     Check identifier
        call sqcIlele(subnam,'ID',0,id,ilast7(jset),' ')
C--     Check pdf set is filled with current parameters
        call sqcParMsg(subnam,'ISET',jset)
      endif

      Bvalij = qnull6                   !null by default

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      iy = iqcIxInside(subnam,ix,jchk)
      if(iy.eq.-1) return
C--   Catch x = 1
      if(iy.eq.0) then
        Bvalij = 0.D0
        return
      endif
C--   Check if inside q-grid or cuts
      it = iqcIqInside(subnam,jq,jchk)
      if(it.eq.0) return
C--   Go ...
      idg    = iqcIdPdfLtoG(jset,id)
      Bvalij = dqcBvalij(idg,iy,it)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double fvalxq(int iset, int id, double x, double qmu2, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define ffvalxq FC_FUNC(fvalxq,FVALXQ)
CXXHFW    double ffvalxq(int*, int*, double*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double fvalxq(int iset, int id, double x, double qmu2, int ichk)
CXXWRP  {
CXXWRP    return ffvalxq(&iset,&id,&x,&qmu2,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ====================================================
      double precision function Fvalxq(jset,id,xx,qq,jchk)
C     ====================================================

C--   Get gluon, quark, antiquark pdf
C--
C--   jset  (in) :  pdf set [1,mset0]
C--   id    (in) :  gluon, q, qbar identifier [-6,6+n]
C--   jchk  (in) :  [-1,1] check input and/or grid boundary
C--
C--                        input    grid
C--           jchk =  1      y        y
C--                =  0      y        n
C--                = -1      n        y

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'FVALXQ ( ISET, ID, X, Q2, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      if(jchk.ne.-1) then
C--     Check pdf set
        call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--     Check status bits  (also check if pdf set is filled)
        call sqcChkflg(jset,ichk,subnam)
C--     Check identifier
        call sqcIlele(subnam,'ID',-6,id,ilast7(jset)-6,' ')
C--     Check pdf set is filled with current parameters
        call sqcParMsg(subnam,'ISET',jset)
C--     Check pdf type
        if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALXQ instead')
      endif

      Fvalxq = qnull6   !null by default

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      yy = dqcXInside(subnam,xx,jchk)
      if(yy.eq.-1.D0) return
      if(yy.eq.0.D0) then
        Fvalxq = 0.D0
        return
      endif
C--   Check if inside q-grid or cuts
      tt = dqcQInside(subnam,qq,jchk)
      if(tt.eq.0.D0) return
C--   Go ...
      idg    = iqcIdPdfLtoG(jset,0)
      Fvalxq = dqcFvalyt(idg,id,yy,tt)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double fvalij(int iset, int id, int ix, int iq, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define ffvalij FC_FUNC(fvalij,FVALIJ)
CXXHFW    double ffvalij(int*, int*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double fvalij(int iset, int id, int ix, int iq, int ichk)
CXXWRP  {
CXXWRP    return ffvalij(&iset,&id,&ix,&iq,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ====================================================
      double precision function Fvalij(jset,id,ix,jq,jchk)
C     ====================================================

C--   Get gluon, quarks or antiquarks
C--
C--   jset  (in) :  pdf set [1,mset0]
C--   id    (in) :  gluon, q, qbar identifier [-6,6+n]
C--   jchk  (in) :  [-1,1] check input and/or grid boundary
C--
C--                        input    grid
C--           jchk =  1      y        y
C--                =  0      y        n
C--                = -1      n        y
C--
C--   Remark: jq > 0 then follow the 4,5,6 convention at iqc,b,t
C--              < 0 then follow the 3,4,5 convention at iqc,b,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'FVALIJ ( ISET, ID, IX, IQ, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      if(jchk .ne. -1) then
C--     Check pdf set
        call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--     Check status bits  (also check if pdf set is filled)
        call sqcChkflg(jset,ichk,subnam)
C--     Check identifier
        call sqcIlele(subnam,'ID',-6,id,ilast7(jset)-6,' ')
C--     Check pdf set is filled with current parameters
        call sqcParMsg(subnam,'ISET',jset)
C--     Check pdf type
        if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALIJ instead')
      endif

      Fvalij = qnull6   !null by default

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      iy = iqcIxInside(subnam,ix,jchk)
      if(iy.eq.-1) return
C--   Catch x = 1
      if(iy.eq.0) then
        Fvalij = 0.D0
        return
      endif
C--   Check if inside q-grid or cuts
      it = iqcIqInside(subnam,jq,jchk)
      if(it.eq.0) return
C--   Go ...
      idg    = iqcIdPdfLtoG(jset,0)
      Fvalij = dqcFvalij(idg,id,iy,it)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double sumfxq(int iset, double *c, int isel, double x, double qmu2, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fsumfxq FC_FUNC(sumfxq,SUMFXQ)
CXXHFW    double fsumfxq(int*, double*, int*, double*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double sumfxq(int iset, double *c, int isel, double x, double qmu2, int ichk)
CXXWRP  {
CXXWRP    return fsumfxq(&iset,c,&isel,&x,&qmu2,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================================================
      double precision function Sumfxq(jset,def,isel,xx,qq,jchk)
C     ==========================================================

C--   Weighted sum of quarks and antiquarks
C--
C--   jset       (in) :  pdf set [1,mset0]
C--   def(-6:6)  (in) :  coefficients of the linear combination
C--   jchk       (in) :  [-1,1] check input and/or grid boundary
C--
C--                        input    grid
C--           jchk =  1      y        y
C--                =  0      y        n
C--                = -1      n        y
C--
C--   NB: def(0) is ignored since it corresponds to the gluon

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension def(-6:6)

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'SUMFXQ ( ISET, C, ISEL, X, Q2, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      if(jchk.ne.-1) then
C--     Check pdf set
        call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--     Check status bits  (also check if pdf set is filled)
        call sqcChkflg(jset,ichk,subnam)
C--     Check pdf set is filled with current parameters
        call sqcParMsg(subnam,'ISET',jset)
C--     Check pdf type
        if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALXQ instead')
      endif

C--   Check isel
      if(isel.ge.13 .and. ilast7(jset).ge.13) then
        call sqcIlele(subnam,'ISEL',13,isel,ilast7(jset),' ')
      else
        call sqcIlele(subnam,'ISEL',0,isel,9,' ')
      endif

      Sumfxq = qnull6   !null by default

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      yy = dqcXInside(subnam,xx,jchk)
      if(yy.eq.-1.D0) return
      if(yy.eq.0.D0) then
        Sumfxq = 0.D0
        return
      endif
C--   Check if inside q-grid or cuts
      tt = dqcQInside(subnam,qq,jchk)
      if(tt.eq.0.D0) return
C--   Go ...
      idg    = iqcIdPdfLtoG(jset,0)
      Sumfxq = dqcFsumyt(idg,def,isel,yy,tt)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double sumfij(int iset, double *c, int isel, int ix, int iq, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fsumfij FC_FUNC(sumfij,SUMFIJ)
CXXHFW    double fsumfij(int*, double*, int*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double sumfij(int iset, double *c, int isel, int ix, int iq, int ichk)
CXXWRP  {
CXXWRP    return fsumfij(&iset,c,&isel,&ix,&iq,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ==========================================================
      double precision function Sumfij(jset,def,isel,ix,jq,jchk)
C     ==========================================================

C--   Weighted sum of quarks and antiquarks
C--
C--   jset       (in) :  pdf set [1,mset0]
C--   def(-6:6)  (in) :  coefficients of the linear combination
C--   isel       (in) :  singlet-nonsinglet selection parameter
C--   jchk       (in) :  [-1,1] check input and/or grid boundary
C--
C--                        input    grid
C--           jchk =  1      y        y
C--                =  0      y        n
C--                = -1      n        y
C--
C--   NB: def(0) is ignored since it corresponds to the gluon
C--
C--   Remark: jq > 0 then follow the 4,5,6 convention at iqc,b,t
C--              < 0 then follow the 3,4,5 convention at iqc,b,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension def(-6:6)
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'SUMFIJ ( ISET, C, ISEL, IX, IQ, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      if(jchk .ne. -1) then
C--     Check pdf set
        call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--     Check status bits  (also check if pdf set is filled)
        call sqcChkflg(jset,ichk,subnam)
C--     Check pdf set is filled with current parameters
        call sqcParMsg(subnam,'ISET',jset)
C--     Check pdf type
        if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALIJ instead')
      endif

C--   Check isel
      if(isel.ge.13 .and. ilast7(jset).ge.13) then
        call sqcIlele(subnam,'ISEL',13,isel,ilast7(jset),' ')
      else
        call sqcIlele(subnam,'ISEL',0,isel,9,' ')
      endif

      Sumfij = qnull6   !null by default

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      iy = iqcIxInside(subnam,ix,jchk)
      if(iy.eq.-1) return
C--   Catch x = 1
      if(iy.eq.0) then
        Sumfij = 0.D0
        return
      endif
C--   Check if inside q-grid or cuts
      it = iqcIqInside(subnam,jq,jchk)
      if(it.eq.0) return
C--   Go ...
      idg    = iqcIdPdfLtoG(jset,0)
      Sumfij = dqcFsumij(idg,def,isel,iy,it)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double fsplne(int iset, int id, double x, int iq);
C-----------------------------------------------------------------------
CXXHFW  #define ffsplne FC_FUNC(fsplne,FSPLNE)
CXXHFW    double ffsplne(int*, int*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double fsplne(int iset, int id, double x, int iq)
CXXWRP  {
CXXWRP    return ffsplne(&iset,&id,&x,&iq);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ===============================================
      double precision function fsplne(jset,id,xx,iq)
C     ===============================================

C--   Spline interpolation in x
C--
C--   jset   (in) : pdf set [1,mset0]
C--   id     (in) : pdf identifier [ifrst7,ilast7]
C--   xx     (in) : value of x
C--   iq     (in) : iq grid point

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'FSPLNE ( ISET, ID, X, IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check pdf set
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Check status bits  (also check if pdf set is filled)
      call sqcChkflg(jset,ichk,subnam)
C--   Check pdf set is filled
      if(.not.Lfill7(jset)) call sqcSetMsg(subnam,'ISET',jset)
C--   Check id
      call sqcIlele(subnam,'ID',ifrst7(jset),id,ilast7(jset),' ')

      fsplne = qnull6   !null by default

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Check if inside x-grid or cuts
      yy = dqcXInside(subnam,xx,1)
C--   Check if inside q-grid or cuts
      it = iqcIqInside(subnam,iq,1)
C--   Go ...
      idg    = iqcIdPdfLtoG(jset,id)
      fsplne = dqcXSplne(idg,yy,it)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double splchk(int iset, int id, int iq);
C-----------------------------------------------------------------------
CXXHFW  #define fsplchk FC_FUNC(splchk,SPLCHK)
CXXHFW    double fsplchk(int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double splchk(int iset, int id, int iq)
CXXWRP  {
CXXWRP    return fsplchk(&iset,&id,&iq);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ============================================
      double precision function splchk(jset,id,jq)
C     ============================================

C--   Check for spline oscillations
C--
C--   jset  (in)  Pdf set [1,mset0]
C--   id    (in)  Si/ns pdf identifier [0-12]
C--   jq    (in)  mu2 grid point

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qpdfs7.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'SPLCHK ( ISET, ID, IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Default value
      splchk = 0.D0
C--   Check pdf set
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Check status bits  (also check if pdf set is filled)
      call sqcChkflg(jset,ichk,subnam)
C--   Check pdf set is filled
      if(.not.Lfill7(jset)) call sqcSetMsg(subnam,'ISET',jset)
C--   Check user input
      call sqcIlele(subnam,'ID',ifrst7(jset),id,ilast7(jset),' ')
C--   Point to correct set
      call sparParTo5(ikeyf7(jset))
C--   Check if inside q-grid or cuts
      it = iqcIqInside(subnam,jq,0)
      if(it.eq.0) return
C--   Go ...
      idg    = iqcIdPdfLtoG(jset,id)
      splchk = dqcSplChk(idg,it)

      return
      end
      
C==   ==================================================================
C==   Lists, tables and plots ==========================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR  void fflist(int iset, double *c, int m, double *x, double *q, double *f, int n, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define ffflist FC_FUNC(fflist,FFLIST)
CXXHFW   void ffflist(int*, double*, int*, double*, double*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fflist(int iset, double *c, int m, double *x, double *q, double *f, int n, int ichk)
CXXWRP  {
CXXWRP    ffflist(&iset,c,&m,x,q,f,&n,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------
            
C     ===========================================
      subroutine FFList(jset,c,isel,x,q,f,n,jchk)
C     ===========================================

C--   Process list of interpolation points.
C--   The number of interpolations can be larger than mpt0 since
C--   this routine autmatically buffers in chunks of mpt0 words.
C--
C--   jset       (in) :  pdf set [1,mset0]
C--   c(-6:6+m)  (in) :  coefficients of a linear combination
C--   isel       (in) :  selection flag
C--   x          (in) :  list of x-points
C--   q          (in) :  list of mu2 points
C--   f         (out) :  list of interpolated pdfs
C--   n          (in) :  number of items in the list
C--   jchk       (in) :  0/1  no/yes check grid boundary
C--
C--   Layout of par(8 + (13+m)*8)
C--   ---------------------------
C--   par(1)  =  Karr(0)
C--   par(2)  =  Karr(1)
C--   par(3)  =  Karr(2)
C--   par(4)  =  Karr(3)
C--   par(5)  =  # pdfs for nf = 3
C--   par(6)  =  # pdfs for nf = 4
C--   par(7)  =  # pdfs for nf = 5
C--   par(8)  =  # pdfs for nf = 6
C--
C--   par(9)  =  first word of Parr(i,j,k)
C--   ..          ..
C--
C--   Layout of Parr(2,13+m,3:6)
C--   --------------------------
C--   Parr(1,j,nf) = base address of pdf j at nf
C--   Parr(2,j,nf) = coefficient  of pdf j at nf

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      
      dimension c(-6:6), par(8+8*mpdf0)
      dimension x(*), q(*), f(*)

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'FFLIST ( ISET, C, ISEL, X, Q, F, N, ICHK )'/

C--   Inline address function in par
      iaP(i,j,k) = int(par(1))+int(par(2))*i+int(par(3))*j+int(par(4))*k

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check pdf set
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Check status bits  (also check if pdf set is filled)
      call sqcChkflg(jset,ichk,subnam)
C--   Check pdf set is filled with current parameters
      call sqcParMsg(subnam,'ISET',jset)
C--   Check selection parameter
      if(isel.le.12) then
        call sqcIlele(subnam,'ISEL',0,isel,9,
     +                                    'Invalid selection parameter')
      else
        call sqcIlele(subnam,'ISEL',13,isel,ilast7(jset),
     +                'Attempt to access nonexisting extra pdf in ISET')
      endif
C--   Check number of interpolation points
      if(n.le.0) call sqcErrMsg(subnam,'N should be larger than zero')
C--   Check pdf type
      if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALXQ instead')

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Setup the par array for sqcPdfLstMpt
      call sqcParForSumPdf(jset,c,isel,par,8+8*mpdf0,np,kerr)

      nmax  = mpt0
*      nmax  = 10
      nlast = 0
      ntodo = min(n,nmax)
      do while( ntodo.gt.0 )
        i1    = nlast+1
        call sqcPdfLstMpt(
     +       subnam,par,idum,x(i1),q(i1),f(i1),ntodo,jchk)
        nlast = nlast+ntodo
        ntodo = min(n-nlast,nmax)
      enddo

      return
      end

C     ===================================================
      subroutine sqcPdfLstMpt(subnam,par,np,x,q,f,n,jchk)
C     ===================================================

C--   Fast interpolation of max mpt0 points
C--
C--   subnam   (in) : subroutine name for error message
C--   par      (in) : parameters passed to fun
C--   np       (in) : number of parameters
C--   x,q      (in) : list of interpolation points
C--   f       (out) : list of interpolated results
C--   n        (in) : number of items in x,q,f <= mpt0
C--   jchk     (in) : error checking flag

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical lmb_eq

      external dqcPdfSum

      dimension par(*), x(*), q(*), f(*)
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
      if(ierr.eq.1) stop 'FFLIST Init: not enough space in ww'
      if(ierr.eq.2) stop 'FFLIST Init: no scratch buffer available'

C--   Fill buffer with function values
      call sqcFillBuffer(dqcPdfSum,stor7,par,np,ww,ierr,jerr)
      if(ierr.eq.1) stop 'FFLIST Fill: ww not initialised'
      if(ierr.eq.2) stop 'FFLIST Fill: evolution parameter change'
      if(ierr.eq.3) stop 'FFLIST Fill: no scratch buffer available'
      if(ierr.eq.4) stop 'FFLIST Fill: error from dqcPdfSum'

C--   Interpolate
      call sqcLstFun(ww,ff,mpt0,nout,ierr)
      if(ierr.eq.1) stop 'FFLIST LstF: ww not initialised'
      if(ierr.eq.2) stop 'FFLIST LstF: evolution parameter change'
      if(ierr.eq.3) stop 'FFLIST LstF: found no buffer to interpolate'

C--   Update f
      do i = 1,nout
        f(ipoint(i)) = ff(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR  void fftabl(int iset, double *c, int isel, double *x, int nx, double *q, int nq, double *table, int m, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define ffftabl FC_FUNC(fftabl,FFTABL)
CXXHFW   void ffftabl(int*, double*, int*, double*, int*, double*, int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fftabl(int iset, double *c, int isel, double *x, int nx, double *q, int nq, double *table, int m, int ichk)
CXXWRP  {
CXXWRP    ffftabl(&iset,c,&isel,x,&nx,q,&nq,table,&m,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------

C     =====================================================
      subroutine FFTabl(jset,c,isel,xx,nx,qq,nq,pdf,m,jchk)
C     =====================================================

C--   Pdf linear combination interpolated on a grid
C--
C--   jset       (in) :  pdf set [1,mset0]
C--   c(-6:6)    (in) :  coefficients of a linear combination
C--   isel       (in) :  number of extra pdfs beyond gluon and quarks
C--   xx         (in) :  input array of x-values
C--   nx         (in) :  number of x-values
C--   qq         (in) :  input array of mu2-values
C--   nq         (in) :  number of mu2-values
C--   pdf       (out) :  pdf(m,n>=nq) 2-dim array of interpolated pdfs
C     m          (in) :  first dimension of pdf
C--   jchk       (in) :  # 0: error if xx,qq out of range
C--                      = 0: return null if xx,qq out of range

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical lmb_eq

      external dqcPdfSum

      dimension c(-6:6), par(8+8*mpdf0)

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)
      
      dimension xx(*), qq(*), pdf(m,*), fff(mxx0*mqq0)

      dimension ww(15+9*(mxx0+mqq0)+27*mxx0*mqq0)
      dimension yy(mxx0*mqq0)     , tt(mxx0*mqq0)
      dimension ipointy(mxx0*mqq0), ipointt(mxx0*mqq0)

      character*80 subnam
      data subnam /'FFTABL ( ISET, C, ISEL, X, NX, Q, NQ, F, M, ICHK )'/

C--   Inline address function for 2-dim array
      iaff(i,j,n) = i + n*(j-1)

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check pdf set
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Check status bits  (also check if pdf set is filled)
      call sqcChkflg(jset,ichk,subnam)
C--   Check pdf set is filled with current parameters
      call sqcParMsg(subnam,'ISET',jset)
C--   Check selection parameter
      if(isel.le.12) then
        call sqcIlele(subnam,'ISEL',0,isel,9,
     +                                    'Invalid selection parameter')
      else
        call sqcIlele(subnam,'ISEL',13,isel,ilast7(jset),
     +                'Attempt to access nonexisting extra pdf in ISET')
      endif
C--   Check number of entries in the table
      call sqcIlele(subnam,'M',nx,m,999999,'M must be >= NX')
      call sqcIlele(subnam,'NX+NQ',1,nx+nq,mxx0+mqq0,
     +                    'NX+NQ cannot exceed MXX0+MQQ0 in qcdnum.inc')
      call sqcIlele(subnam,'NX*NQ',1,nx*nq,mxx0*mqq0,
     +                    'NX*NQ cannot exceed MXX0*MQQ0 in qcdnum.inc')
C--   Check pdf type
      if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALXQ instead')

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Catch x = 1
      if(lmb_eq(xx(nx),1.D0,-aepsi6)) then
        mx = nx-1
        do iq = 1,nq
          pdf(nx,iq) = 0.D0
        enddo
      else
        mx = nx
      endif

C--   Findout x-range
      yma = ygrid2(iymac5)
      xmi = exp(-yma)
      xma = xmaxc2    !exclude x = 1
      call sqcRange(xx,mx,xmi,xma,aepsi6,ixmi,ixma,ierrx)
      if(ierrx.eq.2)
     +    call sqcErrMsg(subnam,'X not in strictly ascending order')
      if(jchk.ne.0 .and. (ixmi.ne.1 .or. ixma.ne.mx))
     +    call sqcErrMsg(subnam,'At least one X(i) out of range')
C--   Findout mu2-range
      tmi = tgrid2(itmic5)
      tma = tgrid2(itmac5)
      qmi = exp(tmi)
      qma = exp(tma)
      call sqcRange(qq,nq,qmi,qma,aepsi6,iqmi,iqma,ierrq)
      if(ierrq.eq.2)
     +    call sqcErrMsg(subnam,'Q not in strictly ascending order')
      if(jchk.ne.0 .and. (iqmi.ne.1 .or. iqma.ne.nq))
     +    call sqcErrMsg(subnam,'At least one Q(i) out of range')

C--   Preset output table
      do iq = 1,nq
        do ix = 1,mx
          pdf(ix,iq) = qnull6
        enddo
      enddo

C--   x or qmu2 completely out of range
      if(ierrx.ne.0 .or. ierrq.ne.0) return

      ny = 0
      do i = ixmi,ixma
        ny          =  ny+1
        yy(ny)      = -log(xx(i))
        ipointy(ny) =  i
      enddo
      nt = 0
      do i = iqmi,iqma
        nt          =  nt+1
        tt(nt)      =  log(qq(i))
        ipointt(nt) =  i
      enddo

C--   Setup the par array for sqcFillBuffer
      call sqcParForSumPdf(jset,c,isel,par,8+8*mpdf0,np,kerr)

C--   Initialise table
      call sqcTabIni(yy,ny,tt,nt,ww,15+9*(ny+nt)+27*ny*nt,nused,ierr)
      if(ierr.eq.1) stop 'FFTABL Init: not enough space in ww'
      if(ierr.eq.2) stop 'FFTABL Init: no scratch buffer available'

C--   Fill buffer with function values
      call sqcFillBuffer(dqcPdfSum,stor7,par,np,ww,ierr,jerr)
      if(ierr.eq.1) stop 'FFTABL Fill: ww not initialised'
      if(ierr.eq.2) stop 'FFTABL Fill: evolution parameter change'
      if(ierr.eq.3) stop 'FFTABL Fill: no scratch buffer available'
      if(ierr.eq.4) stop 'FFTABL Fill: error from dqcPdfSum'

      call sqcTabFun(ww,fff,ierr)
      if(ierr.eq.1) stop 'FFTABL TabF: ww not initialised'
      if(ierr.eq.2) stop 'FFTABL TabF: evolution parameter change'
      if(ierr.eq.3) stop 'FFTABL TabF: found no buffer to interpolate'

C--   Fill output table
      do iy = 1,ny
        ix = ipointy(iy)
        do it = 1,nt
          iq = ipointt(it)
          pdf(ix,iq) = fff(iaff(iy,it,ny))
        enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
CXXHDR  void ftable(int iset, double *c, int m, double *x, int nx, double *q, int nq, double *table, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fftable FC_FUNC(ftable,FTABLE)
CXXHFW   void fftable(int*, double*, int*, double*, int*, double*, int*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void ftable(int iset, double *c, int m, double *x, int nx, double *q, int nq, double *table, int ichk)
CXXWRP  {
CXXWRP    fftable(&iset,c,&m,x,&nx,q,&nq,table,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================================
      subroutine FTable(jset,c,isel,xx,nx,qq,nq,pdf,jchk)
C     ===================================================

C--   Pdf linear combination interpolated on a grid
C--   Routine kept for backward compatibility
C--
C--   jset       (in) :  pdf set [1,mset0]
C--   c(-6:6)    (in) :  coefficients of a linear combination
C--   isel       (in) :  number of extra pdfs beyond gluon and quarks
C--   xx         (in) :  input array of x-values
C--   nx         (in) :  number of x-values
C--   qq         (in) :  input array of mu2-values
C--   nq         (in) :  number of mu2-values
C--   pdf       (out) :  pdf(nx,nq) 2-dim array of interpolated pdfs
C--   jchk       (in) :  # 0: error if xx,qq out of range
C--                      = 0: return null if xx,qq out of range

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical lmb_eq

      external dqcPdfSum

      dimension c(-6:6), par(8+8*mpdf0)

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)
      
      dimension xx(*), qq(*), pdf(nx,*), fff(mxx0*mqq0)

      dimension ww(15+9*(mxx0+mqq0)+27*mxx0*mqq0)
      dimension yy(mxx0*mqq0)     , tt(mxx0*mqq0)
      dimension ipointy(mxx0*mqq0), ipointt(mxx0*mqq0)

      character*80 subnam
      data subnam /'FTABLE ( ISET, C, ISEL, X, NX, Q, NQ, F, ICHK )'/

C--   Inline address function for 2-dim array
      iaff(i,j,n) = i + n*(j-1)

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check pdf set
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Check status bits  (also check if pdf set is filled)
      call sqcChkflg(jset,ichk,subnam)
C--   Check pdf set is filled with current parameters
      call sqcParMsg(subnam,'ISET',jset)
C--   Check selection parameter
      if(isel.le.12) then
        call sqcIlele(subnam,'ISEL',0,isel,9,
     +                                    'Invalid selection parameter')
      else
        call sqcIlele(subnam,'ISEL',13,isel,ilast7(jset),
     +                'Attempt to access nonexisting extra pdf in ISET')
      endif
C--   Check number of entries in the table
      call sqcIlele(subnam,'NX+NQ',1,nx+nq,mxx0+mqq0,
     +                    'NX+NQ cannot exceed MXX0+MQQ0 in qcdnum.inc')
      call sqcIlele(subnam,'NX*NQ',1,nx*nq,mxx0*mqq0,
     +                    'NX*NQ cannot exceed MXX0*MQQ0 in qcdnum.inc')
C--   Check pdf type
      if(itypf7(jset).eq.5) call sqcErrMsg(subnam,
     + 'Cant handle user-defined pdf set (type-5): call BVALXQ instead')

C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

C--   Catch x = 1
      if(lmb_eq(xx(nx),1.D0,-aepsi6)) then
        mx = nx-1
        do iq = 1,nq
          pdf(nx,iq) = 0.D0
        enddo
      else
        mx = nx
      endif

C--   Findout x-range
      yma = ygrid2(iymac5)
      xmi = exp(-yma)
      xma = xmaxc2    !exclude x = 1
      call sqcRange(xx,mx,xmi,xma,aepsi6,ixmi,ixma,ierrx)
      if(ierrx.eq.2) 
     +    call sqcErrMsg(subnam,'X not in strictly ascending order')
      if(jchk.ne.0 .and. (ixmi.ne.1 .or. ixma.ne.mx))
     +    call sqcErrMsg(subnam,'At least one X(i) out of range')
C--   Findout mu2-range
      tmi = tgrid2(itmic5)
      tma = tgrid2(itmac5)
      qmi = exp(tmi)
      qma = exp(tma)
      call sqcRange(qq,nq,qmi,qma,aepsi6,iqmi,iqma,ierrq)
      if(ierrq.eq.2) 
     +    call sqcErrMsg(subnam,'Q not in strictly ascending order')
      if(jchk.ne.0 .and. (iqmi.ne.1 .or. iqma.ne.nq))
     +    call sqcErrMsg(subnam,'At least one Q(i) out of range')

C--   Preset output table
      do iq = 1,nq
        do ix = 1,mx
          pdf(ix,iq) = qnull6
        enddo            
      enddo

C--   x or qmu2 completely out of range
      if(ierrx.ne.0 .or. ierrq.ne.0) return

      ny = 0
      do i = ixmi,ixma
        ny          =  ny+1
        yy(ny)      = -log(xx(i))
        ipointy(ny) =  i
      enddo
      nt = 0
      do i = iqmi,iqma
        nt          =  nt+1
        tt(nt)      =  log(qq(i))
        ipointt(nt) =  i
      enddo

C--   Setup the par array for sqcFillBuffer
      call sqcParForSumPdf(jset,c,isel,par,8+8*mpdf0,np,kerr)

C--   Initialise table
      call sqcTabIni(yy,ny,tt,nt,ww,15+9*(ny+nt)+27*ny*nt,nused,ierr)
      if(ierr.eq.1) stop 'FTABLE Init: not enough space in ww'
      if(ierr.eq.2) stop 'FTABLE Init: no scratch buffer available'

C--   Fill buffer with function values
      call sqcFillBuffer(dqcPdfSum,stor7,par,np,ww,ierr,jerr)
      if(ierr.eq.1) stop 'FTABLE Fill: ww not initialised'
      if(ierr.eq.2) stop 'FTABLE Fill: evolution parameter change'
      if(ierr.eq.3) stop 'FTABLE Fill: no scratch buffer available'
      if(ierr.eq.4) stop 'FTABLE Fill: error from dqcPdfSum'

      call sqcTabFun(ww,fff,ierr)
      if(ierr.eq.1) stop 'FTABLE TabF: ww not initialised'
      if(ierr.eq.2) stop 'FTABLE TabF: evolution parameter change'
      if(ierr.eq.3) stop 'FTABLE TabF: found no buffer to interpolate'

C--   Fill output table
      do iy = 1,ny
        ix = ipointy(iy)
        do it = 1,nt
          iq = ipointt(it)
          pdf(ix,iq) = fff(iaff(iy,it,ny))
        enddo            
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ffplot(string fnam, double(* fun)(int*, double*, bool*), int m, double zmi, double zma, int n, string txt);
C-----------------------------------------------------------------------
CXXHFW  #define fffplotcpp FC_FUNC(ffplotcpp,FFPLOTCPP)
CXXHFW    void fffplotcpp(char*, int*, double (*)(int*, double*, bool*), int*, double*,  double*, int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void ffplot(string fname, double(* fun)(int*, double*, bool*), int m, double zmi, double zma, int n, string txt)
CXXWRP  {
CXXWRP    int ls1 = fname.size();
CXXWRP    char *cfname = new char[ls1+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    int ls2 = txt.size();
CXXWRP    char *ctxt = new char[ls2+1];
CXXWRP    strcpy(ctxt,txt.c_str());
CXXWRP    fffplotcpp(cfname,&ls1,fun,&m,&zmi,&zma,&n,ctxt,&ls2);
CXXWRP    delete[] cfname;
CXXWRP    delete[] ctxt;
CXXWRP  }
C-----------------------------------------------------------------------

C     ======================================================
      subroutine ffplotCPP(file,ls1,fun,m,zmi,zma,n,txt,ls2)
C     ======================================================

      implicit double precision (a-h,o-z)

      external fun
      character*(100) file, txt

      if(ls1.gt.100) stop 'ffplotCPP: input file name > 100 characters'
      if(ls2.gt.100) stop 'ffplotCPP: input txt string > 100 characters'

      call FFPlot(file(1:ls1),fun,m,zmi,zma,n,txt(1:ls2))

      return
      end

C     ===========================================
      subroutine FFPlot(file,fun,m,zmi,zma,n,txt)
C     ===========================================

C--   Write-out a plot file (to be used by gnuplot)
C--
C--   file     (in) : name of the output file
C--   fun(i,z) (in) : functions f_i(z) to plot
C--   m        (in) : number of functions to plot [1-50]
C--   zmi      (in) : lower limit of z
C--   zma      (in) : upper limit of z
C--   n        (in) : number of sample points [2-500]
C--                   n > 0  linear sampling
C--                   n < 0  logarithmic sampling

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'

      parameter(mmax = 50)                             !dont set it > 99
      parameter(nmax = 500)

      logical first
      dimension flist(mmax)
      external fun
      character*(*) file, txt
      character*9 format
      data format/'(  E13.5)'/

      character*80 subnam
      data subnam /'FFPLOT ( filename, fun, m, zmi, zma, n, txt )'/

      nz = abs(n)

C--   Check user input
      if(file.eq.' ') call sqcErrMsg(subnam,'Empty filename')
      call sqcIlele(subnam,'m',1,m ,mmax,' ')
      call sqcIlele(subnam,'n',2,nz,nmax,' ')
      if(zmi.ge.zma) call sqcErrMsg(subnam,
     +              'ZMI greater or equal than ZMA')
C--   Force positive limits when logscale
      if(n.lt.0 .and. (zmi.le.0.D0 .or. zma.le.0.D0))
     +  call sqcErrMsg(subnam,
     +                'Logarithmic sampling only when ZMA > ZMI > 0' )

C--   Open file
      lun = iqcLunFree(10)
      open(unit=lun,file=file,form='formatted',status='unknown',err=500)
C--   Write comment line if there is any
      if(imb_lenoc(txt).ne.0) write(lun,'(''#  '',A)') txt

C--   Set repeat count in the format descriptor
      write(format(2:3),'(I2)') m+1

C--   Go for it
      first  = .true.                           !flag first call to func
      rnull  =  qnull6                              !remember null value
      qnull6 =  0.D0                           !temporarily set null = 0
      if(n.gt.0) then
C--     linear sampling
        bw = (zma-zmi)/(nz-1)
        zz = zmi-bw
        do i = 1,nz
          zz       = zz+bw
          flist(1) = fun(1,zz,first)
          first    = .false.                  !unflag first call to func
          do j = 2,m
            flist(j) = fun(j,zz,first)
          enddo
          write(lun,fmt=format) zz,(flist(j),j=1,m)
        enddo
      else
C--     Logarithmic sampling
        zmil = log(zmi)
        zmal = log(zma)
        bw   = (zmal-zmil)/(nz-1)
        zlog = zmil-bw
        do i = 1,nz
          zlog     = zlog+bw
          zz       = exp(zlog)
          flist(1) = fun(1,zz,first)
          first    = .false.                  !unflag first call to func
          do j = 2,m
            flist(j) = fun(j,zz,first)
          enddo
          write(lun,fmt=format) zz,(flist(j),j=1,m)
        enddo
      endif

      write(lunerr1,'(/'' FFPLOT: write file '',A)') file

      qnull6 = rnull                                 !restore null value
      close(lun)

      return

C--   Error opening file
 500  continue
      call sqcErrmsg(subnam,'Cannot open file '//file)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void fiplot(string fnam, double(* fun)(int*, double*, bool*), int m, double* zval, int n, string txt);
C-----------------------------------------------------------------------
CXXHFW  #define ffiplotcpp FC_FUNC(fiplotcpp,FIPLOTCPP)
CXXHFW    void ffiplotcpp(char*, int*, double (*)(int*, double*, bool*), int*, double*, int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fiplot(string fname, double(* fun)(int*, double*, bool*), int m, double* zval, int n, string txt)
CXXWRP  {
CXXWRP    int ls1 = fname.size();
CXXWRP    char *cfname = new char[ls1+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    int ls2 = txt.size();
CXXWRP    char *ctxt = new char[ls2+1];
CXXWRP    strcpy(ctxt,txt.c_str());
CXXWRP    ffiplotcpp(cfname,&ls1,fun,&m,zval,&n,ctxt,&ls2);
CXXWRP    delete[] cfname;
CXXWRP    delete[] ctxt;
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================================
      subroutine fiplotCPP(file,ls1,fun,m,zval,n,txt,ls2)
C     ===================================================

      implicit double precision (a-h,o-z)

      external fun
      character*(100) file, txt
      dimension zval(*)

      if(ls1.gt.100) stop 'fiplotCPP: input file name > 100 characters'
      if(ls2.gt.100) stop 'fiplotCPP: input txt string > 100 characters'

      call FIPlot(file(1:ls1),fun,m,zval,n,txt(1:ls2))

      return
      end

C     ========================================
      subroutine FIPlot(file,fun,m,zval,n,txt)
C     ========================================

C--   Write-out a plot file (to be used by gnuplot)
C--
C--   file     (in) : name of the output file
C--   fun(i,z) (in) : functions f_i(z) to plot
C--   m        (in) : number of functions to plot [1-50]
C--   zval     (in) : input array of z-values
C--   n        (in) : number of z-values in zval [2-500]
C--                   n > 0  z = zval(i)
C--                   n < 0  z = dble(i) ignore zval array

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'

      parameter(mmax = 50)                             !dont set it > 99
      parameter(nmax = 500)

      logical first
      dimension flist(mmax), zval(*)
      external fun
      character*(*) file, txt
      character*9 format
      data format/'(  E13.5)'/

      character*80 subnam
      data subnam /'FIPLOT ( filename, fun, m, zval, n, txt )'/

      nz = abs(n)

C--   Check user input
      if(file.eq.' ') call sqcErrMsg(subnam,'Empty filename')
      call sqcIlele(subnam,'m',1,m ,mmax,' ')
      call sqcIlele(subnam,'n',2,nz,nmax,' ')
      if(n.gt.0) then
        do i = 2,nz
          if(zval(i-1).ge.zval(i)) call sqcErrMsg(subnam,
     +                'ZVAL not in strictly ascending order')
        enddo
      endif

C--   Open file
      lun = iqcLunFree(10)
      open(unit=lun,file=file,form='formatted',status='unknown',err=500)
C--   Write comment line if there is any
      if(imb_lenoc(txt).ne.0) write(lun,'(''#  '',A)') txt

C--   Set repeat count in the format descriptor
      write(format(2:3),'(I2)') m+1

C--   Go for it
      first  = .true.                           !flag first call to func
      rnull  =  qnull6                              !remember null value
      qnull6 =  0.D0                           !temporarily set null = 0
      if(n.gt.0) then
C--     Use zval array
        do i = 1,nz
          zz       = zval(i)
          flist(1) = fun(1,zz,first)
          first    = .false.                  !unflag first call to func
          do j = 2,m
            flist(j) = fun(j,zz,first)
          enddo
          write(lun,fmt=format) zz,(flist(j),j=1,m)
        enddo
      else
C--     Use z = dble(i)
        do i = 1,nz
          zz       = dble(i)
          flist(1) = fun(1,zz,first)
          first    = .false.                  !unflag first call to func
          do j = 2,m
            flist(j) = fun(j,zz,first)
          enddo
          write(lun,fmt=format) zz,(flist(j),j=1,m)
        enddo
      endif

      write(lunerr1,'(/'' FIPLOT: write file '',A)') file

      qnull6 = rnull                                 !restore null value
      close(lun)

      return

C--   Error opening file
 500  continue
      call sqcErrmsg(subnam,'Cannot open file '//file)

      return
      end

C==   ==================================================================
C==   Pdf table utilities  =============================================
C==   ==================================================================

C     ===============================================
      subroutine sqcPdfLims(idg,iy1,iy2,it1,it2,jmax)
C     ===============================================

C--   Return stor7 pdf table index limits

C--   idg     (in)  stor7 pdf identifier in global format
C--   iy1    (out)  lower y index
C--   iy2    (out)  upper y index
C--   it1    (out)  lower t index
C--   it2    (out)  upper t index
C--   jmax   (out)  upper j index in satellite table(j,id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension imin(6),imax(6)

      call sqcGetLimits(stor7,idg,imin,imax,jmax)
      iy1 = imin(1)            !iymin
      iy2 = imax(1)            !iymax
      it1 = imin(2)            !itmin
      it2 = imax(2)            !itmax

      return
      end

C     ==========================================================
      subroutine sqcParForSumPdf(jset,coef,isel,par,n,nout,ierr)
C     ==========================================================

C--   Setup parameter store for the function dqcSumPdf

C--   jset        (in) : pdf set addressed
C--   coef(-6:6)  (in) : quark coefficients
C--   isel        (in) : selection flag
C--   par        (out) : parameter store
C--   n           (in) : dim of par declared in the calling routine
C--   nout       (out) : number of words used in par
C--   ierr       (out) : 0 = all OK
C--                      1 = not enough space in par
C--
C--   Layout of par(8 + 13*8)
C--   -----------------------
C--   par(1)  =  Karr(0)
C--   par(2)  =  Karr(1)
C--   par(3)  =  Karr(2)
C--   par(4)  =  Karr(3)
C--   par(5)  =  # pdfs for nf = 3
C--   par(6)  =  # pdfs for nf = 4
C--   par(7)  =  # pdfs for nf = 5
C--   par(8)  =  # pdfs for nf = 6
C--
C--   par(9)  =  first word of Parr(i,j,k)
C--   ..          ..
C--
C--   Layout of Parr(2,13,3:6)
C--   --------------------------
C--   Parr(1,j,nf) = base address of pdf j at nf
C--   Parr(2,j,nf) = coefficient  of pdf j at nf

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension coef(-6:6), wt(13), id(13)
      dimension par(*)
      dimension imin(3), imax(3), karr(0:3)

C--   Address function in par
      iaP(i,j,k) = int(par(1))+int(par(2))*i+int(par(3))*j+int(par(4))*k

      ierr = 0

C--   Setup the par array for a maximum of 13 pdfs
      imin(1) = 1
      imax(1) = 2
      imin(2) = 1
      imax(2) = 13
      imin(3) = 3
      imax(3) = 6
      istart  = 9
      call smb_bkmat(imin,imax,karr,3,istart,nout)

C--   Check enough space
      if(nout.gt.n) then
        ierr = 1
        return
      endif

C--   Store partition parameters
      par(1) = dble(karr(0))
      par(2) = dble(karr(1))
      par(3) = dble(karr(2))
      par(4) = dble(karr(3))

      do nf = 3,6
        call sqcElistFF(coef,isel,wt,id,nid,nf)  !nid might be zero
        par(nf+2)   = dble(nid)
        do j = 1,nid
          ia        = iaP(1,j,nf)
          idg       = iqcIdPdfLtoG(jset,id(j))
          par(ia)   = dble(iqcG5ijk(stor7,1,1,idg))
          par(ia+1) = wt(j)
        enddo
      enddo

      return
      end

C     ==================================
      integer function iqcGimmeScratch()
C     ==================================

C--   Returns the global idg of an unfilled scratch buffer in stor7
C--   This buffer is then flagged as filled and contents set to zero
C--
C--   0 = no buffer available
C--
C--   To release the buffer, call sqcReleaseScratch(idg)

      implicit double precision (a-h,o-z)

      logical lqcIsFilled

      include 'qcdnum.inc'
      include 'qstor7.inc'

      iqcGimmeScratch = 0

      do id = ifrst7(0),ilast7(0)
        idg = iqcIdPdfLtoG(0,id)
        if(.not.lqcIsFilled(stor7,idg)) then
          iqcGimmeScratch = idg
          call sqcValidate(stor7,idg)
          call sqcPreset(idg,0.D0)
          return
        endif
      enddo

*      stop 'iqcGimmeScratch: no empty scratch buffer found --> STOP'

      return
      end

C     =================================
      subroutine sqcReleaseScratch(idg)
C     =================================

C--   Release scratch buffer idg (global id) in stor7 by flagging it
C--   as unfilled

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      call sqcInvalidate(stor7,idg)

      return
      end

