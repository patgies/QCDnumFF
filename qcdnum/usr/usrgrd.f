
C--   This is the file usrgrd.f containing the qcdnum grid routines

C--   subroutine gxmake(xmi,iwt,n,nxin,nxout,iosp)
C--   integer function ixfrmx(x)
C--   logical function xxatix(x,ix)
C--   double precision function xfrmix(ix)
C--   subroutine gxcopy(array,n,nx)
C--   subroutine gqmake(qq,ww,n,nqin,nqout)
C--   integer function iqfrmq(q)
C--   logical function qqatiq(q,iq)
C--   double precision function qfrmiq(iq)
C--   subroutine gqcopy(array,n,nq)
C--   subroutine grpars(nx,xmi,xma,nq,qmi,qma,iosp)
C--   subroutine setlim(ixmi,iqmi,iqma,roots)
C--   subroutine getlim(jset,xmi,q2mi,q2ma,roots)
C--   subroutine sqcFilLims(ixmi,itmi,itma)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM grid routines from usrgrd.f                      */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM grid routines from usrgrd.f                        */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM grid routines from usrgrd.f                        */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------
      
C==   ==================================================================
C==   x-Grid routines ==================================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void gxmake(double *xmin, int *iwt, int n, int nxin, int &nxout, int iord);
C-----------------------------------------------------------------------
CXXHFW  #define fgxmake FC_FUNC(gxmake,GXMAKE)
CXXHFW    void fgxmake(double*, int*, int*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void gxmake(double *xmin, int *iwt, int n, int nxin, int &nxout, int iord)
CXXWRP  {
CXXWRP    fgxmake(xmin,iwt,&n,&nxin,&nxout,&iord);
CXXWRP  }
C-----------------------------------------------------------------------

C     ============================================
      subroutine gxmake(xmi,iwx,n,nxin,nxout,iosp)
C     ============================================

C--   Define logarithmic x-grid
C--
C--   xmi(n)   = (in)  list of lowest gridpoints for each subgrid
C--   iwx(n)   = (in)  list of point density weights for each subgrid
C--   n        = (in)  number of subgrids
C--   nxin     = (in)  requested number of grid points
C--   nxout    = (out) generated number of grid points
C--   iosp     = (in)  spline interpolation 2=lin, 3=quad

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'pstor8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension xmi(*)   ,iwx(*)
      dimension yma(mxg0),iwy(mxg0)

      character*80 subnam
      data subnam /'GXMAKE ( XMI, IWT, NGR, NXIN, NXOUT, IORD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check if grid already defined
      if(Lygrid2) call sqcErrMsg2(subnam,
     + 'X-grid already defined',
     + 'To change grid, call QCINIT and start from scratch')

C--   Check user input
C--   1-check if spline order OK
      call sqcIlele(subnam,'IORD',2,iosp,3,
     + 'Only linear (2) or quadratic (3) interpolation is allowed')
C--   2-Check if number of subgrids OK
      call sqcIlele(subnam,'NGR',1,n,mxg0,
     + 'Remark: you can increase mxg0 in qcdnum.inc and recompile')
C--   3-Check if number of gridpoints OK 
      call sqcIlele(subnam,'NXIN', max(iosp,n), nxin, mxx0-11,
     + 'Remark: you can increase mxx0 in qcdnum.inc and recompile')
C--   4-Check if all xmi(i) are in range
      do i = 1,n
        call sqcDltlt(subnam,'XMI(i)',0.D0,xmi(i),1.D0,
     +  'At least one of the XMI(i) outside allowed range')
      enddo
C--   5-Check if all xmi(i) are in ascending order
      if(n.ge.2) then
        do i = 2,n
          if(xmi(i).le.xmi(i-1)) call sqcErrMsg(subnam,
     +    'XMI(i) not in ascending order')
        enddo
      endif
C--   6-Check that all weights are ascending integer multiples
      if(iwx(1).le.0) call sqcErrMsg(subnam, 
     +       'Zero or negative weight encountered')
      do i = 2,n
        if(iwx(i).le.0) call sqcErrMsg(subnam, 
     +     'Zero or negative weight encountered')
        irat = iwx(i)/iwx(i-1)
        if(iwx(i).ne.irat*iwx(i-1)) call sqcErrMsg(subnam, 
     +     'Weights are not ascending integer multiples')
      enddo

C--   Transform x-grid to y-grid
      do i = 1,n
        j = n+1-i
        yma(j) = -log(xmi(i))
        iwy(j) =  iwx(i)
      enddo
C--   Do the work
      call sqcGryDef(yma,iwy,n,nxin,nxout,iosp)
C--   Final check: more than 10 x-grid points required
      if(nxout.le.10) call sqcErrMsg(subnam,
     +   'More than 10 x-grid points required')

C--   Both grids are defined
      if(Ltgrid2) then
C--     Set default limits
        call sqcFilLims(1,1,ntt2)
C--     Initialise base store
        call sqcIniStore(nw,ierr)
        if(ierr.ne.0) call sqcMemMsg(subnam,nw,ierr)
C--     Setup parameter store
        call sparInit(nused)
        if(nused.le.0) call sqcMemMsg(subnam,nused,-4)
C--     Set grid version
        igver2 = igver2+1
C--     Set parameter bits
        ipbits8 = 0
        call smb_sbit1(ipbits8,infbit8)      !no nfmap
        call smb_sbit1(ipbits8,iasbit8)      !no alfas tables
        call smb_sbit1(ipbits8,izcbit8)      !no iz cuts
        call smb_sbit1(ipbits8,ipsbit8)      !no param store
C--     Construct base parameter set
        call sparMakeBase
      endif

C--   Invalidate weight tables
      Lwtini7 = .false.
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    int ixfrmx(double x);
C-----------------------------------------------------------------------
CXXHFW  #define fixfrmx FC_FUNC(ixfrmx,IXFRMX)
CXXHFW    int fixfrmx(double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int ixfrmx(double x)
CXXWRP  {
CXXWRP    return fixfrmx(&x);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================
      integer function ixfrmx(x)
C     ==========================

C--   Gives binnumber ix, given x
C--   ix = 0 if x < xmin, x > 1 or if xgrid not defined.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      logical lmb_eq
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'IXFRMX ( X )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      ixfrmx = 0
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                return
C--   Catch x = 1
      if(lmb_eq(x,1.D0,aepsi6)) then
        ixfrmx = nyy2(0)+1
      endif
C--   Check user input
      if(x.le.0.D0 .or. x.ge.1.D0) return
C--   Go...
      y         = -log(x)
      iy        = iqcFindIy(y)
      if(iqcYhitIy(y,iy).eq.1) then
        ixfrmx = nyy2(0) + 1 - iy
      else
        ixfrmx = nyy2(0) - iy
      endif
 
      return
      end

C-----------------------------------------------------------------------
CXXHDR    int xxatix(double x, int ix);
C-----------------------------------------------------------------------
CXXHFW  #define fxxatix FC_FUNC(xxatix,XXATIX)
CXXHFW    int fxxatix(double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int xxatix(double x, int ix)
CXXWRP  {
CXXWRP    return fxxatix(&x,&ix);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     =============================
      logical function xxatix(x,ix)
C     =============================

C--   True if x is at gridpoint ix.
C--   False of not at gridpoint, x or ix out of range or no xgrid defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      logical lmb_eq
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'XXATIX ( X, IX )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      xxatix = .false.
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                  return
C--   Catch x = 1
      if(lmb_eq(x,1.D0,aepsi6) .and. ix.eq.nyy2(0)+1) then
        xxatix = .true.
      endif
C--   Check user input
      ymax = ygrid2(nyy2(0))
      xmin = exp(-ymax)
      if(x.lt.xmin .or. x.ge.1.D0)   return
      if(ix.lt.1 .or. ix.gt.nyy2(0)) return
C--   Go...
      ihit  = iqcYhitIy(-log(x),nyy2(0)+1-ix)
      if(ihit.eq.1) then
        xxatix = .true.
      else
        xxatix = .false.     
      endif

      return
      end


C-----------------------------------------------------------------------
CXXHDR    double xfrmix(int ix);
C-----------------------------------------------------------------------
CXXHFW  #define fxfrmix FC_FUNC(xfrmix,XFRMIX)
CXXHFW    double fxfrmix(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double xfrmix(int ix)
CXXWRP  {
CXXWRP    return fxfrmix(&ix);
CXXWRP  }
C-----------------------------------------------------------------------

C     ====================================
      double precision function xfrmix(ix)
C     ====================================

C--   Get value of x-grid point ix.
C--   Returns 0.D0 if ix out of range [1,nyy2] or xgrid not defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'XFRMIX ( IX )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      xfrmix = 0.D0
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)               return
C--   Catch x = 1
      if(ix.eq.nyy2(0)+1) then
        xfrmix = 1.D0
        return
      endif
C--   Check user input
      if(ix.lt.1 .or. ix.gt.nyy2(0)) return
C--   Go...     
      iy     = nyy2(0) + 1 - ix
      yy     = ygrid2(iy)
      xfrmix = exp(-yy)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void gxcopy(double *array, int n, int &nx);
C-----------------------------------------------------------------------
CXXHFW  #define fgxcopy FC_FUNC(gxcopy,GXCOPY)
CXXHFW    void fgxcopy(double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void gxcopy(double *array, int n, int &nx)
CXXWRP  {
CXXWRP    fgxcopy(array,&n,&nx);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      subroutine gxcopy(array,n,nx)
C     =============================

C--   Copy x grid to local array
C--
C--   Input      n = dimension of array as declared in the calling routine
C--   Output array = target array declared in the calling routine
C--             nx = number of x grid point copied to array

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension array(*)

      character*80 subnam
      data subnam /'GXCOPY ( XARRAY, N, NX )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check user input
      call sqcIlele(subnam,'N',nyy2(0),n,100000,
     &              'XARRAY not large enough to contain x-grid')

      nx = nyy2(0)     
      do ix = 1,nx
        iy         = nyy2(0) + 1 - ix
        yy         = ygrid2(iy)
        array(ix)  = exp(-yy)
      enddo
      
      return
      end

C==   ==================================================================
C==   Q2-Grid routines =================================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void gqmake(double *qarr, double *wgt, int n, int nqin, int &nqout);
C-----------------------------------------------------------------------
CXXHFW  #define fgqmake FC_FUNC(gqmake,GQMAKE)
CXXHFW    void fgqmake(double*, double*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void gqmake(double *qarr, double *wgt, int n, int nqin, int &nqout)
CXXWRP  {
CXXWRP    fgqmake(qarr,wgt,&n,&nqin,&nqout);
CXXWRP  }
C-----------------------------------------------------------------------

C     =====================================
      subroutine gqmake(qq,ww,n,nqin,nqout)
C     =====================================

C--   Define logarithmic Q2-grid
C--
C--   qq(n)   (in)  List of Q2 values. qq(1) and qq(n) are the grid limits
C--   ww(n)   (in)  Weights: generated point density beween qq(i) and qq(i+1)
C--                 will be proportional to ww(i)   
C--   n       (in)  Number of points in qq and ww (>=2)
C--   nqin    (in)  Requested number of grid points. If <=n then qq will be
C--                 copied to the internal grid. If < 0 loglog spacing
C--   nqout   (out) Number of generated grid points

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'pstor8.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension qq(*),ww(*),qlog(mqq0)
      
      logical loglog

      character*80 subnam
      data subnam /'GQMAKE ( QARR, WGT, N, NQIN, NQOUT )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check if grid already exists
      if(Ltgrid2) call sqcErrMsg2(subnam,
     + 'Q2-grid already defined',
     + 'To change grid, call QCINIT and start from scratch')

C--   Check user input
      call sqcIlele(subnam,'N',2,n,mqq0-5,
     + 'Remark: You can increase mqq0 in qcdnum.inc and recompile')
      call sqcIlele(subnam,'NQIN',n,nqin,mqq0-10,
     + 'Remark: You can increase mqq0 in qcdnum.inc and recompile')
C--   Q2 should be larger than 0.1
      call sqcDltlt(subnam,'QARR(1)',qlimd6,qq(1),qlimu6,
     + 'Remark: the allowed range can be changed by a call to SETVAL')
C--   Spacing should be larger than 0.01 GeV2 --> force ascending order
      do i = 2,n
        if( (qq(i-1)+1.D-2) .ge. qq(i) ) then 
          call sqcErrMsg(subnam, 
     +   'QARR(i) not ascending or spaced by less than 0.01 GeV2')
        endif
        call sqcDltlt(subnam,'QARR(i)',qlimd6,qq(i),qlimu6,
     + 'Remark: these Q2 limits can be changed by a call to SETVAL')
      enddo
C--   Weights should be between 0.1 and 10
      do i = 1,n-1
        call sqcDlele(subnam,'WGT(i)',0.1D0,ww(i),10.D0,
     +        'Weights should be in a reasonable range [0.1,10]')
      enddo
C--   Safety margin for max number of points
      call sqcIlele(subnam,'NQIN',5-mqq0,nqin,mqq0-5,
     + 'Remark: You can increase mqq0 in qcdnum.inc and recompile')

C--   Do the work
      do i = 1,n
        qlog(i) = log(qq(i))
      enddo
C--   Log or loglog, thats the question
C--   Loglog is disabled in the check above since nqin < 0 is not allowed
C--   It is disabled because the loglog option is not fully checked if OK
      if(nqin.gt.0) then       
        nq2    =  nqin
        loglog = .false.
      else
        nq2    = -nqin  
        loglog = .true.
      endif  
      call sqcGrTdef(qlog,ww,n,nq2,loglog,jerr)
      if(jerr.ne.0) then 
         write(lunerr1,*) 'sqcGrTdef jerr = ',jerr,' ---> STOP'
         stop
      endif
      nqout = nq2
C--   Safety margin for max number of points
      call sqcIlele(subnam,'NQOUT',2,nqout,mqq0-5,
     & 'Remark: You can increase mqq0 in qcdnum.inc and recompile')
C--   Final check: more than 10 q-grid points required
      if(nqout.le.10) call sqcErrMsg(subnam,
     +   'More than 10 Q2-grid points required')

C--   Both grids are defined
      if(Lygrid2) then
C--     Set default limits
        call sqcFilLims(1,1,ntt2)
C--     Initialise base store
        call sqcIniStore(nw,ierr)
        if(ierr.ne.0) call sqcMemMsg(subnam,nw,ierr)
C--     Setup parameter store
        call sparInit(nused)
C--     Set grid version
        igver2 = igver2+1
C--     Set parameter bits
        ipbits8 = 0
        call smb_sbit1(ipbits8,infbit8)      !no nfmap
        call smb_sbit1(ipbits8,iasbit8)      !no alfas tables
        call smb_sbit1(ipbits8,izcbit8)      !no iz cuts
        call smb_sbit1(ipbits8,ipsbit8)      !no param store
C--     Construct base parameter set
        call sparMakeBase
      endif

C--   Invalidate weight tables
      Lwtini7 = .false.
C--   Update status bits
      call sqcSetflg(iset,idel,0)
    
      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iqfrmq(double q2);
C-----------------------------------------------------------------------
CXXHFW  #define fiqfrmq FC_FUNC(iqfrmq,IQFRMQ)
CXXHFW    int fiqfrmq(double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int iqfrmq(double q2)
CXXWRP  {
CXXWRP    return fiqfrmq(&q2);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ==========================
      integer function iqfrmq(q)
C     ==========================

C--   Get index iq of first gridpoint below Q2.
C--   iq = 0 id q < qmin, q > qmax or no Q2 grid defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      logical lmb_lt, lmb_gt

      character*80 subnam
      data subnam /'IQFRMQ ( Q2 )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      iqfrmq = 0
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                         return
C--   Check user input
      if(q.le.0.D0)                         return
      t = log(q)
      if(lmb_lt(t,tgrid2(1)   ,aepsi6))   return
      if(lmb_gt(t,tgrid2(ntt2),aepsi6))   return
C--   Get binnumber
      iqfrmq = iqcItfrmt(t)  

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int qqatiq(double q2, int iq);
C-----------------------------------------------------------------------
CXXHFW  #define fqqatiq FC_FUNC(qqatiq,QQATIQ)
CXXHFW    int fqqatiq(double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int qqatiq(double q2, int iq)
CXXWRP  {
CXXWRP    return fqqatiq(&q2,&iq);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     =============================
      logical function qqatiq(q,jq)
C     =============================

C--   True if q at grid point iq
C--   False if q not at iq, q or iq not in range, or if grid not defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      logical lmb_lt, lmb_gt

      character*80 subnam
      data subnam /'QQATIQ ( Q2, IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      qqatiq = .false.
      iq     = abs(jq)
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                         return
C--   Check user input
      if(q.le.0.D0)                         return
      t = log(q)
      if(lmb_lt(t,tgrid2(1)   ,aepsi6))     return
      if(lmb_gt(t,tgrid2(ntt2),aepsi6))     return
      if(iq.lt.1 .or. iq.gt.ntt2)           return
C--   Go...
      ihit = iqcThitit(t,iq)
      if(ihit.eq.1) then
        qqatiq = .true.
      else
        qqatiq = .false.     
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double qfrmiq(int iq);
C-----------------------------------------------------------------------
CXXHFW  #define fqfrmiq FC_FUNC(qfrmiq,QFRMIQ)
CXXHFW    double fqfrmiq(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double qfrmiq(int iq)
CXXWRP  {
CXXWRP    return fqfrmiq(&iq);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ====================================
      double precision function qfrmiq(jq)
C     ====================================

C--   Get value of Q2-grid point iq.
C--   Returns 0.D0 if iq not in range [1,ntt2] or if Q2 grid not defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'QFRMIQ ( IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      qfrmiq = 0.D0
      iq     = abs(jq)
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                  return
C--   Check user input
      if(iq.lt.1 .or. iq.gt.ntt2)    return
C--   Go...     
      qfrmiq = exp(tgrid2(iq))
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void gqcopy(double *array, int n, int &nq);
C-----------------------------------------------------------------------
CXXHFW  #define fgqcopy FC_FUNC(gqcopy,GQCOPY)
CXXHFW    void fgqcopy(double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void gqcopy(double *array, int n, int &nq)
CXXWRP  {
CXXWRP    fgqcopy(array,&n,&nq);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     =============================
      subroutine gqcopy(array,n,nq)
C     =============================

C--   Copy q2 grid to local array
C--
C--   Input      n = dimension of array as declared in the calling routine
C--   Output array = target array declared in the calling routine
C--             nq = number of q2 grid point copied to array

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension array(*)

      character*80 subnam
      data subnam /'GQCOPY ( QARRAY, N, NQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check user input
      call sqcIlele(subnam,'N',ntt2,n,100000,
     &              'QARRAY not large enough to contain Q2-grid')

      nq = ntt2
      do iq = 1,nq     
        array(iq)= exp(tgrid2(iq))
      enddo
      
      return
      stop

      end

C==   ==================================================================
C==   General access to the x-Q2 grid ==================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void grpars(int &nx, double &xmi, double &xma, int &nq, double &qmi, double &qma, int &iord);
C-----------------------------------------------------------------------
CXXHFW  #define fgrpars FC_FUNC(grpars,GRPARS)
CXXHFW    void fgrpars(int*, double*, double*, int*, double*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void grpars(int &nx, double &xmi, double &xma, int &nq, double &qmi, double &qma, int &iord)
CXXWRP  {
CXXWRP    fgrpars(&nx,&xmi,&xma,&nq,&qmi,&qma,&iord);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================================
      subroutine grpars(nx,xmi,xma,nq,qmi,qma,iosp)
C     =============================================

C--   Returns the number of points and the limits of the x-Q2 grid.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'GRPARS ( NX, XMI, XMA, NQ, QMI, QMA, IORD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

      nx   = nyy2(0)
      xmi  = exp(-ymax2(0))
      xma  = 1.D0
      nq   = ntt2
      qmi  = exp(tgrid2(1))
      qma  = exp(tgrid2(ntt2))
      iosp = ioy2
      
      return
      end
      
C==   ==================================================================
C==   Evolution cuts ===================================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void setlim( int ixmin, int iqmin, int iqmax, double dummy);
C-----------------------------------------------------------------------
CXXHFW  #define fsetlim FC_FUNC(setlim,SETLIM)
CXXHFW    void fsetlim(int*, int*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP//---------------------------------------------------------------------
CXXWRP  void setlim( int ixmin, int iqmin, int iqmax, double dummy)
CXXWRP  {
CXXWRP    fsetlim(&ixmin,&iqmin,&iqmax,&dummy);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================================
      subroutine setlim(ixmi,iqmi,iqma,roots)
C     =======================================

C--   Set kinematic cuts

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'pstor8.inc'
      
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'SETLIM ( IXMIN, IQMIN, IQMAX, DUM )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check input
      call sqcIlele(subnam,'IXMIN',0,ixmi,nyy2(0),' ')
      call sqcIlele(subnam,'IQMIN',0,iqmi,ntt2,' ')
      call sqcIlele(subnam,'IQMAX',0,iqma,ntt2,' ')
      
C--   Do the work
      dummy = roots  !avoid compiler warning
C--   Bracket the cuts and force minimum range
      if(ixmi.ge.1 .and. ixmi.le.nyy2(0)) then
        jxmi = ixmi
      else
        jxmi = 1
      endif
      if(nyy2(0)-jxmi.lt.10) call sqcErrMsg(subnam,
     +   'More than 10 x-grid points required after cuts')
      if(iqmi.ge.1 .and. iqmi.le.ntt2) then
        jqmi = iqmi
      else
        jqmi = 1
      endif
      if(iqma.ge.1 .and. iqma.lt.ntt2) then
        jqma = iqma+1                !enlarge q2 range by one grid point
      else
        jqma = ntt2
      endif
      if(jqma-jqmi.lt.10) call sqcErrMsg(subnam,
     +   'More than 10 Q2-grid points required after cuts')
C--   Now properly fill /qgrid2/
      call sqcFilLims(jxmi,jqmi,jqma)

C--   Invalidate parameter store
      call smb_sbit1(ipbits8,ipsbit8)
C--   Invalidate z-cuts
      call smb_sbit1(ipbits8,izcbit8)
C--   Construct base parameter set
      call sparMakeBase

C--   Update status bits
      call sqcSetflg(iset,idel,0)      

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void getlim( int iset, double &xmin, double &qmin, double &qmax, double &dummy);
C-----------------------------------------------------------------------
CXXHFW  #define fgetlim FC_FUNC(getlim,GETLIM)
CXXHFW    void fgetlim(int*, double*, double*, double*, double*);
C-----------------------------------------------------------------------
CXXWRP//---------------------------------------------------------------------
CXXWRP  void getlim( int iset, double &xmin, double &qmin, double &qmax, double &dummy)
CXXWRP  {
CXXWRP    fgetlim(&iset,&xmin,&qmin,&qmax,&dummy);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===========================================
      subroutine getlim(jset,xmi,q2mi,q2ma,roots)
C     ===========================================

C--   Get kinematic cuts

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
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

      character*80 subnam
      data subnam /'GETLIM ( ISET, XMIN, QMIN, QMAX, DUM )'/
      
C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check that jset exists and is filled
      call sqcIlele(subnam,'ISET',0,jset,mset0,'ISET does not exist')
      if(.not.Lfill7(jset)) call sqcSetMsg(subnam,'ISET',jset)

C--   Do the work
      key   = int(dparGetPar(stor7,isetf7(jset),idipver8))
      iymac = int(dparGetPar(pars8,key,idiymac8))
      itmic = int(dparGetPar(pars8,key,iditmic8))
      itmac = int(dparGetPar(pars8,key,iditmac8))
      xmi   = exp(-ygrid2(iymac))
      q2mi  = exp(tgrid2(itmic))
      q2ma  = exp(tgrid2(itmac))
      roots = 0.D0

      return
      end

C     =====================================
      subroutine sqcFilLims(ixmi,itmi,itma)
C     =====================================

C--   Fill the common blocks in qgrid2 (after the grid is defined)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

C--   Fill the common block
      ixmic2 = ixmi
      iymac2 = nyy2(0)-ixmi+1
      itmic2 = itmi
      itmac2 = itma
      xminc2 = exp(-ygrid2(iymac2))
      xmaxc2 = 1.D0 - 2.D0*aepsi6
      ymaxc2 = ygrid2(iymac2)
      qminc2 = exp(tgrid2(itmic2))
      qmaxc2 = exp(tgrid2(itmac2))
      tminc2 = tgrid2(itmic2)
      tmaxc2 = tgrid2(itmac2)

      return
      end


