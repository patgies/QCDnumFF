
C--   This is the file usrsplint.f with user routines and C++ interfaces

C--   isp_SpVers()
C--   ssp_SpInit(nuser)
C--   ssp_Uwrite(i, val)
C--   dsp_Uread(i)

C--   isp_SxMake(istep)
C--   isp_SxUser(xarr, nx)
C--   ssp_SxFill(ias, fun, iq)
C--   ssp_SxFpdf(ias, iset, def, isel, iq)
C--   ssp_SxF123(ias, iset, def, istf, iq)
C--   isp_SqMake(istep)
C--   isp_SqUser(qarr, nq)
C--   ssp_SqFill(ias, fun, ix)
C--   ssp_SqFpdf(ias, iset, def, isel, ix)
C--   ssp_SqF123(ias, iset, def, istf, ix)
C--   dsp_FunS1(ia, z, ichk)

C--   dsp_IntS1(ia, u1, u2)

C--   isp_S2Make(istepx, istepq)
C--   isp_S2User(xarr, nx, qarr, nq)
C--   ssp_S2Fill(ias, fun, rs)
C--   ssp_S2Fpdf(ias, iset, def, isel, rs)
C--   ssp_S2F123(ias, iset, def, istf, rs)
C--   dsp_FunS2(ia, x, q, ichk)

C--   dsp_IntS2(ia, x1, x2, q1, q2, rs, np)

C--   ssp_ExtrapU(ia, n)
C--   ssp_ExtrapV(ia, n)

C--   isp_SplineType(ia)
C--   ssp_SpLims(ia, nu, umi, uma, nv, vmi, vma)
C--   ssp_Unodes(ia, array, n, nus)
C--   ssp_Vnodes(ia, array, n, nvs)
C--   ssp_Nprint(ia)
C--   dsp_RsCut(ia)
C--   dsp_RsMax(ia,rs)

C--   ssp_Mprint()
C--   isp_SpSize(ia)
C--   ssp_Erase(ia)
C--   isp_SpDump(ia, fname)
C--   isp_SpRead(fname)
C--   ssp_SpSetVal(ia, i, val)
C--   dsp_SpGetVal(ia, i)

C=======================================================================
C===  Initialisation  ==================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  SPLINT routines from usrsplint.f                        */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  SPLINT routines from usrsplint.f                          */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  SPLINT routines from usrsplint.f                          */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    int isp_spvers();
C-----------------------------------------------------------------------
CXXHFW  #define fisp_spvers FC_FUNC(isp_spvers,ISP_SPVERS)
CXXHFW    int fisp_spvers();
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_spvers()
CXXWRP    {
CXXWRP      return fisp_spvers();
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================
      integer function isp_SpVers()
C     =============================

      implicit double precision(a-h,o-z)

      include 'splint.inc'

      isp_SpVers = IvDate0

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_spinit(int nuser);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_spinit FC_FUNC(ssp_spinit,SSP_SPINIT)
CXXHFW    void fssp_spinit(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_spinit(int nuser)
CXXWRP    {
CXXWRP      fssp_spinit(&nuser);
CXXWRP    }
C-----------------------------------------------------------------------

C     ============================
      subroutine ssp_SpInit(nuser)
C     ============================

C--   SPLINT workspace initialisation
C--   nuser  (in): user store size (<= 0 --> no store)

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.1) stop
     +   ' SPLINT::SSP_SPINIT: splint memory already initialised'

C--   Go..
      write(6,'(/''  +---------------------------------------+'')')
      write(6,'( ''  | You are using SPLINT version '',I8,'' |'')')
     + ivdate0
      write(6,'( ''  +---------------------------------------+'')')
      write(6,'(/)')

      ia    = iws_WsInit(w, nw0, ntags0,
     +    'Increase NW0 in splint/inc/splint.inc and recompile SPLINT')
      iroot = iws_IaRoot()                        !root address
      iatag = iws_IaFirstTag(w, iroot)            !begin ws tagfield
      w(iatag+Iversion0) = dble(ivdate0)          !store release date
      if(nuser.gt.0) then
        iatab = iws_Wtable (w, 1, nuser, 1)       !book user space
        iausr = iws_BeginTbody(w, iatab)          !begin user space
        w(iatag+IaUbegin0) = dble(iausr)          !store begin uspace
        w(iatag+NwUstore0) = dble(nuser)          !store size uspace
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_uwrite(int i, double val);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_uwrite FC_FUNC(ssp_uwrite,SSP_UWRITE)
CXXHFW    void fssp_uwrite(int*, double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_uwrite(int i, double val)
CXXWRP    {
CXXWRP      fssp_uwrite(&i, &val);
CXXWRP    }
C-----------------------------------------------------------------------

C     ============================
      subroutine ssp_Uwrite(i,val)
C     ============================

C--   Write user store

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.0) stop
     +   ' SPLINT::SSP_UWRITE: splint memory not initialised'
C--   Check user space
      iar = iws_IaRoot()                    !root address
      iat = iws_IaFirstTag(w, iar)          !begin ws tagfield
      iau = int(w(iat+IaUbegin0))           !begin uspace
      nus = int(w(iat+NwUstore0))           !size uspace
      if(nus.eq.0) stop
     +   ' SPLINT::SSP_UWRITE: no user space available'
      if(i.lt.1 .or. i.gt.nus) stop
     +   ' SPLINT::SSP_UWRITE: index I out of range'

      w(iau+i-1) = val

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dsp_uread(int i);
C-----------------------------------------------------------------------
CXXHFW  #define fdsp_uread FC_FUNC(dsp_uread,DSP_UREAD)
CXXHFW    double fdsp_uread(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dsp_uread(int i)
CXXWRP    {
CXXWRP      return fdsp_uread(&i);
CXXWRP    }
C-----------------------------------------------------------------------

C     ======================================
      double precision function dsp_Uread(i)
C     ======================================

C--   Read user store

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.0) stop
     +   ' SPLINT::DSP_UREAD: splint memory not initialised'
C--   Check user space
      iar = iws_IaRoot()                    !root address
      iat = iws_IaFirstTag(w, iar)          !begin ws tagfield
      iau = int(w(iat+IaUbegin0))           !begin uspace
      nus = int(w(iat+NwUstore0))           !size uspace
      if(nus.eq.0) stop
     +   ' SPLINT::DSP_UREAD: no user space available'
      if(i.lt.1 .or. i.gt.nus) stop
     +   ' SPLINT::DSP_UREAD: index I out of range'

      dsp_Uread = w(iau+i-1)

      return
      end

C=======================================================================
C===  1-dim splines  ===================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    int isp_sxmake(int istep);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_sxmake FC_FUNC(isp_sxmake,ISP_SXMAKE)
CXXHFW    int fisp_sxmake(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_sxmake(int istep)
CXXWRP    {
CXXWRP      return fisp_sxmake(&istep);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function isp_SxMake(istep)
C     ==================================

C--   istep   (in) : sample step in QCDNUM x-grid

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension ynodes(maxn0)

C--   Check input
      if(istep.le.0) stop ' SPLINT::ISP_SXMAKE: istep <= 0'
C--   Check QCDNUM initialised
      call GETINT('vers',ivers)
      if(ivers.eq.0) stop ' SPLINT::ISP_SXMAKE: QCDNUM not initialised'
C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.0) stop
     +   ' SPLINT::ISP_SXMAKE: splint memory not initialised'

C--   Automake y-node points
      call sspYnMake(istep, ynodes, nys, ierr)
      if(ierr.ne.0) goto 500

      isp_SxMake = ispS1make(w, ynodes, nys, -1)

      return

 500  continue
      write(6,*) ' '
      write(6,*) ' SPLINT::ISP_SXMAKE: Too many node points'
      write(6,*) '                     MAXN0 in splint.inc exceeded'
      stop

      end

C-----------------------------------------------------------------------
CXXHDR    int isp_sxuser(double* xarr, int nx);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_sxuser FC_FUNC(isp_sxuser,ISP_SXUSER)
CXXHFW    int fisp_sxuser(double *, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_sxuser(double* xarr, int nx)
CXXWRP    {
CXXWRP      return fisp_sxuser(xarr, &nx);
CXXWRP    }
C-----------------------------------------------------------------------

C     =====================================
      integer function isp_SxUser(xarr, nx)
C     =====================================

C--   xarr    (in) : array of proposed x-nodes
C--   nx      (in) : number of points in xarr

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension xarr(*)
      dimension ynodes(maxn0)

C--   Check input
      if(nx.lt.2) stop ' SPLINT::ISP_SXUSER: nx < 2'
C--   Check QCDNUM initialised
      call GETINT('vers',ivers)
      if(ivers.eq.0) stop ' SPLINT::ISP_SXUSER: QCDNUM not initialised'
C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.0) stop
     +   ' SPLINT::ISP_SXUSER: splint memory not initialised'

C--   Setup y-nodes
      call sspYnUser(xarr, nx, ynodes, nys, ierr)
      if(ierr.ne.0) goto 500

      isp_SxUser = ispS1make(w, ynodes, nys, -1)

      return

 500  continue
      write(6,*) ' '
      write(6,*) ' SPLINT::ISP_SXUSER: Too many node points'
      write(6,*) '                     MAXN0 in splint.inc exceeded'
      stop

      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_sxfill(int ias, double (*fun)(int*,int*,bool*),
CXXHDR                    int iq);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_sxfill FC_FUNC(ssp_sxfill,SSP_SXFILL)
CXXHFW    void fssp_sxfill(int*, double (*)(int*,int*,bool*), int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_sxfill(int ias, double (*fun)(int*,int*,bool*),
CXXWRP                    int iq)
CXXWRP    {
CXXWRP      fssp_sxfill(&ias, fun, &iq);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      subroutine ssp_SxFill(ias, fun, iq)
C     ===================================

C--   ias     (in) : address x-spline object
C--   fun     (in) : user defined fun(ix,iq,first)
C--   iq      (in) : fixed iq during fill

      implicit double precision(a-h,o-z)
      logical first

      include 'splint.inc'
      include 'spliws.inc'

      external fun

      dimension fvals(maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_SXFILL: input address IA out of range'
      if(ispSplineType(w,ias).ne.-1) stop
     +    ' SPLINT::SSP_SXFILL: input address IA is not an x-spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_SXFILL: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaOneD(w,ias,iat,iau,nus,iaf,iab,iac,iad)

C--   Initialise
      call smb_vfill(w(iaf),nus,0.D0)
      call smb_vfill(w(iab),nus,0.D0)
      call smb_vfill(w(iac),nus,0.D0)
      call smb_vfill(w(iad),nus,0.D0)

C--   Check iq
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)
      if(iq.le.0 .or. iq.gt.nq)    stop
     +    ' SPLINT::SSP_SXFILL: input iq out of range'

C--   Fill fvals
      first = .true.
      ia = iau-1
      do iy = 1,nus
        ia        = ia+1
        x         = exp(-w(ia))
        ix        = ixfrmx(x)
        fvals(iy) = fun(ix,iq,first)
        first     = .false.
      enddo

      call sspS1fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_sxfpdf(int ias, int iset, double *def, int isel,
CXXHDR                    int iq);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_sxfpdf FC_FUNC(ssp_sxfpdf,SSP_SXFPDF)
CXXHFW    void fssp_sxfpdf(int*, int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_sxfpdf(int ias, int iset, double *def, int isel,
CXXWRP                    int iq)
CXXWRP    {
CXXWRP      fssp_sxfpdf(&ias, &iset, def, &isel, &iq);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      subroutine ssp_SxFpdf(ias, iset, def, isel, iq)
C     ===============================================

C--   Fill x-spline with quark linear combination

C--   ias       (in) : address x-spline object
C--   iset      (in) : pdf set identifier
C--   def(-6:6) (in) : (anti-)quark coefficients
C--   isel      (in) : selection flag
C--   iq        (in) : fixed iq during fill

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension def(-6:6), xx(maxn0), qq(maxn0), fvals(maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_SXFPDF: input address IA out of range'
      if(ispSplineType(w,ias).ne.-1) stop
     +    ' SPLINT::SSP_SXFPDF: input address IA is not an x-spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_SXFPDF: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaOneD(w,ias,iat,iau,nus,iaf,iab,iac,iad)

C--   Initialise
      call smb_vfill(w(iaf),nus,0.D0)
      call smb_vfill(w(iab),nus,0.D0)
      call smb_vfill(w(iac),nus,0.D0)
      call smb_vfill(w(iad),nus,0.D0)

C--   Check iq
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)
      if(iq.le.0 .or. iq.gt.nq)    stop
     +    ' SPLINT::SSP_SXFPDF: input iq out of range'

C--   Make list x points
      q2 = qfrmiq(iq)
      ia = iau-1
      do iy = 1,nus
        ia        = ia+1
        xx(iy)    = exp(-w(ia))
        qq(iy)    = q2
      enddo

C--   Fill fvals
      call FFLIST(iset,def,isel,xx,qq,fvals,nus,1)

      call sspS1fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_sxf123(int ias, int iset, double *def, int istf,
CXXHDR                    int iq);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_sxf123 FC_FUNC(ssp_sxf123,SSP_SXF123)
CXXHFW    void fssp_sxf123(int*, int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_sxf123(int ias, int iset, double *def, int istf,
CXXWRP                    int iq)
CXXWRP    {
CXXWRP      fssp_sxf123(&ias, &iset, def, &istf, &iq);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      subroutine ssp_SxF123(ias, iset, def, istf, iq)
C     ===============================================

C--   Fill x-spline with structure function

C--   ias       (in) : address x-spline object
C--   iset      (in) : pdf set identifier
C--   def(-6:6) (in) : (anti-)quark coefficients
C--   istf      (in) : 1=FL, 2=F2, 3=xF3
C--   iq        (in) : fixed iq during fill

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension def(-6:6), xx(maxn0), qq(maxn0), fvals(maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_SXF123: input address IA out of range'
      if(ispSplineType(w,ias).ne.-1) stop
     +    ' SPLINT::SSP_SXF123: input address IA is not an x-spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_SXF123: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaOneD(w,ias,iat,iau,nus,iaf,iab,iac,iad)

C--   Initialise
      call smb_vfill(w(iaf),nus,0.D0)
      call smb_vfill(w(iab),nus,0.D0)
      call smb_vfill(w(iac),nus,0.D0)
      call smb_vfill(w(iad),nus,0.D0)

C--   Check iq
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)
      if(iq.le.0 .or. iq.gt.nq)    stop
     +    ' SPLINT::SSP_SXF123: input iq out of range'

C--   Make list of x and q2 points
      q2 = qfrmiq(iq)
      ia = iau-1
      do iy = 1,nus
        ia        = ia+1
        xx(iy)    = exp(-w(ia))
        qq(iy)    = q2
      enddo

C--   Fill fvals
      call ZSWITCH(iset)
      call ZMSTFUN(istf, def, xx, qq, fvals, nus, 1)

      call sspS1fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int isp_sqmake(int istep);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_sqmake FC_FUNC(isp_sqmake,ISP_SQMAKE)
CXXHFW    int fisp_sqmake(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_sqmake(int istep)
CXXWRP    {
CXXWRP      return fisp_sqmake(&istep);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function isp_SqMake(istep)
C     ==================================

C--   istep   (in) : sample step in QCDNUM q-grid

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension tnodes(maxn0)

C--   Check input
      if(istep.le.0) stop ' SPLINT::ISP_SQMAKE: istep <= 0'
C--   Check QCDNUM initialised
      call GETINT('vers',ivers)
      if(ivers.eq.0) stop ' SPLINT::ISP_SQMAKE: QCDNUM not initialised'
C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.0) stop
     +   ' SPLINT::ISP_SQMAKE: splint memory not initialised'

C--   Automake t-node points
      call sspTnMake(istep, tnodes, nts, ierr)
      if(ierr.ne.0) goto 500

      isp_SqMake = ispS1make(w, tnodes, nts, +1)

      return

 500  continue
      write(6,*) ' '
      write(6,*) ' SPLINT::ISP_SQMAKE: Too many node points'
      write(6,*) '                     MAXN0 in splint.inc exceeded'
      stop

      end

C-----------------------------------------------------------------------
CXXHDR    int isp_squser(double* qarr, int nq);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_squser FC_FUNC(isp_squser,ISP_SQUSER)
CXXHFW    int fisp_squser(double *, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_squser(double* qarr, int nq)
CXXWRP    {
CXXWRP      return fisp_squser(qarr, &nq);
CXXWRP    }
C-----------------------------------------------------------------------

C     =====================================
      integer function isp_SqUser(qarr, nq)
C     =====================================

C--   qarr    (in) : array of proposed q-nodes
C--   nq      (in) : number of qpoints in qarr

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension qarr(*)
      dimension tnodes(maxn0)

C--   Check input
      if(nq.lt.2) stop ' SPLINT::ISP_SQUSER: nq < 2'
C--   Check QCDNUM initialised
      call GETINT('vers',ivers)
      if(ivers.eq.0) stop ' SPLINT::ISP_SQUSER: QCDNUM not initialised'
C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.0) stop
     +   ' SPLINT::ISP_SQUSER: splint memory not initialised'

C--   Setup t-nodes
      call sspTnUser(qarr, nq, tnodes, nts, ierr)
      if(ierr.ne.0) goto 500

      isp_SqUser = ispS1make(w, tnodes, nts, +1)

      return

 500  continue
      write(6,*) ' '
      write(6,*) ' SPLINT::ISP_SQUSER: Too many node points'
      write(6,*) '                     MAXN0 in splint.inc exceeded'
      stop

      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_sqfill(int ias, double (*fun)(int*,int*,bool*),
CXXHDR                    int ix);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_sqfill FC_FUNC(ssp_sqfill,SSP_SQFILL)
CXXHFW    void fssp_sqfill(int*, double (*)(int*,int*,bool*), int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_sqfill(int ias, double (*fun)(int*,int*,bool*),
CXXWRP                    int ix)
CXXWRP    {
CXXWRP      fssp_sqfill(&ias, fun, &ix);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      subroutine ssp_SqFill(ias, fun, ix)
C     ===================================

C--   ias     (in) : address q-spline object
C--   fun     (in) : user defined fun(ix,iq,first)
C--   ix      (in) : fixed ix during fill

      implicit double precision(a-h,o-z)
      logical first

      include 'splint.inc'
      include 'spliws.inc'

      external fun

      dimension fvals(maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_SQFILL: input address IA out of range'
      if(ispSplineType(w,ias).ne.1)  stop
     +    ' SPLINT::SSP_SQFILL: input address IA is not an q-spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_SQFILL: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaOneD(w,ias,iat,iau,nus,iaf,iab,iac,iad)

C--   Initialise
      call smb_vfill(w(iaf),nus,0.D0)
      call smb_vfill(w(iab),nus,0.D0)
      call smb_vfill(w(iac),nus,0.D0)
      call smb_vfill(w(iad),nus,0.D0)

C--   Check iq
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)
      if(ix.le.0 .or. ix.gt.nx+1)    stop
     +    ' SPLINT::SSP_SQFILL: input ix out of range'

C--   Fill fvals
      first = .true.
      ia = iau-1
      do iy = 1,nus
        ia        = ia+1
        q         = exp(w(ia))
        iq        = iqfrmq(q)
        fvals(iy) = fun(ix,iq,first)
        first     = .false.
      enddo

      call sspS1fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_sqfpdf(int ias, int iset, double *def, int isel,
CXXHDR                    int ix);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_sqfpdf FC_FUNC(ssp_sqfpdf,SSP_SQFPDF)
CXXHFW    void fssp_sqfpdf(int*, int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_sqfpdf(int ias, int iset, double *def, int isel,
CXXWRP                    int ix)
CXXWRP    {
CXXWRP      fssp_sqfpdf(&ias, &iset, def, &isel, &ix);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      subroutine ssp_SqFpdf(ias, iset, def, isel, ix)
C     ===============================================

C--   Fill q-spline with quark linear combination

C--   ias       (in) : address x-spline object
C--   iset      (in) : pdf set identifier
C--   def(-6:6) (in) : (anti-)quark coefficients
C--   isel      (in) : selection flag
C--   ix        (in) : fixed ix during fill

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension def(-6:6), xx(maxn0), qq(maxn0), fvals(maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_SQFPDF: input address IA out of range'
      if(ispSplineType(w,ias).ne.-1) stop
     +    ' SPLINT::SSP_SQFPDF: input address IA is not an x-spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_SQFPDF: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaOneD(w,ias,iat,iau,nus,iaf,iab,iac,iad)

C--   Initialise
      call smb_vfill(w(iaf),nus,0.D0)
      call smb_vfill(w(iab),nus,0.D0)
      call smb_vfill(w(iac),nus,0.D0)
      call smb_vfill(w(iad),nus,0.D0)

C--   Check ix
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)
      if(ix.le.0 .or. ix.gt.nx+1)    stop
     +    ' SPLINT::SSP_SQFPDF: input ix out of range'

C--   Make list q2 points
      x  = xfrmix(ix)
      ia = iau-1
      do iy = 1,nus
        ia        = ia+1
        xx(iy)    = x
        qq(iy)    = exp(w(ia))
      enddo

C--   Fill fvals
      call FFLIST(iset,def,isel,xx,qq,fvals,nus,1)

      call sspS1fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_sqf123(int ias, int iset, double *def, int istf,
CXXHDR                    int ix);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_sqf123 FC_FUNC(ssp_sqf123,SSP_SQF123)
CXXHFW    void fssp_sqf123(int*, int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_sqf123(int ias, int iset, double *def, int istf,
CXXWRP                    int ix)
CXXWRP    {
CXXWRP      fssp_sqf123(&ias, &iset, def, &istf, &ix);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      subroutine ssp_SqF123(ias, iset, def, istf, ix)
C     ===============================================

C--   Fill q-spline with structure function

C--   ias       (in) : address q-spline object
C--   iset      (in) : pdf set identifier
C--   def(-6:6) (in) : (anti-)quark coefficients
C--   istf      (in) : 1=FL, 2=F2, 3=xF3
C--   ix        (in) : fixed ix during fill

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension def(-6:6), xx(maxn0), qq(maxn0), fvals(maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_SQF123: input address IA out of range'
      if(ispSplineType(w,ias).ne.1) stop
     +    ' SPLINT::SSP_SQF123: input address IA is not a q-spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_SQF123: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaOneD(w,ias,iat,iau,nus,iaf,iab,iac,iad)

C--   Initialise
      call smb_vfill(w(iaf),nus,0.D0)
      call smb_vfill(w(iab),nus,0.D0)
      call smb_vfill(w(iac),nus,0.D0)
      call smb_vfill(w(iad),nus,0.D0)

C--   Check iq
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)
      if(ix.le.0 .or. ix.gt.nx+1)    stop
     +    ' SPLINT::SSP_SQF123: input ix out of range'

C--   Make list of x and q2 points
      x  = xfrmix(ix)
      ia = iau-1
      do iy = 1,nus
        ia        = ia+1
        xx(iy)    = x
        qq(iy)    = exp(w(ia))
      enddo

C--   Fill fvals
      call ZSWITCH(iset)
      call ZMSTFUN(istf, def, xx, qq, fvals, nus, 1)

      call sspS1fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dsp_funs1(int ia, double z, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fdsp_funs1 FC_FUNC(dsp_funs1,DSP_FUNS1)
CXXHFW    double fdsp_funs1(int*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dsp_funs1(int ia, double z, int ichk)
CXXWRP    {
CXXWRP      return fdsp_funs1(&ia, &z, &ichk);
CXXWRP    }
C-----------------------------------------------------------------------

C     ================================================
      double precision function dsp_FunS1(ia, z, ichk)
C     ================================================

C--   Compute 1-dim spline function

      implicit double  precision (a-h,o-z)
      logical lmb_lt, lmb_le, lmb_gt

      include 'splint.inc'
      include 'spliws.inc'

      dsp_FunS1 = 0.D0
C--   Check input ia
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused)         stop
     +    ' SPLINT::DSP_FUNS1: input address IA out of range'
      if(abs(ispSplineType(w,ia)).ne.1)    stop
     +    ' SPLINT::DSP_FUNS1: input address IA is not a 1-dim spline'
      call sspSpLims(w,ia,nu,umi,uma,nv,vmi,vma,ndim,nactive)
C--   Check input z
      if(lmb_le(z,0.D0,-deps0))            stop
     +  ' SPLINT::DSP_FUNS1: input coordinate <= 0'
      u = ndim*log(z)
C--   Check input u
      if( lmb_lt(u,umi,-deps0) .or. lmb_gt(u,uma,-deps0) ) then
        if(ichk.eq.1) then
          stop ' SPLINT::DSP_FUNS1: input coordinate out of range'
        elseif(ichk.eq.0) then
          dsp_FunS1 = 0.D0
          return
        endif
      endif
C--   Spline function
      dsp_FunS1 = dspS1fun(w,ia,u)

      return
      end

C=======================================================================
C===  1-dim integration ================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    double dsp_ints1(int ia, double u1, double u2);
C-----------------------------------------------------------------------
CXXHFW  #define fdsp_ints1 FC_FUNC(dsp_ints1,DSP_INTS1)
CXXHFW    double fdsp_ints1(int*,double*,double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dsp_ints1(int ia, double u1, double u2)
CXXWRP    {
CXXWRP      return fdsp_ints1(&ia, &u1, &u2);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      double precision function dsp_IntS1(ia, u1, u2)
C     ===============================================

C--   Compute 1-dim spline integral

      implicit double  precision (a-h,o-z)
      logical lmb_lt, lmb_gt, lmb_ge

      include 'splint.inc'
      include 'spliws.inc'

      dsp_IntS1 = 0.D0
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused)         stop
     +    ' SPLINT::DSP_INTS1: input address IA out of range'
      itype = ispSplineType(w,ia)
      if(abs(itype).ne.1)                  stop
     +    ' SPLINT::DSP_INTS1: input address IA is not a 1-dim spline'
      call sspSpLims(w,ia,nu,umi,uma,nv,vmi,vma,ndim,nactive)

      if(lmb_ge(u1,u2,-deps0)) then
        dsp_IntS1 = 0.D0
        return
      endif

      if(itype.eq.-1) then
        xmi = exp(-uma)
        xma = exp(-umi)
        if(lmb_lt(u1,xmi,-deps0) .or. lmb_gt(u1,xma,-deps0)) stop
     +    ' SPLINT::DSP_INTS1: lower integration limit out of range'
        if(lmb_lt(u2,xmi,-deps0) .or. lmb_gt(u2,xma,-deps0)) stop
     +    ' SPLINT::DSP_INTS1: upper integration limit out of range'
        ymi = -log(u2)
        yma = -log(u1)
        dsp_IntS1 = dspSpIntY(w,ia,ymi,yma)

      elseif(itype.eq.1) then
        qmi = exp(umi)
        qma = exp(uma)
        if(lmb_lt(u1,qmi,-deps0) .or. lmb_gt(u1,qma,-deps0)) stop
     +    ' SPLINT::DSP_INTS1: lower integration limit out of range'
        if(lmb_lt(u2,qmi,-deps0) .or. lmb_gt(u2,qma,-deps0)) stop
     +    ' SPLINT::DSP_INTS1: upper integration limit out of range'
        tmi = log(u1)
        tma = log(u2)
        dsp_IntS1 = dspSpIntT(w,ia,tmi,tma)

      endif

      return
      end

C=======================================================================
C===  2-dim splines  ===================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    int isp_s2make(int istepx, int istepq);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_s2make FC_FUNC(isp_s2make,ISP_S2MAKE)
CXXHFW    int fisp_s2make(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_s2make(int istepx, int istepq)
CXXWRP    {
CXXWRP      return fisp_s2make(&istepx, &istepq);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===========================================
      integer function isp_S2Make(istepx, istepq)
C     ===========================================

C--   istepx  (in) : sample step in QCDNUM x-grid
C--   istepq  (in) : sample step in QCDNUM q-grid

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension ynodes(maxn0), tnodes(maxn0)

C--   Check input
      if(istepx.le.0) stop ' SPLINT::ISP_S2MAKE: istepx <= 0'
      if(istepq.le.0) stop ' SPLINT::ISP_S2MAKE: istepq <= 0'
C--   Check QCDNUM initialised
      call GETINT('vers',ivers)
      if(ivers.eq.0) stop ' SPLINT::ISP_S2MAKE: QCDNUM not initialised'
C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.0) stop
     +   ' SPLINT::ISP_S2MAKE: splint memory not initialised'

C--   Automake y-node points
      call sspYnMake(istepx, ynodes, nys, ierr)
      if(ierr.ne.0) goto 500

C--   Automake t-node points
      call sspTnMake(istepq, tnodes, nts, ierr)
      if(ierr.ne.0) goto 500

      isp_S2make = ispS2make(w, ynodes, nys, tnodes, nts)

      return

 500  continue
      write(6,*) ' '
      write(6,*) ' SPLINT::ISP_S2MAKE: Too many node points'
      write(6,*) '                     MAXN0 in splint.inc exceeded'
      stop

      end

C-----------------------------------------------------------------------
CXXHDR    int isp_s2user(double *xarr, int nx, double *qarr, int nq);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_s2user FC_FUNC(isp_s2user,ISP_s2user)
CXXHFW    int fisp_s2user(double*, int*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_s2user(double *xarr, int nx, double *qarr, int nq)
CXXWRP    {
CXXWRP      return fisp_s2user(xarr, &nx, qarr, &nq);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      integer function isp_S2User(xarr, nx, qarr, nq)
C     ===============================================

C--   xarr,nx  (in) : node points in x
C--   qarr,nq  (in) : node points in q

      implicit double precision(a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension xarr(*), qarr(*)
      dimension ynodes(maxn0), tnodes(maxn0)

C--   Check input
      if(nx.lt.2) stop ' SPLINT::ISP_S2USER: nx < 2'
      if(nq.lt.2) stop ' SPLINT::ISP_S2USER: nq < 2'
C--   Check QCDNUM initialised
      call GETINT('vers',ivers)
      if(ivers.eq.0) stop ' SPLINT::ISP_S2USER: QCDNUM not initialised'
C--   Check memory initialised
      if(iws_IsaWorkspace(w).eq.0) stop
     +   ' SPLINT::ISP_S2USER: splint memory not initialised'

C--   Setup y-nodes
      call sspYnUser(xarr, nx, ynodes, nys, ierr)
      if(ierr.ne.0) goto 500

C--   Setup t-nodes
      call sspTnUser(qarr, nq, tnodes, nts, ierr)
      if(ierr.ne.0) goto 500

      isp_S2user = ispS2make(w, ynodes, nys, tnodes, nts)

      return

 500  continue
      write(6,*) ' '
      write(6,*) ' SPLINT::ISP_S2USER: Too many node points'
      write(6,*) '                     MAXN0 in splint.inc exceeded'
      stop

      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_s2fill(int ias, double (*fun)(int*,int*,bool*),
CXXHDR                    double rs);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_s2fill FC_FUNC(ssp_s2fill,SSP_S2FILL)
CXXHFW    void fssp_s2fill(int*, double (*)(int*,int*,bool*), double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_s2fill(int ias, double (*fun)(int*,int*,bool*),
CXXWRP                    double rs)
CXXWRP    {
CXXWRP      fssp_s2fill(&ias, fun, &rs);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      subroutine ssp_S2Fill(ias, fun, rs)
C     ===================================

C--   ias     (in) : address xq-spline object
C--   fun     (in) : user defined fun(ix,iq,first)
C--   rs      (in) : roots cut (0 = no cut)

      implicit double precision(a-h,o-z)
      logical first, lmb_le

      include 'splint.inc'
      include 'spliws.inc'

      external fun

      dimension fvals(maxn0,maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_S2FILL: input address IA out of range'
      if(ispSplineType(w,ias).ne.2)  stop
     +    ' SPLINT::SSP_S2FILL: input address IA is not a 2-dim spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_S2FILL: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaTwoD(w,ias,iat,iau,nus,iav,nvs,iaF,iaC)

C--   Initialise
      call smb_vfill(w(iau+nus),nus,0.D0)
      call smb_vfill(w(iav+nvs),nvs,0.D0)
      ia1 = iws_BeginTbody(w,iaF)
      ia2 = iws_EndTbody(w,iaF)
      n   = ia2-ia1+1
      call smb_vfill(w(ia1),n,0.D0)
      ia1 = iws_BeginTbody(w,iaC)
      ia2 = iws_EndTbody(w,iaC)
      n   = ia2-ia1+1
      call smb_vfill(w(ia1),n,0.D0)

C--   Set rs cut
      if(lmb_le(rs,0.D0,-deps0)) then
        rscut = 0.D0
        sslog = 0.D0
      else
        rscut = rs
        sslog = log(rs*rs)
      endif
      call sspRangeYT(w,ias,sslog)

C--   Store rscut
      w(iat+IrsCut0) = rscut

C--   Fill fvals
      first = .true.
      do iv = 1,nvs
        q  = exp(w(iav+iv-1))
        iq = iqfrmq(q)
        nu = int(w(iav+nvs+iv-1))
        do iu = 1,nu
          x            =  exp(-w(iau+iu-1))
          ix           =  ixfrmx(x)
          fvals(iu,iv) =  fun(ix,iq,first)
          first        = .false.
        enddo
      enddo

      call sspS2fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_s2fpdf(int ias, int iset, double *def, int isel,
CXXHDR                    double rs);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_s2fpdf FC_FUNC(ssp_s2fpdf,SSP_S2FPDF)
CXXHFW    void fssp_s2fpdf(int*, int*, double*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_s2fpdf(int ias, int iset, double *def, int isel,
CXXWRP                    double rs)
CXXWRP    {
CXXWRP      fssp_s2fpdf(&ias, &iset, def, &isel, &rs);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      subroutine ssp_S2Fpdf(ias, iset, def, isel, rs)
C     ===============================================

C--   Fill xq-spline with quark linear combination

C--   ias       (in) : address xq-spline object
C--   iset      (in) : pdf set identifier
C--   def(-6:6) (in) : (anti-)quark coefficients
C--   isel      (in) : selection flag
C--   rs        (in) : roots cut (0 = no cut)

      implicit double precision(a-h,o-z)
      logical lmb_le

      include 'splint.inc'
      include 'spliws.inc'

      dimension def(-6:6), fvals(maxn0,maxn0)
      dimension xx(maxn0*maxn0), qq(maxn0*maxn0), ff(maxn0*maxn0)
      dimension ju(maxn0*maxn0), jv(maxn0*maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_S2FPDF: input address IA out of range'
      if(ispSplineType(w,ias).ne.2)  stop
     +    ' SPLINT::SSP_S2FPDF: input address IA is not a 2-dim spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_S2FPDF: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaTwoD(w,ias,iat,iau,nus,iav,nvs,iaF,iaC)

C--   Initialise
      call smb_vfill(w(iau+nus),nus,0.D0)
      call smb_vfill(w(iav+nvs),nvs,0.D0)
      ia1 = iws_BeginTbody(w,iaF)
      ia2 = iws_EndTbody(w,iaF)
      n   = ia2-ia1+1
      call smb_vfill(w(ia1),n,0.D0)
      ia1 = iws_BeginTbody(w,iaC)
      ia2 = iws_EndTbody(w,iaC)
      n   = ia2-ia1+1
      call smb_vfill(w(ia1),n,0.D0)

C--   Set rs cut
      if(lmb_le(rs,0.D0,-deps0)) then
        rscut = 0.D0
        sslog = 0.D0
      else
        rscut = rs
        sslog = log(rs*rs)
      endif
      call sspRangeYT(w,ias,sslog)

C--   Store rscut
      w(iat+IrsCut0) = rscut

C--   Make list of x and q2 points
      nn = 0
      do iv = 1,nvs
        q  = exp(w(iav+iv-1))
        nu = int(w(iav+nvs+iv-1))
        do iu = 1,nu
          x      = exp(-w(iau+iu-1))
          nn     = nn+1
          xx(nn) = x
          qq(nn) = q
          ju(nn) = iu
          jv(nn) = iv
        enddo
      enddo

      call FFLIST(iset, def, isel, xx, qq, ff, nn, 1)

C--   Copy 1-dim ff to 2-dim fvals
      do i = 1,nn
        iu           = ju(i)
        iv           = jv(i)
        fvals(iu,iv) = ff(i)
      enddo

      call sspS2fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_s2f123(int ias, int iset, double *def, int istf,
CXXHDR                    double rs);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_s2f123 FC_FUNC(ssp_s2f123,SSP_S2F123)
CXXHFW    void fssp_s2f123(int*, int*, double*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_s2f123(int ias, int iset, double *def, int istf,
CXXWRP                    double rs)
CXXWRP    {
CXXWRP      fssp_s2f123(&ias, &iset, def, &istf, &rs);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      subroutine ssp_S2F123(ias, iset, def, istf, rs)
C     ===============================================

C--   Fill xq-spline with structure function

C--   ias       (in) : address xq-spline object
C--   iset      (in) : pdf set identifier
C--   def(-6:6) (in) : (anti-)quark coefficients
C--   istf      (in) : 1=FL, 2=F2, 3=xF3, 4=FL'
C--   rs        (in) : roots cut (0 = no cut)

      implicit double precision(a-h,o-z)
      logical lmb_le

      include 'splint.inc'
      include 'spliws.inc'

      dimension def(-6:6), fvals(maxn0,maxn0)
      dimension xx(maxn0*maxn0), qq(maxn0*maxn0), ff(maxn0*maxn0)
      dimension ju(maxn0*maxn0), jv(maxn0*maxn0)

C--   Check input address
      nused = iws_WordsUsed(w)
      if(ias.le.0 .or. ias.gt.nused) stop
     +    ' SPLINT::SSP_S2F123: input address IA out of range'
      if(ispSplineType(w,ias).ne.2)  stop
     +    ' SPLINT::SSP_S2F123: input address IA is not a 2-dim spline'
      if(ispReadOnly(w,ias).eq.1)    stop
     +    ' SPLINT::SSP_S2F123: Cannot fill because spline is read-only'

C--   Get addresses
      call sspGetIaTwoD(w,ias,iat,iau,nus,iav,nvs,iaF,iaC)

C--   Initialise
      call smb_vfill(w(iau+nus),nus,0.D0)
      call smb_vfill(w(iav+nvs),nvs,0.D0)
      ia1 = iws_BeginTbody(w,iaF)
      ia2 = iws_EndTbody(w,iaF)
      n   = ia2-ia1+1
      call smb_vfill(w(ia1),n,0.D0)
      ia1 = iws_BeginTbody(w,iaC)
      ia2 = iws_EndTbody(w,iaC)
      n   = ia2-ia1+1
      call smb_vfill(w(ia1),n,0.D0)

C--   Set rs cut
      if(lmb_le(rs,0.D0,-deps0)) then
        rscut = 0.D0
        sslog = 0.D0
      else
        rscut = rs
        sslog = log(rs*rs)
      endif
      call sspRangeYT(w,ias,sslog)

C--   Store rscut
      w(iat+IrsCut0) = rscut

C--   Make list of x and q2 points
      nn = 0
      do iv = 1,nvs
        q  = exp(w(iav+iv-1))
        nu = int(w(iav+nvs+iv-1))
        do iu = 1,nu
          x      = exp(-w(iau+iu-1))
          nn     = nn+1
          xx(nn) = x
          qq(nn) = q
          ju(nn) = iu
          jv(nn) = iv
        enddo
      enddo

      call ZSWITCH(iset)
      call ZMSTFUN(istf, def, xx, qq, ff, nn, 1)

C--   Copy 1-dim ff to 2-dim fvals
      do i = 1,nn
        iu           = ju(i)
        iv           = jv(i)
        fvals(iu,iv) = ff(i)
      enddo

      call sspS2fill(w, ias, fvals)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dsp_funs2(int ia, double x, double q, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fdsp_funs2 FC_FUNC(dsp_funs2,DSP_FUNS2)
CXXHFW    double fdsp_funs2(int*,double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dsp_funs2(int ia, double x, double q, int ichk)
CXXWRP    {
CXXWRP      return fdsp_funs2(&ia, &x, &q, &ichk);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================================
      double precision function dsp_FunS2(ia, x, q, ichk)
C     ===================================================

C--   Compute 2-dim spline function

      implicit double  precision (a-h,o-z)
      logical lmb_lt, lmb_gt

      include 'splint.inc'
      include 'spliws.inc'

      dsp_FunS2 = 0.D0
C--   Check ia
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +    ' SPLINT::DSP_FUNS2: input address IA out of range'
      if(ispSplineType(w,ia).ne.2) stop
     +    ' SPLINT::DSP_FUNS2: input address IA is not a 2-dim spline'
      call sspSpLims(w,ia,nu,umi,uma,nv,vmi,vma,ndim,nactive)
C--   Check input x and q
      xmi = exp(-uma)
      xma = exp(-umi)
      qmi = exp( vmi)
      qma = exp( vma)
      icx = 0
      icq = 0
      if(lmb_lt(x,xmi,-deps0) .or. lmb_gt(x,xma,-deps0)) icx = 1
      if(lmb_lt(q,qmi,-deps0) .or. lmb_gt(q,qma,-deps0)) icq = 1
      if(icx.eq.1 .and. ichk.eq.1) stop
     +    ' SPLINT::DSP_FUNS2: x-coordinate out of range'
      if(icq.eq.1 .and. ichk.eq.1) stop
     +    ' SPLINT::DSP_FUNS2: q-coordinate out of range'
C--   Go for it
      if((icx.eq.1 .or. icq.eq.1) .and. ichk.eq.0) then
        dsp_FunS2 = 0.D0
      else
        y = -log(x)
        t =  log(q)
        dsp_FunS2 = dspS2fun(w,ia,y,t)
      endif

      return
      end

C=======================================================================
C===  2-dim integration ================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    double dsp_ints2(int ia, double x1, double x2, double q1,
CXXHDR                     double q2, double rs, int np);
C-----------------------------------------------------------------------
CXXHFW  #define fdsp_ints2 FC_FUNC(dsp_ints2,DSP_INTS2)
CXXHFW    double fdsp_ints2(int*,double*,double*,double*,double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dsp_ints2(int ia, double x1, double x2, double q1,
CXXWRP                     double q2, double rs, int np)
CXXWRP    {
CXXWRP      return fdsp_ints2(&ia, &x1, &x2, &q1, &q2, &rs, &np);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================================
      double precision function dsp_IntS2(ia, x1, x2, q1, q2, rs, np)
C     ===============================================================

C--   Compute 2-dim spline integral

      implicit double  precision (a-h,o-z)
      logical lmb_lt, lmb_gt, lmb_ge

      include 'splint.inc'
      include 'spliws.inc'

      dsp_IntS2 = 0.D0
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused)         stop
     +    ' SPLINT::DSP_INTS2: input address IA out of range'
      itype = ispSplineType(w,ia)
      if(abs(itype).ne.2)                  stop
     +    ' SPLINT::DSP_INTS2: input address IA is not a 2-dim spline'
      call sspSpLims(w,ia,nu,umi,uma,nv,vmi,vma,ndim,nactive)

      if(lmb_ge(x1,x2,-deps0) .or. lmb_ge(q1,q2,-deps0)) then
        dsp_IntS2 = 0.D0
        return
      endif

      xmi = exp(-uma)
      xma = exp(-umi)
      qmi = exp( vmi)
      qma = exp( vma)

      if(lmb_lt(x1,xmi,-deps0) .or. lmb_gt(x1,xma,-deps0)) stop
     +    ' SPLINT::DSP_INTS2: lower limit x1 out of range'
      if(lmb_lt(x2,xmi,-deps0) .or. lmb_gt(x2,xma,-deps0)) stop
     +    ' SPLINT::DSP_INTS2: upper limit x2 out of range'
      if(lmb_lt(q1,qmi,-deps0) .or. lmb_gt(q1,qma,-deps0)) stop
     +    ' SPLINT::DSP_INTS2: lower limit q1 out of range'
      if(lmb_lt(q2,qmi,-deps0) .or. lmb_gt(q2,qma,-deps0)) stop
     +    ' SPLINT::DSP_INTS2: upper limit q2 out of range'

      ymi       = -log(x2)
      yma       = -log(x1)
      tmi       =  log(q1)
      tma       =  log(q2)
      dsp_IntS2 =  dspSpIntYT(w,ia,ymi,yma,tmi,tma,rs,np,ierr)

      if(ierr.eq.1) stop
     +    ' SPLINT::DSP_INTS2: RS not compatible with RScut on spline'
      return
      end

C=======================================================================
C==  Spline extrapolation  =============================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    void ssp_extrapu(int ia, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_extrapu FC_FUNC(ssp_extrapu,SSP_EXTRAPU)
CXXHFW    void fssp_extrapu(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_extrapu(int ia, int n)
CXXWRP    {
CXXWRP      fssp_extrapu(&ia, &n);
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================
      subroutine ssp_ExtrapU(ia, n)
C     =============================

C--   Set extrapolation degree u-coordinate

C--   ia     (in): address of spline
C--   n      (in): extrapolation degree [0,3]

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +   ' SPLINT::SSP_EXTRAPU: input address IA out of range'
      if(ispSplineType(w,ia).eq.0) stop
     +   ' SPLINT::SSP_EXTRAPU: input address IA is not a spline'
C--   Check input
      if(n.lt.0 .or. n.gt.3)       stop
     +   ' SPLINT::SSP_EXTRAPU: extrapolation degree not in range [0,3]'
C--   Setit
      iatag            = iws_IaFirstTag(w,ia)
      w(iatag+NedegU0) = dble(n)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_extrapv(int ia, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_extrapv FC_FUNC(ssp_extrapv,SSP_EXTRAPV)
CXXHFW    void fssp_extrapv(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_extrapv(int ia, int n)
CXXWRP    {
CXXWRP      fssp_extrapv(&ia, &n);
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================
      subroutine ssp_ExtrapV(ia, n)
C     =============================

C--   Set extrapolation degree v-coordinate

C--   ia     (in): address of spline
C--   n      (in): extrapolation degree [0,3]

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +   ' SPLINT::SSP_EXTRAPV: input address IA out of range'
      if(ispSplineType(w,ia).eq.0) stop
     +   ' SPLINT::SSP_EXTRAPV: input address IA is not a spline'
C--   Check input
      if(n.lt.0 .or. n.gt.3)       stop
     +   ' SPLINT::SSP_EXTRAPV: extrapolation degree not in range [0,3]'
C--   Setit
      iatag            = iws_IaFirstTag(w,ia)
      w(iatag+NedegV0) = dble(n)

      return
      end

C=======================================================================
C===  Access to spline object  =========================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    int isp_splinetype(int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_splinetype FC_FUNC(isp_splinetype,ISP_SPLINETYPE)
CXXHFW    int fisp_splinetype(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_splinetype(int ia)
CXXWRP    {
CXXWRP      return fisp_splinetype(&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      integer function isp_SplineType(ia)
C     ===================================

C--   Return type +-1,2  (0 = not a spline)

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::ISP_SPLINETYPE: input address IA out of range'

      isp_SplineType = ispSplineType(w,ia)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_splims(int ia, int &nu, double &umi, double &uma,
CXXHDR                    int &nv, double &vmi, double &vma, int &nb);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_splims FC_FUNC(ssp_splims,SSP_SPLIMS)
CXXHFW    void fssp_splims(int*, int*, double*, double*, int*, double*,
CXXHFW                     double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_splims(int ia, int &nu, double &umi, double &uma,
CXXWRP                    int &nv, double &vmi, double &vma, int &nb)
CXXWRP    {
CXXWRP      fssp_splims(&ia, &nu, &umi, &uma, &nv, &vmi, &vma, &nb);
CXXWRP    }
C-----------------------------------------------------------------------

C     =========================================================
      subroutine ssp_SpLims(ia, nu, umi, uma, nv, vmi, vma, nb)
C     =========================================================

C--   Return node limits converted to x,q

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::SSP_SPLIMS: input address IA out of range'
      if(ispSplineType(w,ia).ne.0) then
        call sspSpLims(w, ia, nu, u1, u2, nv, v1, v2, ndim, nb)
        if(ndim.eq.-1) then
          umi = exp(-u2)
          uma = exp(-u1)
          vmi = 0.D0
          vma = 0.D0
        elseif(ndim.eq.1) then
          umi = exp( u1)
          uma = exp( u2)
          vmi = 0.D0
          vma = 0.D0
        elseif(ndim.eq.2) then
          umi = exp(-u2)
          uma = exp(-u1)
          vmi = exp( v1)
          vma = exp( v2)
        else
          stop ' SPLINT::SSP_SPLIMS: problem with ndim'
        endif
      else
        stop ' SPLINT::SSP_SPLIMS: input address IA is not a spline'
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_unodes(int ia, double* array, int n, int &nus);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_unodes FC_FUNC(ssp_unodes,SSP_UNODES)
CXXHFW    void fssp_unodes(int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_unodes(int ia, double* array, int n, int &nus)
CXXWRP    {
CXXWRP      fssp_unodes(&ia, array, &n, &nus);
CXXWRP    }
C-----------------------------------------------------------------------

C     ========================================
      subroutine ssp_Unodes(ia, array, n, nus)
C     ========================================

C--   Return array with u-nodes converted from y --> x or t --> q

C--   ia     (in): address of spline
C--   array (out): copy of u-nodes
C--   n      (in): dimension of array
C--   nus   (out): number of unodes copied

C--   The routine sets the remainder of array to zero

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension array(*)

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::SSP_UNODES: input address IA out of range'
      if(ispSplineType(w,ia).ne.0) then
        itag       = iws_IaFirstTag(w,ia)
        nus        = int(w(itag+NwUtab0))
        ndim       = int(w(itag+IdimSp0))
        if(nus.gt.n) stop
     +      ' SPLINT::SSP_UNODES: insufficient output array size'
        itabu = int(w(itag+IsUtab0))+ia
        iau   = iws_BeginTbody(w,itabu)-1
        if(ndim.eq.-1 .or. ndim.eq.2) then
          do i = 1,nus
            array(nus+1-i) = exp(-w(iau+i))
          enddo
          do i = nus+1,n
            array(i)       = 0.D0
          enddo
        else
          do i = 1,nus
            array(i)       = exp( w(iau+i))
          enddo
          do i = nus+1,n
            array(i)       = 0.D0
          enddo
        endif
      else
        stop ' SPLINT::SSP_UNODES: input address IA is not a spline'
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_vnodes(int ia, double* array, int n, int &nvs);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_vnodes FC_FUNC(ssp_vnodes,SSP_VNODES)
CXXHFW    void fssp_vnodes(int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_vnodes(int ia, double* array, int n, int &nvs)
CXXWRP    {
CXXWRP      fssp_vnodes(&ia, array, &n, &nvs);
CXXWRP    }
C-----------------------------------------------------------------------

C     ========================================
      subroutine ssp_Vnodes(ia, array, n, nvs)
C     ========================================

C--   Return array with v-nodes converted from t --> q

C--   ia     (in): address of spline
C--   array (out): copy of v-nodes (all zero for 1-dim spline)
C--   n      (in): dimension of array
C--   nvs   (out): number of vnodes copied (0 for 1-dim spline)

C--   The routine sets the remainder of array to zero

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dimension array(*)

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::SSP_VNODES: input address IA out of range'
      if(ispSplineType(w,ia).ne.0) then
        itag       = iws_IaFirstTag(w,ia)
        nvs        = int(w(itag+NwVtab0))
        if(nvs.gt.n) then
          stop ' SPLINT::SSP_VNODES: insufficient output array size'
        elseif(nvs.eq.0) then
          do i = 1,n
            array(i) = 0.D0
          enddo
        else
          itabv = int(w(itag+IsVtab0))+ia
          iav   = iws_BeginTbody(w,itabv)-1
          do i = 1,nvs
            array(i) = exp(w(iav+i))
          enddo
          do i = nvs+1,n
            array(i) = 0.D0
          enddo
        endif
      else
        stop ' SPLINT::SSP_VNODES: input address IA is not a spline'
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_nprint(int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_nprint FC_FUNC(ssp_nprint,SSP_NPRINT)
CXXHFW    void fssp_nprint(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_nprint(int ia)
CXXWRP    {
CXXWRP      fssp_nprint(&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     =========================
      subroutine ssp_Nprint(ia)
C     =========================

C--   Print list of node points

C--   ia     (in): address of spline

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::SSP_NPRINT: input address IA out of range'
      if(ispSplineType(w,ia).eq.0) stop
     +   ' SPLINT::SSP_NPRINT: input address IA is not a spline'

      itag       = iws_IaFirstTag(w,ia)
      idim       = int(w(itag+IdimSp0))
      nx         = 0
      nq         = 0
      inx        = 0
      inq        = 0
      if(idim.eq.-1) then
        ixtab    = int(w(itag+IsUtab0)) + ia
        nx       = int(w(itag+NwUtab0))
        iax      = iws_BeginTbody(w,ixtab)
        nq       = 0
        iaq      = 0
      elseif(idim.eq.1) then
        nx       = 0
        iax      = 0
        iqtab    = int(w(itag+IsUtab0)) + ia
        nq       = int(w(itag+NwUtab0))
        iaq      = iws_BeginTbody(w,iqtab)-1
      elseif(idim.eq.2) then
        ixtab    = int(w(itag+IsUtab0)) + ia
        nx       = int(w(itag+NwUtab0))
        iax      = iws_BeginTbody(w,ixtab)
        inx      = iax+nx
        iqtab    = int(w(itag+IsVtab0)) + ia
        nq       = int(w(itag+NwVtab0))
        iaq      = iws_BeginTbody(w,iqtab)-1
        inq      = iaq+nq
      endif
      imax = max(nx,nq)
      write(6,
     +   '(/''  N   IX     XNODE    NQMA      IQ     QNODE    NXMI'')')
      do i = 1,imax
C--     Both column x and colomn q
        if(i.le.nx .and. i.le.nq) then
          xx = exp(-w(iax+nx-i))
          ix = ixfrmx(xx)
          qq = exp(w(iaq+i))
          iq = iqfrmq(qq)
          jx = 1
          jq = nq
          if(idim.eq.2) then
            jq = int(w(inx+nx-i))
            jy = int(w(inq+i))
            if(jy.eq.0) then
              jx = 0
            else
              jx = nx-jy+1
            endif
          endif
          write(6,'(I3,I5,E13.5,I5,I8,E13.5,I5)') i,ix,xx,jq,iq,qq,jx
C--     Only column x
        elseif(i.le.nx .and. i.gt.nq) then
          xx = exp(-w(iax+nx-i))
          ix = ixfrmx(xx)
          jq = nq
          if(idim.eq.2) jq = int(w(inx+nx-i))
          write(6,
     +   '(I3,I5,E13.5,I5,''       -       -         -'')') i,ix,xx,jq
C--     Only column q
        elseif(i.gt.nx .and. i.le.nq) then
          qq = exp(w(iaq+i))
          iq = iqfrmq(qq)
          jx = 1
          if(idim.eq.2) then
            jy = int(w(inq+i))
            if(jy.eq.0) then
              jx = 0
            else
             jx = nx-jy+1
            endif
          endif
          write(6,
     +   '(I3,''    -       -         -   '',I5,E13.5,I5)') i,iq,qq,jx
        endif
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dsp_rscut(int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fdsp_rscut FC_FUNC(dsp_rscut,DSP_RSCUT)
CXXHFW    double fdsp_rscut(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dsp_rscut(int ia)
CXXWRP    {
CXXWRP      return fdsp_rscut(&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     =======================================
      double precision function dsp_RsCut(ia)
C     =======================================

C--   Get roots cut

C--   ia     (in): address of spline

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      dsp_RsCut = 0.D0

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::DSP_RSCUT: input address IA out of range'
      if(ispSplineType(w,ia).ne.2) stop
     +   ' SPLINT::DSP_RSCUT: input address IA is not a 2-dim spline'

C--   Get rscut
      itag      = iws_IaFirstTag(w,ia)
      dsp_RsCut = w(itag+IrsCut0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dsp_rsmax(int ia, double rs);
C-----------------------------------------------------------------------
CXXHFW  #define fdsp_rsmax FC_FUNC(dsp_rsmax,DSP_RSMAX)
CXXHFW    double fdsp_rsmax(int*,double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dsp_rsmax(int ia, double rs)
CXXWRP    {
CXXWRP      return fdsp_rsmax(&ia, &rs);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==========================================
      double precision function dsp_RsMax(ia,rs)
C     ==========================================

C--   Get max roots cut

C--   ia     (in): address of spline

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      logical lmb_le

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::DSP_RSMAX: input address IA out of range'
      if(ispSplineType(w,ia).ne.2) stop
     +   ' SPLINT::DSP_RSMAX: input address IA is not a 2-dim spline'

C--   Get max rscut
      if(lmb_le(rs,0.D0,-deps0)) then
        dsp_RsMax = 0.D0
      else
        dsp_RsMax = dspRsMax(w,ia,log(rs*rs))
      endif

      return
      end

C=======================================================================
C===  Memory ===========================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    void ssp_mprint();
C-----------------------------------------------------------------------
CXXHFW  #define fssp_mprint FC_FUNC(ssp_mprint,SSP_MPRINT)
CXXHFW    void fssp_mprint();
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_mprint()
CXXWRP    {
CXXWRP      fssp_mprint();
CXXWRP    }
C-----------------------------------------------------------------------

C     =======================
      subroutine ssp_Mprint()
C     =======================

C--   Print memory tree

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      call sws_WsTree(w)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int isp_spsize(int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_spsize FC_FUNC(isp_spsize,ISP_SPSIZE)
CXXHFW    int fisp_spsize(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_spsize(int ia)
CXXWRP    {
CXXWRP      return fisp_spsize(&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================
      integer function isp_SpSize(ia)
C     ===============================

C--   Return object size

C--   ia (in): address of spline; 0 --> memory size; 1 --> words used

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      ntot = iws_SizeOfW(w)
      if(ia.eq.0) then
        isp_SpSize = ntot
      elseif(ia.eq.1) then
        isp_SpSize = iws_WordsUsed(w)+1         !count also trailer word
      elseif(ia.lt.0 .or. ia.gt.ntot) then
        stop ' SPLINT::ISP_SPSIZE: input address IA out of range'
      elseif(ispSplineType(w,ia).ne.0) then
        isp_SpSize = iws_ObjectSize(w,ia)
      else
        isp_SpSize = 0
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_erase(int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_erase FC_FUNC(ssp_erase,SSP_ERASE)
CXXHFW    void fssp_erase(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_erase(int ia)
CXXWRP    {
CXXWRP      fssp_erase(&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ========================
      subroutine ssp_Erase(ia)
C     ========================

C--   Wipe memory from ia onwards

C--   ia (in): address of spline;  ia = 1 --> wipe all splines
C--
C--   NB: 1st table-set may be user store and not 1st spline in memory
C--       Therefore the address of the 1st spline is stored in tag field

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) then
        stop  ' SPLINT::SSP_ERASE: input address IA out of range'
      endif
C--   Get address of first spline in memory
      iaR  = iws_IaRoot()
      itag = iws_IaFirstTag(w,iaR)
      ias1 = int(w(itag+IaSfirst0))
C--   No splines in memory
      if(ias1.eq.0) return
C--   Re-set address if necessary
      if(ia.eq.1) then
        ja = ias1
      else
        ja = ia
      endif
C--   Wipe
      if(ispSplineType(w,ja).ne.0) then
        call sws_WsWipe(w,ja)
      else
        stop ' SPLINT::SSP_ERASE: input address IA is not a spline'
      endif
C--   No first spline in memory after wipe
      if(ja.eq.ias1) w(itag+IaSfirst0) = 0.D0

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_spdump(int ia, string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_spdumpcpp FC_FUNC(ssp_spdumpcpp,SSP_SPDUMPCPP)
CXXHFW    void fssp_spdumpcpp(int*,char*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_spdump(int ia, string fname)
CXXWRP    {
CXXWRP      int ls = fname.size();
CXXWRP      char *cfname = new char[ls+1];
CXXWRP      strcpy(cfname,fname.c_str());
CXXWRP      fssp_spdumpcpp(&ia, cfname, &ls);
CXXWRP      delete[] cfname;
CXXWRP    }
C-----------------------------------------------------------------------

C     =======================================
      subroutine ssp_SpDumpCPP(ia, fname, ls)
C     =======================================

      implicit double precision (a-h,o-z)

      character*(100) fname

      if(ls.gt.100) stop
     +            'SPLINT::SSP_SPDUMP: input file name > 100 characters'

      call ssp_SpDump(ia, fname)

      return
      end

C     ================================
      subroutine ssp_SpDump(ia, fname)
C     ================================

      implicit double precision (a-h,o-z)

      character*(*) fname

      include 'splint.inc'
      include 'spliws.inc'

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::SSP_SPDUMP: input address IA out of range'
      if(ispSplineType(w,ia).eq.0) stop
     +      ' SPLINT::SSP_SPDUMP: input address IA is not a spline'

C--   Dump spline
      call sws_TsDump(fname,keyval0,w,ia,ierr)

      if(ierr.ne.0) stop
     +      ' SPLINT::SSP_SPDUMP: cannot open or write output file'

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int isp_spread(string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fisp_spreadcpp FC_FUNC(isp_spreadcpp,ISP_SPREADCPP)
CXXHFW    int fisp_spreadcpp(char*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int isp_spread(string fname)
CXXWRP    {
CXXWRP      int ls = fname.size();
CXXWRP      char *cfname = new char[ls+1];
CXXWRP      strcpy(cfname,fname.c_str());
CXXWRP      int ia = fisp_spreadcpp(cfname, &ls);
CXXWRP      delete[] cfname;
CXXWRP      return ia;
CXXWRP    }
C-----------------------------------------------------------------------

C     =========================================
      integer function isp_SpReadCPP(fname, ls)
C     =========================================

      implicit double precision (a-h,o-z)

      character*(100) fname

      if(ls.gt.100) stop
     +            'SPLINT::SSP_SPREAD: input file name > 100 characters'

      isp_SpReadCPP = isp_SpRead(fname)

      return
      end

C     ==================================
      integer function isp_SpRead(fname)
C     ==================================

      implicit double precision (a-h,o-z)

      character*(*) fname

      include 'splint.inc'
      include 'spliws.inc'

C--   Read spline
      ia = iws_TsRead(fname,keyval0,w,ierr)
      if(ierr.eq.-1) then
        stop  ' SPLINT::ISP_SPREAD: cannot open or read input file'
      elseif(ierr.eq.-2) then
        stop  ' SPLINT::ISP_SPREAD: incompatible or obsolete input file'
      endif
C--   Set readonly flag in spline tag field
      iatag = iws_IaFirstTag(w, ia)
      w(iatag+IrOnly0) = dble(1)
C--   Store address of first spline
      iaR = iws_iaRoot()
      iat = iws_IaFirstTag(w, iaR)
      if(int(w(iat+IaSfirst0)).eq.0) w(iat+IaSfirst0) = dble(ia)
C--   Return spline address
      isp_SpRead = ia

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void ssp_spsetval(int ia, int i, double val);
C-----------------------------------------------------------------------
CXXHFW  #define fssp_spsetval FC_FUNC(ssp_spsetval,SSP_SPSETVAL)
CXXHFW    void fssp_spsetval(int*,int*,double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void ssp_spsetval(int ia, int i, double val)
CXXWRP    {
CXXWRP      fssp_spsetval(&ia, &i, &val);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      subroutine ssp_SpSetVal(ia, i, val)
C     ===================================

C--   Store value in a spline ia

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::SSP_SPSETVAL: input address IA out of range'
      if(ispSplineType(w,ia).eq.0) stop
     +      ' SPLINT::SSP_SPSETVAL: input address IA is not a spline'
C--   Get storage address
      ja = ispIaFromI(w, ia, i)
C--   Index out of range
      if(ja.eq.0) stop
     +      ' SPLINT::SSP_SPSETVAL: index I not in range'
C--   Set value
      w(ja) = val

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dsp_spgetval(int ia, int i);
C-----------------------------------------------------------------------
CXXHFW  #define fdsp_spgetval FC_FUNC(dsp_spgetval,DSP_SPGETVAL)
CXXHFW    double fdsp_spgetval(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dsp_spgetval(int ia, int i)
CXXWRP    {
CXXWRP      return fdsp_spgetval(&ia, &i);
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================================
      double precision function dsp_SpGetVal(ia, i)
C     =============================================

C--   Get value from a spline ia

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      include 'spliws.inc'

C--   Return storage space available
      if(ia.eq.0 .and. i.eq.0) then
        dsp_SpGetVal = dble(nusr0)
        return
      endif

C--   Check ia is in range
      nused = iws_WordsUsed(w)
      if(ia.le.0 .or. ia.gt.nused) stop
     +      ' SPLINT::SSP_SPGETVAL: input address IA out of range'
      if(ispSplineType(w,ia).eq.0) stop
     +      ' SPLINT::SSP_SPGETVAL: input address IA is not a spline'

C--   Get storage address
      ja = ispIaFromI(w, ia, i)
C--   Index out of range
      if(ja.eq.0) stop
     +      ' SPLINT::SSP_SPGETVAL: index I not in range'
C--   Get value
      dsp_SpGetVal = w(ja)

      return
      end
