***********************************************************************
*                                                                     *
*     ZMSTF package   Author: Michiel Botje, h24@nikhef.nl            *
*                                                                     *
*     QCDNUM add-on that calculates zero-mass structure functions     *
*     F2, FL, xF3 and FL' from QCDNUM evolved pdfs. The package       *
*     contains fast routines for general use and slow routines just   *
*     to illustrate how to use the QCDNUM convolution engine in       *
*     prototype code for structure function calculations.             *
*                                                                     *
***********************************************************************

C--   file zmstfs.f containing the zmstf structure function routines

C--   integer function izmvers()
C--   subroutine zmwords(ntot,nused)

C--   subroutine zmdefq2(a,b)
C--   subroutine zmabval(a,b)
C--   double precision function zmqfrmu(qmu2)
C--   double precision function zmufrmq(q2)
C--   logical function LzmRvar(epsi)
C--   logical function LzmQvar(epsi)

C--   subroutine zswitch(jset)

C--   double precision function zmstfij(istf,def,ix,iq,ichk)
C--   subroutine zmstfun(istf,qvec,x,q,f,n,ichk)
C--   subroutine zgetstf(istf,iord,iscale,qvec,x,qmu,f,n,ichk)
C--   double precision function zfacL(ix,iq,nf,ithrs)
C--   double precision function zfacL2(ix,iq,nf,ithrs)
C--   subroutine zmSteer(istf,ityp,iord,iscale,idmat,ifound)
C--   subroutine zmSetupL0(istf,ityp,idwt,iwrk)
C--   subroutine zmSetupL1(istf,ityp,idwt,iwrk)
C--   subroutine zmSetupL2(istf,ityp,idwt,iwrk)

C--   subroutine zmslowf(istf,def,x,q,f,n,ichk)
C--   subroutine zselect(isel,lquark,lgloun,ichk)
C--   double precision function dzmF2ij(ix,iq)
C--   double precision function dzmFLij(ix,iq)
C--   double precision function dzmF3ij(ix,iq)

C=====================================================================

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  ZMSTF routines                                          */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  ZMSTF routines                                            */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  ZMSTF routines                                            */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    int izmvers();
C-----------------------------------------------------------------------
CXXHFW  #define fizmvers FC_FUNC(izmvers,IZMVERS)
CXXHFW    int fizmvers();
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int izmvers()
CXXWRP    {
CXXWRP      return fizmvers();
CXXWRP    }
C-----------------------------------------------------------------------

C     ==========================
      integer function izmvers()
C     ==========================

      implicit double precision(a-h,o-z)

      include 'zmstf.inc'

      izmvers = iivers

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void zmwords(int &ntotal, int &nused);
C-----------------------------------------------------------------------
CXXHFW  #define fzmwords FC_FUNC(zmwords,ZMWORDS)
CXXHFW    void fzmwords(int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zmwords(int &ntotal, int &nused)
CXXWRP  {
CXXWRP    fzmwords(&ntotal,&nused);
CXXWRP  }
C-----------------------------------------------------------------------
 
C     ============================== 
      subroutine zmwords(ntot,nused)
C     ==============================

C--   Number of words available, and used

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'zmstf.inc'
      include 'zmstore.inc'
      
      ntot = nzmstor
      
      if(izini.eq.12345) then
        nused = nzused
      else
        nused = 0
      endif
      
      return        
      end

C=======================================================================
C==   Scale variations  ================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    void zmdefq2(double a, double b);
C-----------------------------------------------------------------------
CXXHFW  #define fzmdefq2 FC_FUNC(zmdefq2,ZMDEFQ2)
CXXHFW    void fzmdefq2(double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zmdefq2(double a, double b)
CXXWRP  {
CXXWRP    fzmdefq2(&a,&b);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================
      subroutine zmdefq2(a,b)
C     =======================

C--   Q2 = a*mu2 + b

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'

      call setUmsg('ZMDEFQ2')
      
      if(izini.ne.12345) stop
     + 'ZMDEFQ2: please first call ZMFILLW or ZMREADW --> STOP'

C--   Limit allowed range
      if(a.lt.0.1D0 .or. a.gt.10.D0) then
         stop 'ZMDEFQ2: Coefficient A outside range [0.1,10] --> STOP'
      endif
      if(abs(b).gt.100.D0) then
        stop 
     +      'ZMDEFQ2: Coefficient B outside range [-100,100] --> STOP'
      endif                             
      
      ascale = a
      bscale = b
      
C--   Tolerance
      call GetVal('epsi',epsi)
C--   Yes/no scale change thats the question
      jscale = 0
      if(abs(a-1.D0).gt.epsi) jscale = 1
      if(abs(b)     .gt.epsi) jscale = 1

      call clrUmsg
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void zmabval(double &a, double &b);
C-----------------------------------------------------------------------
CXXHFW  #define fzmabval FC_FUNC(zmabval,ZMABVAL)
CXXHFW    void fzmabval(double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zmabval(double &a, double &b)
CXXWRP  {
CXXWRP    fzmabval(&a,&b);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     =======================
      subroutine zmabval(a,b)
C     =======================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'
      
      if(izini.ne.12345) stop
     + 'ZMABVAL: please first call ZMFILLW or ZMREADW --> STOP'
      
      a = ascale
      b = bscale
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    double zmqfrmu(double qmu2);
C-----------------------------------------------------------------------
CXXHFW  #define fzmqfrmu FC_FUNC(zmqfrmu,ZMQFRMU)
CXXHFW    double fzmqfrmu(double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double zmqfrmu(double qmu2)
CXXWRP  {
CXXWRP    return fzmqfrmu(&qmu2);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================================
      double precision function zmqfrmu(qmu2)
C     =======================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'
      
      if(izini.ne.12345) 
     +    stop 'ZMQFRMU: ZMSTF not initialized --> STOP'
      
      zmqfrmu = ascale*qmu2 + bscale
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    double zmufrmq(double Q2);
C-----------------------------------------------------------------------
CXXHFW  #define fzmufrmq FC_FUNC(zmufrmq,ZMUFRMQ)
CXXHFW    double fzmufrmq(double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double zmufrmq(double Q2)
CXXWRP  {
CXXWRP    return fzmufrmq(&Q2);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     =====================================
      double precision function zmufrmq(q2)
C     =====================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'
      
      if(izini.ne.12345) stop
     + 'ZMUFRMQ: please first call ZMFILLW or ZMREADW --> STOP'
      
      zmufrmq = (q2-bscale)/ascale
      
      return
      end

C     ==============================
      logical function LzmRvar(epsi)
C     ==============================

C--   True if renormalization scale is varied (within tolerance epsi)

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'

      if(izini.ne.12345) stop
     + 'LZMRVAR: please first call ZMFILLW or ZMREADW --> STOP'

C--   Renormalization scale
      call getabr(ar,br)
C--   Renor scale and fact scale equal or not thats the question
      if(abs(ar-1.D0).gt.epsi .or.
     +   abs(br     ).gt.epsi      ) then
        LzmRvar = .true.
      else
        LzmRvar = .false.
      endif
      
      return
      end
      
C     ==============================
      logical function LzmQvar(epsi)
C     ==============================

C--   True if Q2 scale is varied (within tolerance epsi)

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'
      
      if(izini.ne.12345) stop
     + 'LZMQVAR: please first call ZMFILLW or ZMREADW --> STOP'

C--   Q2 scale and fact scale equal or not thats the question
      if(abs(ascale-1.D0).gt.epsi .or.
     +   abs(bscale     ).gt.epsi      ) then
        LzmQvar = .true.
      else
        LzmQvar = .false.
      endif      
      
      return
      end

C=======================================================================
C==   Switch pdf set  ==================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    void zswitch(int iset);
C-----------------------------------------------------------------------
CXXHFW  #define fzswitch FC_FUNC(zswitch,ZSWITCH)
CXXHFW    void fzswitch(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zswitch(int iset)
CXXWRP  {
CXXWRP    fzswitch(&iset);
CXXWRP  }
C-----------------------------------------------------------------------

C     ========================
      subroutine zswitch(jset)
C     ========================

C--   Pdf set identifier

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'

      dimension pars(13)

      call setUmsg('ZSWITCH')

      if(izini.ne.12345) stop 
     + 'ZSWITCH: please first call ZMFILLW or ZMREADW --> STOP'

      call getInt('mset',mset)
      
      if(jset.lt.1 .or. jset.gt.mset) stop
     + 'ZSWITCH: iset not in range [1,mset] --> STOP'

      call cpypar(pars,13,jset)
      ityp = int(pars(13))

      if(ityp.eq.2) stop
     + 'ZSWITCH: cannot handle polarised pdfs --> STOP'
      if(ityp.eq.3) stop
     + 'ZSWITCH: cannot handle fragmentation functions --> STOP'
      if(ityp.eq.5) stop
     + 'ZSWITCH: cannot handle user-defined pdf set --> STOP'
       
      izpdf = jset

      call clrUmsg

      return
      end

C=======================================================================
C==   Fast structure functions   =======================================
C=======================================================================

C     ======================================================
      double precision function zmstfij(jstf,def,ix,iq,ichk)
C     ======================================================

C--   Compute structure function on grid point

C--   jstf        (in): 1=FL, 2=F2, 3=xF3, 4=FL'
C--   def(-6:6)   (in): quark linear combination coefficients
C--   ix          (in): x-grid point
C--   iq          (in): mu2-grid point
C--   ichk        (in): 0   return null value if ix or iq outside grid
C--                     1   fatal error if any ix or iq outside grid

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'

      logical lquark, lgluon, LzmQvar
      common /qgflags/ lquark(3),lgluon(3)
      common /pdfdefs/ qvec(-6:6)

      dimension def(-6:6)

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('ZMSTFIJ')
      
C--   Check ZMSTF initialized
      if(izini.ne.12345) stop
     + 'ZMSTFIJ: please first call ZMFILLW or ZMREADW --> STOP'

C--   Decode jstf
      iset = jstf/10
      if(iset.ne.0) call ZSWITCH(iset)
      istf = jstf-10*iset

C--   Check ix, iq in range
      call grpars(nx, x1, x2, nq, q1, q2, io)
      if(ix.eq.nx+1) then
        zmstfij = 0.D0                                      !catch x = 1
      elseif( (ix.le.0.or.ix.gt.nx) .or.
     +        (iq.le.0.or.iq.gt.nq) ) then
        if(ichk.eq.1) then
          stop 'ZMSTFIJ: IX or IQ out of range'
        else
          zmstfij = 0.D0
          call clrUmsg
          return
        endif
      endif
C--   Tolerance
      call GetVal('epsi',epsi)
C--   No Q2 scale variation, thank you
      if(LzmQvar(epsi)) then
        stop 'ZMSTFIJ: You cant vary Q2 scale --> STOP'
      endif

C--   Start scoping (not really necessary here)
      call IdScope(0.D0,izpdf)
C--   Copy to common block
      do i = -6,6
        qvec(i) = def(i)
      enddo
C--   Select what to calculate
      do iord = 1,3
        lquark(iord) = .true.
        lgluon(iord) = .true.
      enddo
C--   Go for stf
      if(istf.eq.1)     then
        zmstfij = dzmFLij(ix,iq)
      elseif(istf.eq.2) then
        zmstfij = dzmF2ij(ix,iq)
      elseif(istf.eq.3) then
        zmstfij = dzmF3ij(ix,iq)
      elseif(istf.eq.4) then
        zmstfij = dzmFpij(ix,iq)
      else
        zmstfij = 0.D0
        stop 'ZMSTFIJ: invalid structure function label --> STOP'
      endif
C--   End scoping
      call IdScope(0.D0,-izpdf)
      call clrUmsg
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void zmstfun(int istf, double *def, double *x, double *Q2, double *f, int n, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fzmstfun FC_FUNC(zmstfun,ZMSTFUN)
CXXHFW    void fzmstfun(int*, double*, double*, double*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zmstfun(int istf, double *def, double *x, double *Q2, double *f, int n, int ichk)
CXXWRP  {
CXXWRP    fzmstfun(&istf,def,x,Q2,f,&n,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================================
      subroutine zmstfun(jstf,qvec,x,q,f,n,ichk)
C     ==========================================

C--   Fast structure function calculation
C--
C--   jstf       (in)   1,2,3,4 = FL, F2, xF3, FL'
C--   qvec(-6:6) (in)   quark linear combination
C--   x(i)       (in)   list of x values
C--   q(i)       (in)   list of Q2 values
C--   f(i)       (out)  list of interpolated structure functions
C--   n          (in)   number of items in x, q and f
C--   ichk       (in)   0   return null value if x,q outside grid
C--                     1   fatal error if any x or q outside grid

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'
      
      logical LzmRvar, LzmQvar

      dimension qvec(-6:6),x(*),q(*),f(*)

      save jchkrem
      data jchkrem/999999999/
      
C--   Set this the same as mpt0 in qcdnum.inc (not critical)
      parameter( nzpnts = 5000 )
      dimension xx(nzpnts),qm(nzpnts)

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('ZMSTFUN')
      
C--   Check ZMSTF initialized      
      if(izini.ne.12345) stop
     + 'ZMSTFUN: please first call ZMFILLW or ZMREADW --> STOP'

C--   Max number of interpolations
      call getint('mpt0',nmax)
      nmax = min(nmax,nzpnts)

C--   Decode jstf
      iset = jstf/10
      if(iset.ne.0) call ZSWITCH(iset)
      istf = jstf-10*iset
               
C--   Check input
      if(istf.le.0 .or. istf.gt.4) then
        stop 'ZMSTFUN: ISTF not in range [1-4] --> STOP'
      endif
      if(n.le.0) then
        call getInt('lunq',lun)
        write(lun,'(
     +    '' ZMSTFUN: N = number of points .le. zero --> STOP'')')
        stop
      endif

C--   Force evolution parameters of set izpdf
      call IdScope(0.D0,izpdf)
C--   QCD order
      call getord(iord)
C--   Tolerance
      call GetVal('epsi',epsi)
C--   No renor and Q2 scale variation, thank you
      if(LzmRvar(epsi) .and. LzmQvar(epsi)) then
        stop 'ZMSTFUN: You cant vary both Q2 and muR2 scales --> STOP'
      endif
C--   jscale is set in common block zmscale.inc by s/r zmDefQ2
      iscale = jscale
      
C--   Fill output array f in batches of nmax words
      ipt = 0
      jj  = 0      
      do i = 1,n
        ipt     = ipt+1
        xx(ipt) = x(i)
        qm(ipt) = zmufrmq(q(i))
        if(ipt.eq.nmax) then
          call zgetstf
     +             (istf,iord,iscale,qvec,xx,qm,f(jj*nmax+1),ipt,ichk)
          ipt = 0
          jj  = jj+1
        endif
      enddo
C--   Flush remaining ipt points
      if(ipt.ne.0) then
        call zgetstf
     +             (istf,iord,iscale,qvec,xx,qm,f(jj*nmax+1),ipt,ichk)
      endif

C--   End scoping
      call IdScope(0.D0,-izpdf)
      call clrUmsg

      return
      end
      
C     ========================================================
      subroutine zgetstf(istf,iord,iscale,qvec,x,qmu,f,n,ichk)
C     ========================================================

C--   Fast structure function calculation
C--
C--   istf       (in)   1,2,3,4 = FL, F2, xF3, FL'
C--   iord       (in)   1,2,3   = LO, NLO, NNLO
C--   iscale     (in)   if 0 then Q2 = mu2, otherwise Q2 # mu2
C--   qvec(-6:6) (in)   quark linear combination
C--   x(i)       (in)   list of x values
C--   qmu(i)     (in)   list of mu2 values
C--   f(i)       (out)  list of interpolated structure functions
C--   n          (in)   number of items in x, q and f
C--   ichk       (in)   0   return null value if x,q outside grid
C--                     1   fatal error if any x or q outside grid

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'
      
      dimension qvec(-6:6),x(*),qmu(*),f(*)
      dimension iweit(4), idmat(5,0:2)

      save isel
      dimension isel(6)
C--               GL QQ SI NS NP VA+NM
C--                1  2  3  4  5   6
      data isel /  0, 7, 1, 6, 2,  5  /
      
      external zfacL, zfacL2

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(isel)

C--   Pdf type unpolarized
      jset     = izpdf

C--   Initialize
      call fastIni(x,qmu,n,ichk)

C--   Wipe scratch buffers
      call fastClr(0)
C--   Loop over all possible pdf types (GL, QQ, SI, NS, ....)
      do ityp = 1,6
        call zmSteer(istf,ityp,iord,iscale,idmat,ifound)
C--     Pdf found for this structure function
        if(ifound.ne.0) then
C--       Copy ityp = (GL,QQ,SI,NS,NP,VA+NM) to scratch buffer 1
          call fastSns(jset,qvec,isel(ityp),1)
C--       Convolute with perturbative expansion of coeff functions
          do i = 1,3
            if(idmat(i,0) .ne. 0) then
              iweit(i) = 1000*jsetw+idmat(i,0)
            else
              iweit(i) = 0
            endif
          enddo
          iweit(4) = idmat(4,0)
          call fastFxK(zmstor,iweit,1,2)
C--       Accumulate structure function
          call fastCpy(2,3,1)
C--       Scale dependent terms
          if(idmat(5,1).ne.0) then
C--         L term
            do i = 1,3
              if(idmat(i,1) .ne. 0) then
                iweit(i) = 1000*jsetw+idmat(i,1)
              else
                iweit(i) = 0
              endif
            enddo
            iweit(4) = idmat(4,1)
            call fastFxK(zmstor,iweit,1,2)
            call fastKin(2,zfacL)
            call fastCpy(2,3,1)
          endif 
          if(idmat(5,2).ne.0) then
C--         L2 term
            do i = 1,3
              if(idmat(i,2) .ne. 0) then
                iweit(i) = 1000*jsetw+idmat(i,2)
              else
                iweit(i) = 0
              endif
            enddo
            iweit(4) = idmat(4,2)
            call fastFxK(zmstor,iweit,1,2)
            call fastKin(2,zfacL2)
            call fastCpy(2,3,1)
          endif
        endif
      enddo

C--   Interpolate
      call fastFxq(3,f,n)

      return
      end
            
C     ===============================================
      double precision function zfacL(ix,iq,nf,ithrs)
C     ===============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      jx    = ix           !avoid compiler warning
      jf    = nf           !avoid compiler warning
      jthrs = ithrs        !avoid compiler warning
      qmu   = qfrmiq(iq)
      qsq   = zmqfrmu(qmu)
      zfacL = log(qsq/qmu)
      
      return
      end
      
C     ================================================
      double precision function zfacL2(ix,iq,nf,ithrs)
C     ================================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      fac    = zfacL(ix,iq,nf,ithrs)
      zfacL2 = fac*fac
      
      return
      end
      
C     ======================================================
      subroutine zmSteer(istf,ityp,iord,iscale,idmat,ifound)
C     ======================================================

C--   Figures out if there is something to do for a given
C--   combination of structure function and pdf-type.
C--   This depends on the order of the calculation and 
C--   also if the Q2 scale is different from the factorization
C--   scale. Therefore this info is passed too.
C--
C--   Input:   istf         1=FL, 2=F2, 3=xF3, 4=FL'
C--            ityp         1=gluon, 2=all quarks, 3=singlet, 
C--                         4=all non-singlet, 5=ns+, 6=va+ns-
C--            iord         1=LO, 2=NLO, 3=NNLO
C--            iscale       if 0 then Q2 = mu2, otherwise Q2 # mu2
C--   Output   idmat(5,0:2) j=0,1,2  L^0, L^1, L^2 terms
C--                         i=1,2,3  LO, NLO, NNLO wtable ids
C--                         i=4      Leading power of alphas
C--                         i=5      Set to 0 if no L^j terms
C--            ifound       Set to zero if there is nothing to do 

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      dimension idmat(5,0:2),idwt(4)
      
      ifound = 0
      do j = 0,2
        do i = 1,5
          idmat(i,j) = 0
        enddo
      enddo    
      
      call zmSetupL0(istf,ityp,idwt,iwrk)
      ifound = max(ifound,iwrk)
      do i = 1,4
        idmat(i,0) = idwt(i)
      enddo  
      idmat(5,0) = iwrk
      if(iscale.eq.0) return
      if(iord.eq.1)   return
      
      call zmSetupL1(istf,ityp,idwt,iwrk)
      ifound = max(ifound,iwrk)
      do i = 1,4
        idmat(i,1) = idwt(i)
      enddo
      idmat(5,1) = iwrk
      if(iord.eq.2) return
      
      call zmSetupL2(istf,ityp,idwt,iwrk)
      ifound = max(ifound,iwrk)
      do i = 1,4
        idmat(i,2) = idwt(i)
      enddo
      idmat(5,2) = iwrk
      
      return
      end

C     =========================================
      subroutine zmSetupL0(istf,ityp,idwt,iwrk)
C     =========================================

C--   Setup weight indices and pdf mask
C--
C--   istf     (in)  :  1,2,3,4 = FL, F2, xF3, FL'
C--   ityp     (in)  :  1-6 = GL, QQ, SI, NS, NP, VA+NM 
C--   idwt(4)  (out) :  id(1)-id(3)  weight table id for LO, NLO, 
C--                                  NNLO (0=no table)
C--                     id(4)        power of alphas at LO (0 or 1)
C--   iwrk     (out) :  0 = no tables for combination istf, ityp

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dimension idwt(4)
      dimension itable(5,6,4)
C--           i =   LO  NLO NNLO   as iwrk
      data itable /  0, 202, 210,   0,   1, !1=GL      k=1 FL
     +               0, 103,   0,   0,   1, !2=QQ 
     +               0,   0, 209,   0,   1, !3=SI 
     +               0,   0,   0,   0,   0, !4=NS 
     +               0,   0, 207,   0,   1, !5=NP 
     +               0,   0, 208,   0,   1, !6=VA+NM

     +               0, 201, 206,   0,   1, !1=GL      k=2 F2
     +             101, 102,   0,   0,   1, !2=QQ 
     +               0,   0, 205,   0,   1, !3=SI 
     +               0,   0,   0,   0,   0, !4=NS
     +               0,   0, 203,   0,   1, !5=NP 
     +               0,   0, 204,   0,   1, !6=VA+NM      

     +               0,   0,   0,   0,   0, !1=Gl      k=3 xF3
     +               0,   0,   0,   0,   0, !2=QQ
     +               0,   0,   0,   0,   0, !3=SI
     +             101, 104,   0,   0,   1, !4=NS
     +               0,   0, 211,   0,   1, !5=NP 
     +               0,   0, 212,   0,   1, !6=VA+NM     

     +             202, 210, 213,   1,   1, !1=GL      k=4 FL'
     +             103,   0,   0,   1,   1, !2=QQ 
     +               0, 209, 214,   1,   1, !3=SI
     +               0,   0, 215,   1,   1, !4=NS   
     +               0, 207,   0,   1,   1, !5=NP 
     +               0, 208,   0,   1,   1 /!6=VA+NM
     
C--  Remark on FL': in NNLO one should have a table for
C--                 NPlus (table 215 in fact) and another
C--                 one for NMin. The NMin coefficient function
C--                 is not available so we use the NPlus table
C--                 instead. Thus 215 is assigned to all NS

      do i = 1,4
        idwt(i) = itable(i,ityp,istf)
      enddo
      
      iwrk = itable(5,ityp,istf)
      
      return
      end
      
C     =========================================
      subroutine zmSetupL1(istf,ityp,idwt,iwrk)
C     =========================================

C--   Setup weight indices and pdf mask for L term
C--
C--   istf     (in)  :  1,2,3,4 = FL, F2, xF3, FL'
C--   ityp     (in)  :  1-6 = GL, QQ, SI, NS, NP, VA+NM 
C--   idwt(4)  (out) :  id(1)-id(3)  weight table id for LO, NLO, 
C--                                  NNLO (0=no table)
C--                     id(4)        power of alphas at LO (0 or 1)
C--   iwrk     (out) :  0 = no tables for combination istf, ityp


C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dimension idwt(4)
      dimension itable(5,6,4)
C--           i =   LO  NLO NNLO   as iwrk
      data itable /  0,   0, 226,   0,   1, !1=GL      k=1 FL
     +               0,   0,   0,   0,   0, !2=QQ
     +               0,   0, 225,   0,   1, !3=SI
     +               0,   0, 227,   0,   1, !4=NS
     +               0,   0,   0,   0,   0, !5=NP
     +               0,   0,   0,   0,   0, !6=VA+NM
     
     +               0, 217, 219,   0,   1, !1=GL      k=2 F2
     +               0, 216,   0,   0,   1, !2=QQ
     +               0,   0, 218,   0,   1, !3=SI
     +               0,   0,   0,   0,   0, !4=NS
     +               0,   0, 220,   0,   1, !5=NP
     +               0,   0, 221,   0,   1, !6=VA+NM
 
     +               0,   0,   0,   0,   0, !1=Gl      k=3 xF3
     +               0,   0,   0,   0,   0, !2=QQ
     +               0,   0,   0,   0,   0, !3=SI
     +               0, 216,   0,   0,   1, !4=NS
     +               0,   0, 228,   0,   1, !5=NP
     +               0,   0, 229,   0,   1, !6=VA+NM

     +               0, 226, 232,   1,   1, !1=GL      k=4 FL'
     +               0,   0,   0,   1,   0, !2=QQ
     +               0, 225, 231,   1,   1, !3=SI
     +               0, 227,   0,   1,   1, !4=NS
     +               0,   0, 233,   1,   1, !5=NP
     +               0,   0, 234,   1,   1 /!6=VA+NM


      do i = 1,4
        idwt(i) = itable(i,ityp,istf)
      enddo
      
      iwrk = itable(5,ityp,istf)
      
      return
      end
      
C     =========================================
      subroutine zmSetupL2(istf,ityp,idwt,iwrk)
C     =========================================

C--   Setup weight indices and pdf mask for L2 term
C--
C--   istf     (in)  :  1,2,3,4 = FL, F2, xF3, FL'
C--   ityp     (in)  :  1-6 = GL, QQ, SI, NS, NP, VA+NM 
C--   idwt(4)  (out) :  id(1)-id(3)  weight table id for LO, NLO, 
C--                                  NNLO (0=no table)
C--                     id(4)        power of alphas at LO (0 or 1)
C--   iwrk     (out) :  0 = no tables for combination istf, ityp

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dimension idwt(4)
      dimension itable(5,6,4)
C--           i =   LO  NLO NNLO   as iwrk
      data itable /  0,   0,   0,   0,   0, !1=GL      k=1 FL
     +               0,   0,   0,   0,   0, !2=QQ
     +               0,   0,   0,   0,   0, !3=SI
     +               0,   0,   0,   0,   0, !4=NS
     +               0,   0,   0,   0,   0, !5=NP
     +               0,   0,   0,   0,   0, !6=VA+NM
     
     +               0,   0, 223,   0,   1, !1=GL      k=2 F2
     +               0,   0,   0,   0,   0, !2=QQ
     +               0,   0, 222,   0,   1, !3=SI
     +               0,   0, 224,   0,   1, !4=NS
     +               0,   0,   0,   0,   0, !5=NP
     +               0,   0,   0,   0,   0, !6=VA+NM

     +               0,   0,   0,   0,   0, !1=GL      k=3 xF3
     +               0,   0,   0,   0,   0, !2=QQ
     +               0,   0,   0,   0,   0, !3=SI
     +               0,   0, 230,   0,   1, !4=NS
     +               0,   0,   0,   0,   0, !5=NP
     +               0,   0,   0,   0,   0, !6=VA+NM

     +               0,   0, 236,   1,   1, !1=GL      k=4 FL'
     +               0,   0,   0,   1,   0, !2=QQ
     +               0,   0, 235,   1,   1, !3=SI
     +               0,   0,   0,   1,   0, !4=NS
     +               0,   0, 237,   1,   1, !5=NP
     +               0,   0, 238,   1,   1 /!6=VA+NM
 
      do i = 1,4
        idwt(i) = itable(i,ityp,istf)
      enddo
      
      iwrk = itable(5,ityp,istf)
      
      return
      end
      
C=====================================================================

C--   Below is some extra code to illustrate how to prototype a 
C--   structure function calculation using the QCDNUM slow convolution
C--   engine. The resulting code is straight forward but much slower
C--   than the fast code given above.
C--
C--   The weight table identifiers in common 'zmwidee.inc' are set-up
C--   by the weight filling routines which can be found in zmweits.f
C--
C--   Q2 scale variation w.r.t. mu2 is not implemented here.

C     =========================================
      subroutine zmslowf(istf,def,x,q,f,n,isel)
C     =========================================

C--   Slow structure function calculation
C--
C--   isel selects the terms to be included
C--
C--   101,102,103 -->  LO,NLO,NNLO quark+gluon term
C--   201,202,203 -->  LO,NLO,NNLO quark term only
C--   301,302,303 -->  LO,NLO,NNLO gluon term only
C--   other       -->  Full structure function (all terms)
C--
C--   When isel >  0  do grid boundary check (set ichk = 1)
C--   When isel <= 0  no grid boundary check (set ichk = 0)

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'

      external dzmFLij,dzmF2ij,dzmF3ij,dzmFpij
      
      logical lquark, lgluon
      common /qgflags/ lquark(3),lgluon(3)
      common /pdfdefs/ qvec(-6:6)

      dimension def(-6:6),x(*),q(*),f(*)

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('ZMSLOWF')
      
      if(izini.ne.12345) stop
     + 'ZMSLOWF: please first call ZMFILLW or ZMREADW --> STOP'

C--   Force evolution parameters of set izpdf
      call IdScope(0.D0,izpdf)

C--   Copy to common block
      do i = -6,6
        qvec(i) = def(i)
      enddo
C--   Select what to calculate
      call zselect(isel,lquark,lgluon,ichk)
C--   Go for all points
      if(istf.eq.1)     then
        call stfunxq(dzmFLij,x,q,f,n,ichk)
      elseif(istf.eq.2) then
        call stfunxq(dzmF2ij,x,q,f,n,ichk)
      elseif(istf.eq.3) then
        call stfunxq(dzmF3ij,x,q,f,n,ichk)
      elseif(istf.eq.4) then
        call stfunxq(dzmFpij,x,q,f,n,ichk)
      else
        stop 'ZMSLOWF: invalid structure function label --> STOP'
      endif

      call clrUmsg

      return
      end
      
C     ===========================================
      subroutine zselect(isel,lquark,lgluon,ichk)
C     ===========================================

C--   Select quark and/or gluon terms according to isel
C--
C--   101,102,103 -->  LO,NLO,NNLO quark+gluon term
C--   201,202,203 -->  LO,NLO,NNLO quark term only
C--   301,302,303 -->  LO,NLO,NNLO gluon term only
C--   other       -->  Full structure function (all terms)
C--
C--   When isel >  0  do grid boundary check (set ichk = 1)
C--   When isel <= 0  no grid boundary check (set ichk = 0)

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      logical lquark(3), lgluon(3), lset
      
      if(isel.le.0) then
        ichk = 0          !no grid boundary check
      else
        ichk = 1          !yes grid boundary check
      endif
      
      jsel = abs(isel)
      lset = .false.
      
      do iord = 1,3
        lquark(iord) = .false.
        lgluon(iord) = .false.
        if    (jsel.eq.100+iord) then
          lquark(iord) = .true.
          lgluon(iord) = .true.
          lset         = .true.
        elseif(jsel.eq.200+iord) then
          lquark(iord) = .true.
          lset         = .true.
        elseif(jsel.eq.300+iord) then
          lgluon(iord) = .true.
          lset         = .true.
        endif        
      enddo
      if(.not.lset) then
        do iord = 1,3
          lquark(iord) = .true.
          lgluon(iord) = .true.
        enddo
      endif    
*mb      write(6,*) ' '
*mb      write(6,*) 'lquark =',lquark,jsel
*mb      write(6,*) 'lgluon =',lgluon,jsel
      
      return
      end

C     ========================================
      double precision function dzmF2ij(ix,iq)
C     ========================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmwidee.inc'

C--   If(lquark(iord)) add quark term
C--   If(lqluon(iord)) add gluon term
      logical lquark, lgluon
      common /qgflags/ lquark(3),lgluon(3)
      common /pdfdefs/ qvec(-6:6)

      dimension evec(12),evpars(13)
      
      jset = izpdf
      
C--   Default if things go wrong
      call getval('null',qnull)
C--   Tolerance in floating point comparison
      call getval('epsi',epsi)
C--   Get evolution parameters
      call cpypar(evpars,13,jset)
C--   QCD order
      iord = int(evpars(1))
C--   Number of flavors
      nfl = nfrmiq(jset,iq,ithrs)
C--   Get linear combination
      call efromqq(qvec,evec,nfl)
C--   Global identifier (jsetw in common block zmstore)
      ig0 = 1000*jsetw

C--   LO
      term = 0.D0
      if(lquark(1)) then
        do id = 1,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       Singlet and ns+
          if(abs(evec(id)).gt.epsi) then
*            df2  = FcrossK(zmstor,ig0+idwtLO,jset,id1,ix,iq)
            df2  = BVALIJ(jset,id,ix,iq,-1)
            term = term + evec(id) * df2
          endif
C--       Valence and ns-
          if(abs(evec(id+6)).gt.epsi) then
*            df2  = FcrossK(zmstor,ig0+idwtLO,jset,id6,ix,iq)
            df2  = BVALIJ(jset,id+6,ix,iq,-1)
            term = term + evec(id+6) * df2
          endif
        enddo
      endif
C--   End of LO
      dzmF2ij = term
      if(iord.le.1) return

C--   NLO
      term = 0.D0
C--   Gluon
      idg = ipdftab(jset,0)
      if(lgluon(2).and.abs(evec(1)).gt.epsi) then
        df2  = FcrossK(zmstor,ig0+idC2G1,jset,idg,ix,iq)
        term = term + evec(1) * df2
      endif
C--   Quarks
      if(lquark(2)) then
        do id = 1,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       Singlet and ns+
          if(abs(evec(id)).gt.epsi) then
            df2  = FcrossK(zmstor,ig0+idC2Q1,jset,id1,ix,iq)
            term = term + evec(id) * df2
          endif
C--       Valence and ns-
          if(abs(evec(id+6)).gt.epsi) then
            df2  = FcrossK(zmstor,ig0+idC2Q1,jset,id6,ix,iq)
            term = term + evec(id+6) * df2
          endif
        enddo
      endif
C--   End of NLO
      dzmF2ij = dzmF2ij + altabn(jset,iq,-1,ierr) * term
      if(iord.le.2) return

C--   NNLO
      term = 0.D0
C--   Gluon and singlet
      if(abs(evec(1)).gt.epsi) then
        idg = ipdftab(jset,0)
        ids = ipdftab(jset,1)
        df2g  = 0.D0
        if(lgluon(3)) df2g  = FcrossK(zmstor,ig0+idC2G2,jset,idg,ix,iq)
        df2s  = 0.D0
        if(lquark(3)) df2s  = FcrossK(zmstor,ig0+idC2S2,jset,ids,ix,iq)
        term  = term + evec(1) * (df2g+df2s)
      endif
      if(lquark(3)) then
C--     Valence
        if(abs(evec(7)).gt.epsi) then
          id7   = ipdftab(jset,7)
          df2v  = FcrossK(zmstor,ig0+idC2M2,jset,id7,ix,iq)
          term  = term + evec(7) * df2v
        endif
C--     Nonsinglet
        do id = 2,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       NS+
          if(abs(evec(id)).gt.epsi) then
            df2  = FcrossK(zmstor,ig0+idC2P2,jset,id1,ix,iq)
            term = term + evec(id) * df2
          endif
C--       NS-
          if(abs(evec(id+6)).gt.epsi) then
            df2  = FcrossK(zmstor,ig0+idC2M2,jset,id6,ix,iq)
            term = term + evec(id+6) * df2
          endif
        enddo
      endif  
C--   End of NNLO
      dzmF2ij = dzmF2ij + altabn(jset,iq,-2,ierr) * term

      return
      end

C     ========================================
      double precision function dzmFLij(ix,iq)
C     ========================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmwidee.inc'

C--   If(lquark(iord)) add quark term
C--   If(lqluon(iord)) add gluon term
      logical lquark, lgluon
      common /qgflags/ lquark(3),lgluon(3)
      common /pdfdefs/ qvec(-6:6)

      dimension evec(12),evpars(13)
      
      jset = izpdf

C--   Default if things go wrong
      call getval('null',qnull)
C--   Tolerance in floating point comparison
      call getval('epsi',epsi)
C--   Get evolution parameters
      call cpypar(evpars,13,jset)
C--   QCD order
      iord = int(evpars(1))
C--   Number of flavors
      nfl = nfrmiq(jset,iq,ithrs)
C--   Get linear combination
      call efromqq(qvec,evec,nfl)
C--   Global identifier (jsetw in common block zmstore)
      ig0 = 1000*jsetw

C--   LO
      dzmFLij = 0.D0
      if(iord.le.1) return

C--   NLO
      term = 0.D0
C--   Gluon
      if(lgluon(2).and.abs(evec(1)).gt.epsi) then
        idg  = ipdftab(jset,0)
        dfL  = FcrossK(zmstor,ig0+idCLG1,jset,idg,ix,iq)
        term = term + evec(1) * dfL
      endif
C--   Quarks
      if(lquark(2)) then
        do id = 1,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       Singlet and ns+
          if(abs(evec(id)).gt.epsi) then
            dfL  = FcrossK(zmstor,ig0+idCLQ1,jset,id1,ix,iq)
            term = term + evec(id) * dfL
          endif
C--       Valence and ns-
          if(abs(evec(id+6)).gt.epsi) then
            dfL  = FcrossK(zmstor,ig0+idCLQ1,jset,id6,ix,iq)
            term = term + evec(id+6) * dfL
          endif
        enddo
      endif  
C--   End of NLO
      dzmFLij = dzmFLij + altabn(jset,iq,-1,ierr) * term
      if(iord.le.2) return

C--   NNLO
      term = 0.D0
C--   Gluon and singlet
      if(abs(evec(1)).gt.epsi) then
        idg = ipdftab(jset,0)
        ids = ipdftab(jset,1)
        dfLg  = 0.D0
        if(lgluon(3)) dfLg  = FcrossK(zmstor,ig0+idCLG2,jset,idg,ix,iq)
        dfLs  = 0.D0
        if(lquark(3)) dfLs  = FcrossK(zmstor,ig0+idCLS2,jset,ids,ix,iq)
        term  = term + evec(1) * (dfLg+dfLs)
      endif
      if(lquark(3)) then
C--     Valence
        if(abs(evec(7)).gt.epsi) then
          id7   = ipdftab(jset,7)
          dfLv  = FcrossK(zmstor,ig0+idCLM2,jset,id7,ix,iq)
          term  = term + evec(7) * dfLv
        endif
C--     Nonsinglet
        do id = 2,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       NS+
          if(abs(evec(id)).gt.epsi) then
            dfL  = FcrossK(zmstor,ig0+idCLP2,jset,id1,ix,iq)
            term = term + evec(id) * dfL
          endif
C--       NS-
          if(abs(evec(id+6)).gt.epsi) then
            dfL  = FcrossK(zmstor,ig0+idCLM2,jset,id6,ix,iq)
            term = term + evec(id+6) * dfL
          endif
        enddo
      endif  
C--   End of NNLO
      dzmFLij = dzmFLij + altabn(jset,iq,-2,ierr) * term

      return
      end

C     ========================================
      double precision function dzmF3ij(ix,iq)
C     ========================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmwidee.inc'

C--   If(lquark(iord)) add quark term
C--   If(lqluon(iord)) add gluon term
      logical lquark, lgluon
      common /qgflags/ lquark(3),lgluon(3)
      common /pdfdefs/ qvec(-6:6)

      dimension evec(12),evpars(13)
      
      jset = izpdf

C--   Default if things go wrong
      call getval('null',qnull)
C--   Tolerance in floating point comparison
      call getval('epsi',epsi)
C--   Get evolution parameters
      call cpypar(evpars,13,jset)
C--   QCD order
      iord = int(evpars(1))
C--   Number of flavors
      nfl = nfrmiq(jset,iq,ithrs)
C--   Get linear combination
      call efromqq(qvec,evec,nfl)
C--   Global identifier (jsetw in common block zmstore)
      ig0 = 1000*jsetw

C--   LO
      term = 0.D0
      if(lquark(1)) then
C--     Valence
        if(abs(evec(7)).gt.epsi) then
          id7  = ipdftab(jset,7)
          df3  = FcrossK(zmstor,ig0+idwtLO,jset,id7,ix,iq)
          term = term + evec(7) * df3
        endif
C--     Nonsinglet
        do id = 2,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       NS+
          if(abs(evec(id)).gt.epsi) then
            df3  = FcrossK(zmstor,ig0+idwtLO,jset,id1,ix,iq)
            term = term + evec(id) * df3
          endif
C--       NS-
          if(abs(evec(id+6)).gt.epsi) then
            df3  = FcrossK(zmstor,ig0+idwtLO,jset,id6,ix,iq)
            term = term + evec(id+6) * df3
          endif
        enddo
      endif  
C--   End of LO
      dzmF3ij = term
      if(iord.le.1) return

C--   NLO
      term = 0.D0
      if(lquark(2)) then
C--     Valence
        if(abs(evec(7)).gt.epsi) then
          id7  = ipdftab(jset,7)
          df3  = FcrossK(zmstor,ig0+idC3Q1,jset,id7,ix,iq)
          term = term + evec(7) * df3
        endif
C--     Nonsinglet
        do id = 2,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       NS+
          if(abs(evec(id)).gt.epsi) then
            df3  = FcrossK(zmstor,ig0+idC3Q1,jset,id1,ix,iq)
            term = term + evec(id) * df3
          endif
C--       NS-
          if(abs(evec(id+6)).gt.epsi) then
            df3  = FcrossK(zmstor,ig0+idC3Q1,jset,id6,ix,iq)
            term = term + evec(id+6) * df3
          endif
        enddo
      endif  
C--   End of NLO
      dzmF3ij = dzmF3ij + altabn(jset,iq,-1,ierr) * term
      if(iord.le.2) return

C--   NNLO
      term = 0.D0
      if(lquark(3)) then
C--     Valence
        if(abs(evec(7)).gt.epsi) then
          id7   = ipdftab(jset,7)
          df3v  = FcrossK(zmstor,ig0+idC3M2,jset,id7,ix,iq)
          term  = term + evec(7) * df3v
        endif
C--     Nonsinglet
        do id = 2,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       NS+
          if(abs(evec(id)).gt.epsi) then
            df3  = FcrossK(zmstor,ig0+idC3P2,jset,id1,ix,iq)
            term = term + evec(id) * df3
          endif
C--       NS-
          if(abs(evec(id+6)).gt.epsi) then
            df3  = FcrossK(zmstor,ig0+idC3M2,jset,id6,ix,iq)
            term = term + evec(id+6) * df3
          endif
        enddo
      endif  
C--   End of NNLO
      dzmF3ij = dzmF3ij + altabn(jset,iq,-2,ierr) * term

      return
      end
      
C     ========================================
      double precision function dzmFpij(ix,iq)
C     ========================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmwidee.inc'

C--   If(lquark(iord)) add quark term
C--   If(lqluon(iord)) add gluon term
      logical lquark, lgluon
      common /qgflags/ lquark(3),lgluon(3)
      common /pdfdefs/ qvec(-6:6)

      dimension evec(12),evpars(13)
      
      jset = izpdf
      
C--   Default if things go wrong
      call getval('null',qnull)
C--   Tolerance in floating point comparison
      call getval('epsi',epsi)
C--   Get evolution parameters
      call cpypar(evpars,13,jset)
C--   QCD order
      iord = int(evpars(1))
C--   Number of flavors
      nfl = nfrmiq(jset,iq,ithrs)
C--   Get linear combination
      call efromqq(qvec,evec,nfl)
C--   Global identifier (jsetw in common block zmstore)
      ig0 = 1000*jsetw

C--   LO
      term = 0.D0
C--   Gluon
      if(lgluon(1).and.abs(evec(1)).gt.epsi) then
        idg  = ipdftab(jset,0)
        dfp  = FcrossK(zmstor,ig0+idCLG1,jset,idg,ix,iq)
        term = term + evec(1) * dfp
      endif
C--   Quarks
      if(lquark(1)) then
        do id = 1,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       Singlet and ns+
          if(abs(evec(id)).gt.epsi) then
            dfp  = FcrossK(zmstor,ig0+idCLQ1,jset,id1,ix,iq)
            term = term + evec(id) * dfp
          endif
C--       Valence and ns-
          if(abs(evec(id+6)).gt.epsi) then
            dfp  = FcrossK(zmstor,ig0+idCLQ1,jset,id6,ix,iq)
            term = term + evec(id+6) * dfp
          endif
        enddo
      endif  
C--   End of LO
      dzmFpij = altabn(jset,iq,1,ierr) * term
      if(iord.le.1) return

C--   NLO
      term = 0.D0
C--   Gluon and singlet
      if(abs(evec(1)).gt.epsi) then
        idg   = ipdftab(jset,0)
        ids   = ipdftab(jset,1)
        dfpg  = 0.D0
        if(lgluon(2)) dfpg  = FcrossK(zmstor,ig0+idCLG2,jset,idg,ix,iq)
        dfps  = 0.D0
        if(lquark(2)) dfps  = FcrossK(zmstor,ig0+idCLS2,jset,ids,ix,iq)
        term  = term + evec(1) * (dfpg+dfps)
      endif
      if(lquark(2)) then
C--     Valence
        if(abs(evec(7)).gt.epsi) then
          id7   = ipdftab(jset,7)
          dfpv  = FcrossK(zmstor,ig0+idCLM2,jset,id7,ix,iq)
          term  = term + evec(7) * dfpv
        endif
C--     Nonsinglet
        do id = 2,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       NS+
          if(abs(evec(id)).gt.epsi) then
            dfp  = FcrossK(zmstor,ig0+idCLP2,jset,id1,ix,iq)
            term = term + evec(id) * dfp
          endif
C--       NS-
          if(abs(evec(id+6)).gt.epsi) then
            dfp  = FcrossK(zmstor,ig0+idCLM2,jset,id6,ix,iq)
            term = term + evec(id+6) * dfp
          endif
        enddo
      endif  
C--   End of NLO
      dzmFpij = dzmFpij + altabn(jset,iq,2,ierr) * term
      if(iord.le.2) return

C--   NNLO
      term = 0.D0
C--   Gluon and singlet
      if(abs(evec(1)).gt.epsi) then
        idg   = ipdftab(jset,0)
        ids   = ipdftab(jset,1)
        dfpg  = 0.D0
        if(lgluon(3)) dfpg  = FcrossK(zmstor,ig0+idCLG3,jset,idg,ix,iq)
        dfps  = 0.D0
        if(lquark(3)) dfps  = FcrossK(zmstor,ig0+idCLS3,jset,ids,ix,iq)
        term  = term + evec(1) * (dfpg+dfps)
      endif
      if(lquark(3)) then
C--     Valence
        if(abs(evec(7)).gt.epsi) then
          id7   = ipdftab(jset,7)
          dfpv  = FcrossK(zmstor,ig0+idCLN3,jset,id7,ix,iq)
          term  = term + evec(7) * dfpv
        endif
C--     Nonsinglet
        do id = 2,nfl
          id1 = ipdftab(jset,id)
          id6 = ipdftab(jset,id+6)
C--       NS+
          if(abs(evec(id)).gt.epsi) then
            dfp  = FcrossK(zmstor,ig0+idCLN3,jset,id1,ix,iq)
            term = term + evec(id) * dfp
          endif
C--       NS-
          if(abs(evec(id+6)).gt.epsi) then
            dfp  = FcrossK(zmstor,ig0+idCLN3,jset,id6,ix,iq)
            term = term + evec(id+6) * dfp
          endif
        enddo
      endif  
C--   End of NNLO
      dzmFpij = dzmFpij + altabn(jset,iq,3,ierr) * term

      return
      end
      
