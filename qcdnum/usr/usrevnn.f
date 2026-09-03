

C--   This is the file usrevnn.f containing the evolution toolbox routines

C--   subroutine EvFillA (w,id,func)
C--   double precision function EvGetAA(w,id,iq,nf,ithresh)
C--   subroutine EvDglap(w,subr,idw,ida,idf,n,iq0,eps)
C--   subroutine CpyParW(w,array,n,kset)
C--   subroutine UseParW(w,kset)
C--   integer function KeyParW(w,kset)
C--   integer function KeyGrpW(w,kset,igroup)
C--   double precision function EvPdfij(w,id,ix,iq,jchk)
C--   subroutine EvPlist(w,id,x,qmu2,pdf,m,jchk)
C--   subroutine EvTable(w,id,x,nx,q,nq,table,jchk)
C--   subroutine EvPCopy(w,id,def,n,jset)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM toolbox evolution from usrevnn.f                 */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM toolbox evolution from usrevnn.f                   */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM toolbox evolution from usrevnn.f                   */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    void evfilla(double *w, int id, double (*func)(int*,int*,int*));
C-----------------------------------------------------------------------
CXXHFW  #define fevfilla FC_FUNC(evfilla,EVFILLA)
CXXHFW    void fevfilla(double*,int*,double(*)(int*,int*,int*));
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void evfilla(double *w, int id, double (*func)(int*,int*,int*))
CXXWRP  {
CXXWRP    fevfilla(w, &id, func);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================
      subroutine EvFillA (w,id,func)
C     ==============================

C--   Fill a type-6 table with expansion coefficients
C--
C--   w     (in)   Store containing type-6 table
C--   id    (in)   Type-6 table identifier in global format
C--   func  (in)   function declared external in the calling routine

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     6   / , iflg /   0,    0  /
      data icma  /    1,     6   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam /'EVFILLA ( W, ID, FUNC )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check global identifier
      igl = iqcSjekId(subnam,'ID',w,id,icmi,icma,iflg,lint)
C--   Point to base set
      call sparParTo5(1)

C--   Do the work
      call sqcEvFillA(w,igl,func)

C--   Store parameter version
      ipver = int(dparGetPar(pars8,1,idipver8))
      ksetw = igl/1000
      call sparSetPar(w,ksetw,idipver8,dble(ipver))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double evgetaa(double *w, int id, int iq, int &nf, int &ithresh);
C-----------------------------------------------------------------------
CXXHFW  #define fevgetaa FC_FUNC(evgetaa,EVGETAA)
CXXHFW    double fevgetaa(double*,int*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double evgetaa(double *w, int id, int iq, int &nf, int &ithresh)
CXXWRP  {
CXXWRP    return fevgetaa(w, &id, &iq, &nf, &ithresh);
CXXWRP  }
C-----------------------------------------------------------------------

C     =====================================================
      double precision function EvGetAA(w,id,iq,nf,ithresh)
C     =====================================================

C--   Get expansion coefficient from a type-6 table
C--
C--   w        (in)   Store containing type-6 table
C--   id       (in)   Type-6 table identifier in global format
C--   iq       (in)   Qgrid point (may be < 0)
C--   nf      (out)   Number of flavours at iq
C--   ithresh (out)   Threshold indicator 0,+-1

C--   for iq < 0 nf = 3,4,5 instead of 4,5,6 at iqc,b,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     6   / , iflg /   0,    0  /
      data icma  /    1,     6   /

      logical lint

      dimension w(*)

      character*20 etxt
      character*80 emsg

      character*80 subnam
      data subnam /'EVGETAA ( W, ID, IQ, NF, ITHRESH )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check global identifier
      igl = iqcSjekId(subnam,'ID',w,id,icmi,icma,iflg,lint)

C--   Check iq
      call sqcIlele(subnam,'IQ',1,abs(iq),ntt2,' ')

C--   Point to proper parameter set
      ksetw = igl/1000
      ipver = int(dparGetPar(w,ksetw,idipver8))
      if(ipver.le.0) then
        call smb_itoch(igl,etxt,ltxt)
        write(emsg,
     +      '(''Table id = '',A,'' does not exist or is empty'')')
     +          etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif
      call sparParTo5(ipver)

C--   Do the work
      EvGetAA = sqcEvGetAA(w,igl,iq,nf,ithresh)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void evdglap(double *w, int idw, int ida, int idf, double *start, int m, int n, int *iqlim, int &nf, double &eps);
C-----------------------------------------------------------------------
CXXHFW  #define fevdglap FC_FUNC(evdglap,EVDGLAP)
CXXHFW    void fevdglap(double*,int*,int*,int*,double*,int*,int*,int*,int*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void evdglap(double *w, int idw, int ida, int idf, double *start, int m, int n, int *iqlim, int &nf, double &eps)
CXXWRP  {
CXXWRP    fevdglap(w, &idw, &ida, &idf, start, &m, &n, iqlim, &nf, &eps);
CXXWRP  }
C-----------------------------------------------------------------------

C     ========================================================
      subroutine EvDglap(w,idw,ida,idf,start,m,n,iqlim,nf,eps)
C     ========================================================

C--   n-fold evolution routine
C--
C--   w          (in)    store with weight, coefficient and pdf tables
C--   idw(i,j,k) (in)    Pij weight table identifiers of order k
C--   ida(i,j,k) (in)    aij coefficient table identifiers of order k
C--   idf(i)     (in)    Pdf identifiers
C--   start(i,j) (inout) either start pdf_i(xj) or discontinuity at iq0
C--   m          (in)    idw(m,m,.), ida(m,m,.) and start(m,.)
C--   n          (in)    abs(n) = number of pdfs to evolve simultaneously
C--                      n > 0  = select input mode
C--                      n < 0  = select transfer mode
C--   iqlim(2)   (inout) iq0,iq1 on entry, and iq0,iq1 on exit
C--   nf         (out)   number of flavours; nf < 0 flags end of range
C--   eps        (out)   max deviation at midpoint of xbins
C--
C--   n > 0 input    mode: start = input on entry, evolved pdf on exit
C--   n < 0 transfer mode: start = discontinuity on entry, zero on exit

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

      dimension icmiw(2), icmaw(2), iflgw(2)
C--                isign  itype              pzi  cfil
      data icmiw  /    1,     1   /, iflgw /   1,    1  /
      data icmaw  /    1,     4   /
      dimension icmia(2), icmaa(2), iflga(2)
C--                isign  itype              pzi  cfil
      data icmia  /    1,     6   /, iflga /   1,    1  /
      data icmaa  /    1,     6   /
      dimension icmif(2), icmaf(2), iflgf(2)
C--                isign  itype              pzi  cfil
      data icmif  /    1,     5   /, iflgf /   0,    0  /
      data icmaf  /    1,     5   /

      logical lint

      dimension w(*),idw(m,m,*),ida(m,m,*),idf(*)
      dimension start(m,*)
      dimension iqlim(2)

      character*80 subnam
      data subnam /
     +'EVDGLAP ( W, IDW, IDA, IDF, START, M, N, IQLIM, NF, EPS )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check number of coupled evolutions
      nabs = abs(n)
      call sqcIlele(subnam,'N',1,nabs,mce0,
     + 'To many coupled evolutions, please increase MCE0 in qcdnum.inc')

C--   Check dimension M
      if(m.lt.nabs) call sqcErrMsg(subnam,'Invalid dimension M < N' )

C--   Check order
      call sqcIlele(subnam,'IORD',1,iord6,nnopt6(0),
     + 'Please call SETORD(IORD) with correct value of IORD')

C--   Check Pij and Alfa table identifiers
      nwt = 0
      do i = 1,nabs
        do j = 1,nabs
          do k = 1,iord6
            ierr =
     +      iqcSjekId(subnam,'IDW',w,idw(i,j,k),icmiw,icmaw,iflgw,lint)
            ierr =
     +      iqcSjekId(subnam,'IDA',w,ida(i,j,k),icmia,icmaa,iflga,lint)
            if(idw(i,j,k).ne.0 .and. ida(i,j,k).eq.0) then
              call sqcErrMsg(subnam,
     +        'Found IDW(i,j,k) without associated IDA(i,j,k)')
            endif
            if(idw(i,j,k).ne.0) nwt = nwt+1
          enddo
        enddo
      enddo

      if(nwt.eq.0) call sqcErrMsg(subnam,'Array IDW is empty')

C--   Check Pdf table identifiers
      do i = 1,nabs
        ierr =
     +      iqcSjekId(subnam,'IDF',w,idf(i),icmif,icmaf,iflgf,lint)
      enddo

C--   Check that all pdfs are in the same set
      ksetf = iqcGetSetNumber(idf(1))
      do i = 1,nabs
        if(iqcGetSetNumber(idf(i)).ne.ksetf) then
          call sqcErrMsg(subnam,
     +        'Not all IDF identifiers are in the same table set')
        endif
      enddo

C--   Check that all alfas tables are in the same set
      klast = 0
      do i = 1,nabs
        do j = 1,nabs
          do k = 1,iord6
            if(ida(i,j,k).ne.0) then
              kseta   = iqcGetSetNumber(ida(i,j,k))
              if(klast.ne.0 .and. kseta.ne.klast) then
                call sqcErrMsg(subnam,
     +             'Not all IDA identifiers are in the same table set')
              endif
              klast = kseta
            endif
          enddo
        enddo
      enddo

C--   Check parameter version of alfas table
      ipvera = int(dparGetPar(w,kseta,idipver8))
      ipver0 = int(dparGetPar(pars8,1,idipver8))
      if(ipvera.ne.ipver0) call sqcErrMsg(subnam,
     +        'IDA tables not created with current parameters')

C--   Point to base set
      call sparParTo5(1)

C--   Set evolution range
      it1 = iqlim(1)
      it2 = iqlim(2)
C--   First check that starting point is within cuts
      if(it1.lt.itmic5 .or. it1.gt.itmac5) then
        nf = -1
        return
      endif
C--   Set end point and iz range forward evolution
      if(it2.ge.it1) then
C--     Forward evolution it2 >= it1
        iz1   = izfit5( it1)
        nf    = itfiz5(-iz1)
        itlim = itfiz5(iz25(nf))
        it2   = min(it2,itlim,itmac5)
        iz2   = izfit5(-it2)
      else
C--     Backward evolution it2 < it1 
        iz1   = izfit5(-it1)
        nf    = itfiz5(-iz1)
        itlim = itfiz5(iz15(nf))
        it2   = max(it2,itlim,itmic5)
        iz2   = izfit5(it2)
      endif

C--   Set start values
      do i = 1,nabs
        id = idf(i)
        do ig = 1,nyg2
          ia = iqcG5ijk(w,1,-ig,id)-1
          if(n.gt.0) then
            do iy = 1,nyy2(0)
              ix    = nyy2(0)+1-iy
              ia    = ia+1
              w(ia) = start(i,ix) !input mode: set start
            enddo
          else
            do iy = 1,nyy2(0)
              ix    = nyy2(0)+1-iy
              ia    = ia+1
              w(ia) = w(ia)+start(i,ix) !transfer mode: add discontinuity
            enddo
          endif
        enddo
      enddo

C--   Number of perturbative terms
      nopt = nnopt6(iord6)

C--   Do the work
      call sqcNNallg(w,idw,w,ida,w,idf,m,nopt,nf,iz1,iz2,nabs,eps)

C--   Pass end values
      do i = 1,nabs
        id = idf(i)
        ia = iqcG5ijk(w,1,iz2,id)-1
        if(n.gt.0) then
          do iy = 1,nyy2(0)
            ix          = nyy2(0)+1-iy
            ia          = ia+1
            start(i,ix) = w(ia) !input mode: return endvalue
          enddo
        else
          do iy = 1,nyy2(0)
            ix          = nyy2(0)+1-iy
            start(i,ix) = 0.D0  !transfer mode: return zero
          enddo
        endif
      enddo

C--   Pass evolution range
      iqlim(1) = itfiz5(iz1)
      iqlim(2) = itfiz5(iz2)

C--   Flag end of range
      if(iqlim(2).eq.itmic5 .or. iqlim(2).eq.itmac5) nf = -nf

C--   Check max deviation
      if(dflim6.gt.0.D0 .and. eps.gt.dflim6) call sqcErrMsg(subnam, 
     +          'Possible spline oscillation detected')

C--1  Find current key
      ksetw   = iqcGetSetNumber(idf(1))
      keypdf  = int(dparGetPar(w,ksetw,idipver8))
      if(keypdf.lt.0 .or. keypdf.gt.mpl0) stop 'EVDGLAP: invalid key'
C--2  Get base_key
      keybase = iparGetGroupKey(pars8,1,6)
C--3  Do nothing if pdf_key = base_key
      if(keypdf.ne.keybase) then
        call sparCountDn(keypdf)
        call sparBaseToKey(keybase)
        call sparCountUp(keybase)
        call sparParAtoB(pars8,keybase,w,ksetw)
      endif

C--   Store the type of evolution (type = 4)
      call sparSetPar(w,ksetw,idievtyp8,4.D0)
      call sparSetPar(w,ksetw,idnfheavy8,6.D0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void cpyparw(double *w, double *array, int n, int kset);
C-----------------------------------------------------------------------
CXXHFW  #define fcpyparw FC_FUNC(cpyparw,CPYPARW)
CXXHFW    void fcpyparw(double*,double*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void cpyparw(double *w, double *array, int n, int kset)
CXXWRP  {
CXXWRP    fcpyparw(w, array, &n, &kset);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================
      subroutine CpyParW(w,array,n,kset)
C     ==================================

C--   Copy parameter list to a local array

C--   w      (in) : workspace
C--   array (out) : array with parameter values
C--   n      (in) : dimension of pars declared in the calling routine
C--   kset   (in) : table set identifier in w

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'CPYPARW ( W, ARRAY, N, ISET )'/

      dimension w(*), array(*)

      logical lqcIsetExists

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check n in range
      call sqcIlele(subnam,'N',13,n,9999,' ')

C--   Get the key
      key    = 0                                 !avoid compiler warning
      ievtyp = 0                                 !avoid compiler warning
      if(w(1).eq.0.D0) then
C--     Internal memory
        call sqcIlele(subnam,'ISET',0,kset,mset0,' ')
        if(kset.eq.0) then
          key    = int(dparGetPar(pars8,1,idipver8))
          ievtyp = 0
        elseif(isetf7(kset).ne.0) then
          key    = int(dparGetPar(stor7,isetf7(kset),idipver8))
          ievtyp = int(dparGetPar(stor7,isetf7(kset),idievtyp8))
        else
          call sqcSetMsg(subnam,'ISET',kset)
        endif
      elseif(lqcIsetExists(w,kset)) then
C--     Toolbox workspace
        key    = int(dparGetPar(w,kset,idipver8))
        ievtyp = int(dparGetPar(w,kset,idievtyp8))
      else
C--     Pdf set does not exist
        call sqcSetMsg(subnam,'ISET',kset)
      endif

C--   Copy params to array
      call sparListPar(key,array,ierr)
      array(13) = ievtyp

      if(ierr.ne.0) then
C--     Kset does not exist or has no parameters
        call sqcSetMsg(subnam,'ISET',kset)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void useparw(double *w, int kset);
C-----------------------------------------------------------------------
CXXHFW  #define fuseparw FC_FUNC(useparw,USEPARW)
CXXHFW    void fuseparw(double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void useparw(double *w, int kset)
CXXWRP  {
CXXWRP    fuseparw(w, &kset);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================
      subroutine UseParW(w,kset)
C     ==========================

C--   Copy parameters of kset back into qpars6

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'USEPARW ( W, ISET )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Get the key
      key = 0                                    !avoid compiler warning
      if(w(1).eq.0.D0) then
C--     Internal memory
        call sqcIlele(subnam,'ISET',0,kset,mset0,' ')
        if(kset.eq.0) then
          return                                 !base set nothing to do
        elseif(isetf7(kset).ne.0) then
          key = int(dparGetPar(stor7,isetf7(kset),idipver8))
        else
          call sqcSetMsg(subnam,'ISET',kset)
        endif
      elseif(lqcIsetExists(w,kset)) then
C--     Toolbox workspace
        key = int(dparGetPar(w,kset,idipver8))
      else
C--     Pdf set does not exist
        call sqcSetMsg(subnam,'ISET',kset)
      endif

C--   Pdf set not filled
      if(key.eq.0) call sqcSetMsg(subnam,'ISET',kset)

C--   Reset /qpars6/ and remake base
      call sparRemakeBase(key)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int keyparw(double *w, int kset);
C-----------------------------------------------------------------------
CXXHFW  #define fkeyparw FC_FUNC(keyparw,KEYPARW)
CXXHFW    int fkeyparw(double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int keyparw(double *w, int kset)
CXXWRP  {
CXXWRP    return fkeyparw(w, &kset);
CXXWRP  }
C-----------------------------------------------------------------------

C     ================================
      integer function KeyParW(w,kset)
C     ================================

C--   Returns the version number of a parameter set

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'KEYPARW ( W, ISET )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Get the key
      key = 0                                    !avoid compiler warning
      if(w(1).eq.0.D0) then
C--     Internal memory
        call sqcIlele(subnam,'ISET',0,kset,mset0,' ')
        if(kset.eq.0)  then
          key = int(dparGetPar(pars8,1,idipver8))
        elseif(isetf7(kset).ne.0) then
          key = int(dparGetPar(stor7,isetf7(kset),idipver8))
        else
          call sqcSetMsg(subnam,'ISET',kset)
        endif
      elseif(lqcIsetExists(w,kset)) then
C--     Toolbox workspace
        key = int(dparGetPar(w,kset,idipver8))
      else
C--     Pdf set does not exist
        call sqcSetMsg(subnam,'ISET',kset)
      endif

C--   Pdf set not filled
      if(key.eq.0) call sqcSetMsg(subnam,'ISET',kset)

      KeyParW = key

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int keygrpw(double *w, int kset, int igroup);
C-----------------------------------------------------------------------
CXXHFW  #define fkeygrpw FC_FUNC(keygrpw,KEYGRPW)
CXXHFW    int fkeygrpw(double*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int keygrpw(double *w, int kset, int igroup)
CXXWRP  {
CXXWRP    return fkeygrpw(w, &kset, &igroup);
CXXWRP  }
C-----------------------------------------------------------------------


C     =======================================
      integer function KeyGrpW(w,kset,igroup)
C     =======================================

C--   Returns the version number of a parameter set

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'KEYPARW ( W, ISET )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Get the key
      key = 0                                    !avoid compiler warning
      kgr = 0                                    !avoid compiler warning
      if(w(1).eq.0.D0) then
C--     Internal memory
        call sqcIlele(subnam,'ISET',0,kset,mset0,' ')
        if(kset.eq.0)  then
          key  = int(dparGetPar(pars8,1,idipver8))
          kgr  = iparGetGroupKey(pars8,1,igroup)
        elseif(isetf7(kset).ne.0) then
          key = int(dparGetPar(stor7,isetf7(kset),idipver8))
          kgr = iparGetGroupKey(pars8,key,igroup)
        else
          call sqcSetMsg(subnam,'ISET',kset)
        endif
      elseif(lqcIsetExists(w,kset)) then
C--     Toolbox workspace
        key = int(dparGetPar(w,kset,idipver8))
        kgr = iparGetGroupKey(pars8,key,igroup)
      else
C--     Pdf set does not exist
        call sqcSetMsg(subnam,'ISET',kset)
      endif

C--   Pdf set not filled
      if(key.eq.0) call sqcSetMsg(subnam,'ISET',kset)

      KeyGrpW = kgr

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double evpdfij(double *w, int id, int ix, int jq, int jchk);
C-----------------------------------------------------------------------
CXXHFW  #define fevpdfij FC_FUNC(evpdfij,EVPDFIJ)
CXXHFW    double fevpdfij(double*,int*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double evpdfij(double *w, int id, int ix, int jq, int jchk)
CXXWRP  {
CXXWRP    return fevpdfij(w, &id, &ix, &jq, &jchk);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================================
      double precision function EvPdfij(w,id,ix,jq,jchk)
C     ==================================================

C--   Get pdf at gridpoint
C--
C--   w      (in)  store (not aligned)
C--   id     (in)  pdf table identifier in global format
C--   ix     (in)  x-grid index
C--   jq     (in)  q-grid index can be < 0
C--   jchk   (in)  0 =   check input return null outside grid
C--                1 =   check input fatal error outside grid
C--               -1 = nocheck input fatal error outside grid

C--   for iq < 0 nf = 3,4,5 instead of 4,5,6 at iqc,b,t

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

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /   -1,     5   / , iflg /   0,    1  /
      data icma  /    1,     5   /

      logical lint

      dimension w(*)

      character*80 subnam
      data subnam /'EVPDFIJ ( W, ID, IX, IQ, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      iq = abs(jq)

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

      if(jchk.eq.-1 .and. .not.Lscopechek6) then
C--     Do not check pdf id
        igl = abs(id)
      else
C--     Check pdf id (also that pdf is filled with current params)
        igl = iqcSjekId(subnam,'ID',w,id,icmi,icma,iflg,lint)
      endif
*      if(lint) stop 'EVPDFIJ: cannot access pdfs in internal memory'

C--   Catch x = 1
      if(ix.eq.nyy2(0)+1) then
        EvPdfij = 0.D0
        return
      endif

C--   Point to the correct set
      ksetw  = abs(igl)/1000
      if(lint) then
        ipver = int(dparGetPar(stor7,ksetw,idipver8))
      else
        ipver = int(dparGetPar(  w  ,ksetw,idipver8))
      endif
      call sparParTo5(ipver)

C--   Ranges (taking cuts into account)
      ixmi  = nyy2(0)+1-iymac5
      iqmi  = itmic5
      iqma  = itmac5

C--   Check or nocheck thas the question
      if(jchk.ne.0) then
        call sqcIlele(subnam,'IX',ixmi,ix,nyy2(0),' ')
        call sqcIlele(subnam,'IQ',iqmi,iq,iqma   ,' ')
      else
        EvPdfij = qnull6
        if(ix.lt.ixmi .or. ix.gt.nyy2(0))    return
        if(iq.lt.iqmi .or. iq.gt.iqma)       return
      endif

C--   Do the work
      iy      = nyy2(0)+1-ix
      it      = jq
      if(Lint) then
        EvPdfij = dqcEvPdfij(stor7,igl,iy,it)
      else
        EvPdfij = dqcEvPdfij(  w  ,igl,iy,it)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void evplist(double *w, int id, double *x, double *q, double *f, int n, int jchk);
C-----------------------------------------------------------------------
CXXHFW  #define fevplist FC_FUNC(evplist,EVPLIST)
CXXHFW    void fevplist(double*,int*,double*,double*,double*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void evplist(double *w, int id, double *x, double *q, double *f, int n, int jchk)
CXXWRP  {
CXXWRP    fevplist(w, &id, x, q, f, &n, &jchk);
CXXWRP  }
C-----------------------------------------------------------------------

C     =====================================
      subroutine EvPlist(w,id,x,q,f,n,jchk)
C     =====================================

C--   Process an arbitrarily long list of interpolation points
C--
C--   w          (in)  Workspace
C--   id         (in)  Global identifier of the pdf to be interpolated
C--   x(n),q(n)  (in)  List of interpolation points
C--   f(n)      (out)  List of interpolated results
C--   n          (in)  Number of interpolation points
C--   jchk       (in)  Yes/no check boundaries

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /   -1,     5   / , iflg /   0,    1  /
      data icma  /    1,     5   /

      logical lint

      dimension w(*), x(*), q(*), f(*)

      character*80 subnam
      data subnam /'EVPLIST ( W, ID, X, QMU2, PDF, N, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check Pdf identifier
      igl = iqcSjekId(subnam,'ID',w,id,icmi,icma,iflg,lint)
*      if(lint) stop 'EVPLIST: cannot access pdfs in internal memory'

C--   Check user input
      if(n.le.0) call sqcErrMsg(subnam,'N should be larger than zero')

C--   Point to the correct set
      ksetw   = abs(igl)/1000
      if(lint) then
        ipver = int(dparGetPar(stor7,ksetw,idipver8))
        call sparParTo5(ipver)
        call sqcInterpList(subnam,stor7,igl,x,q,f,n,jchk)
      else
        ipver = int(dparGetPar(  w  ,ksetw,idipver8))
        call sparParTo5(ipver)
        call sqcInterpList(subnam,  w  ,igl,x,q,f,n,jchk)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void evtable(double *w, int id, double *xx, int nx, double *qq, int nq, double *pdf, int jchk);
C-----------------------------------------------------------------------
CXXHFW  #define fevtable FC_FUNC(evtable,EVTABLE)
CXXHFW    void fevtable(double*,int*,double*,int*,double*,int*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void evtable(double *w, int id, double *xx, int nx, double *qq, int nq, double *pdf, int jchk)
CXXWRP  {
CXXWRP    fevtable(w, &id, xx, &nx, qq, &nq, pdf, &jchk);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================================
      subroutine EvTable(w,id,xx,nx,qq,nq,pdf,jchk)
C     =============================================

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

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /   -1,     5   / , iflg /   0,    1 /
      data icma  /    1,     5   /

      logical lint

      dimension w(*)
      dimension xx(nx), qq(nq), pdf(nx,nq),fff(nx*nq)

      character*80 subnam
      data subnam /'EVTABLE ( W, ID, X, NX, Q, NQ, TABLE, ICHK)'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check Pdf identifier
      igl = iqcSjekId(subnam,'ID',w,id,icmi,icma,iflg,lint)
*      if(lint) stop 'EVTABLE: cannot access pdfs in internal memory'

C--   Check user input
      if(nx.le.0) call sqcErrMsg(subnam,'NX .le. 0 not allowed')
      if(nq.le.0) call sqcErrMsg(subnam,'NQ .le. 0 not allowed')

C--   Point to the correct set
      ksetw  = abs(igl)/1000
      if(lint) then
        ipver = int(dparGetPar(stor7,ksetw,idipver8))
      else
        ipver = int(dparGetPar(  w  ,ksetw,idipver8))
      endif
      call sparParTo5(ipver)
      yma     = ygrid2(iymac5)
      xmi     = exp(-yma)
      xma     = xmaxc2                                    !exclude x = 1
      tmi     = tgrid2(itmic5)
      tma     = tgrid2(itmac5)
      qmi     = exp(tmi)
      qma     = exp(tma)

C--   Findout x-range
      call sqcRange(xx,nx,xmi,xma,aepsi6,ixmi,ixma,ierrx)
      if(ierrx.eq.2) 
     +    call sqcErrMsg(subnam,'X not in strictly ascending order')
      if(jchk.ne.0 .and. (ixmi.ne.1 .or. ixma.ne.nx))
     +    call sqcErrMsg(subnam,'At least one X(i) out of range')

C--   Findout mu2-range
      call sqcRange(qq,nq,qmi,qma,aepsi6,iqmi,iqma,ierrq)
      if(ierrq.eq.2) 
     +    call sqcErrMsg(subnam,'Q not in strictly ascending order')
      if(jchk.ne.0 .and. (iqmi.ne.1 .or. iqma.ne.nq))
     +    call sqcErrMsg(subnam,'At least one Q(i) out of range')

C--   Preset output table
      do iq = 1,nq
        do ix = 1,nx
          pdf(ix,iq) = qnull6
        enddo            
      enddo

C--   x or qmu2 completely out of range
      if(ierrx.ne.0 .or. ierrq.ne.0) return

C--   Now do the work
      nxx = ixma-ixmi+1
      nqq = iqma-iqmi+1
      if(lint) then
        call sqcEvTable(stor7,igl,xx(ixmi),nxx,qq(iqmi),nqq,fff)
      else
        call sqcEvTable(  w  ,igl,xx(ixmi),nxx,qq(iqmi),nqq,fff)
      endif

C--   Now copy linear store fff(i) to pdf(ix,iq)
      i = 0
      do iq = iqmi,iqma
        do ix = ixmi,ixma
          i = i+1
          pdf(ix,iq) = fff(i)
        enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void evpcopy(double *w, int id, double *def, int n, int jset);
C-----------------------------------------------------------------------
CXXHFW  #define fevpcopy FC_FUNC(evpcopy,EVPCOPY)
CXXHFW    void fevpcopy(double*,int*,double*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void evpcopy(double *w, int id, double *def, int n, int jset)
CXXWRP  {
CXXWRP    fevpcopy(w, &id, def, &n, &jset);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine EvPCopy(w,id,def,n,jset)
C     ===================================

C--   Copy EVDGLAP pdfs from toolbox workspace to internal memory
C--
C--   w            (in) : local workspace
C--   id(0:12+n)   (in) : list of global identifiers to copy (0=gluon)
C--   def(-6:6,12) (in) : def(i,j) is amount of flavour i in pdf j
C--   n            (in) : number of extra tables beyond 13
C--   jset         (in) : pdf set number in internal memory [1,mset0]

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

      character*60 emsg
      character*10 etxt

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension icmi(2), icma(2), iflg(2)
C--               isign  itype              pzi  cfil
      data icmi  /    1,     5   / , iflg /   0,    1 /
      data icma  /    1,     5   /

      logical lint

      dimension w(*)
      dimension id(0:12+n), def(-6:6,12)

      character*80 subnam
      data subnam /'EVPCOPY ( W, ID, DEF, N, ISET )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check n is in range
      call sqcIlele(subnam,'N',0,n,mpdf0-13,' ')

C--   Check pdf set in range
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')

C--   Now book jset if it does not already exist
      ntab  = 13+n     !number of pdfs
      ifrst = 0        !first pdf id 0=gluon
      noalf = 0        !book also alfas tables
      call sqcPdfBook(jset,ntab,ifrst,noalf,nwlast,ierr)

C--   Handle error code from sqcPdfBook
      if(ierr.ge.-3) then
        call sqcMemMsg(subnam,nwlast,ierr)
      elseif(ierr.eq.-4) then
        call sqcNtbMsg(subnam,'ISET',jset)
      elseif(ierr.eq.-5) then
        call sqcErrMsg(subnam,'ISET exists but has no pointer tables')
      else
        stop 'EVPCOPY: unkown error code from sqcPdfBook'
      endif

C--1  Get keys
      kset1   = abs(id(0))/1000
      key1    = int(dparGetPar(w    ,kset1,idipver8))
      kset2   = isetf7(jset)
      key2    = int(dparGetPar(stor7,kset2,idipver8))
C--2  Handle counters
      if(key1.ne.key2) then
        call sparCountDn(key2)
        call sparCountUp(key1)
        call sparParAtoB(pars8,key1,stor7,kset2)
      endif

C--   Get max number of flavours
      nfmax  = int(dparGetPar(pars8,key1,idnfmax8))

C--   Check Pdf table identifiers
C--   For gluon, quark, antiquark pdfs
      do i = 0,2*nfmax
        ierr =
     +      iqcSjekId(subnam,'ID(i)',w,id(i),icmi,icma,iflg,lint)
      enddo
C--   For extra Pdfs beyond gluon and quarks
      do i = 13,n
        ierr =
     +      iqcSjekId(subnam,'ID(i)',w,id(i),icmi,icma,iflg,lint)
      enddo

C--   Check that all Pdfs are in the same set
C--   For gluon, quark, antiquark pdfs
      do i = 0,2*nfmax
        kset = abs(id(i))/1000
        if(kset.ne.kset1) call sqcErrMsg(subnam,
     +                   'Not all input ID(i) are in the same set')
      enddo
C--   For extra Pdfs beyond gluon and quarks
      do i = 13,n
        kset = abs(id(i))/1000
        if(kset.ne.kset1) call sqcErrMsg(subnam,
     +                   'Not all input ID(i) are in the same set')
      enddo

C--   Now copy the pdfs from w to stor7
      id0 = iqcIdPdfLtoG(jset,0)        !global gluon id in stor7
      call sqcEvPCopy(w,id,def,n,id0,nfmax,ierr)

C--   Input pdfs not linearly independent
      if(ierr.gt.0) then
        call smb_itoch(2*ierr,etxt,ltxt)
        write(emsg,'(''First '',A,
     &  '' input pdfs not linearly independent'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif

C--   Validate the pdfs
      do i = ifrst7(jset),ilast7(jset)
        idglobal = iqcIdPdfLtoG(jset,i)
        call sqcValidate(stor7,idglobal)
      enddo

C--   Store the type of evolution and max heavy flavour
      call sparSetPar(stor7,kset2,idievtyp8,4.D0)
      call sparSetPar(w,ksetw,idnfheavy8,dble(nfmax))
C--   Set filled
      Lfill7(jset) = .true.
C--   Set parameter key for jset
      ikeyf7(jset) = key1

C--   Update status bits
      call sqcSetflg(iset,idel,jset)

      return
      end
