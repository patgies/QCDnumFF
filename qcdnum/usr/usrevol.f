
C--   This is the file usrevol.f containing the evolution routines

C--   double precision function rfromf(f2)
C--   double precision function ffromr(r2)
C--   function asfunc(r2,nfout,ierr)
C--   double precision function AlTabN(jset,jq,n,ierr)
C--   subroutine EvolFG(kset,func,def,iq0,epsi)
C--   subroutine EvSGNS(ityp,func,isns,n,it0,epsi,ierr)
C--   subroutine PdfCpy(jset1,jset2)
C--   subroutine ExtPdf(func,jset,n,offset,epsi)
C--   subroutine UsrPdf(func,jset,n,offset,epsi)
C--   integer function  NpTabs(jset)
C--   integer function  Ievtyp(jset)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM evolution routines from usrevol.f                */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM evolution routines from usrevol.f                  */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM evolution routines from usrevol.f                  */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------


C==   ==================================================================
C==   Scale transformations ============================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    double rfromf(double fscale2);
C-----------------------------------------------------------------------
CXXHFW  #define frfromf FC_FUNC(rfromf,RFROMF)
CXXHFW    double frfromf(double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double rfromf(double fscale2)
CXXWRP  {
CXXWRP    return frfromf(&fscale2);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ====================================
      double precision function rfromf(f2)
C     ====================================

C--   Convert factorization scale f2 to renormalization scale r2.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*80 subnam
      data subnam /'RFROMF ( F2 )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Transform
      rfromf = aar6*f2 + bbr6

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double ffromr(double rscale2);
C-----------------------------------------------------------------------
CXXHFW  #define fffromr FC_FUNC(ffromr,FFROMR)
CXXHFW    double fffromr(double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double ffromr(double rscale2)
CXXWRP  {
CXXWRP    return fffromr(&rscale2);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ====================================
      double precision function ffromr(r2)
C     ====================================

C--   Convert renormalization scale r2 to factorization scale f2.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*80 subnam
      data subnam /'FFROMR ( R2 )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Transform
      ffromr = (r2-bbr6)/aar6

      return
      end
      
C==   ==================================================================
C==   Alpha-s ==========================================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    double asfunc(double r2, int &nf, int &ierr);
C-----------------------------------------------------------------------
CXXHFW  #define fasfunc FC_FUNC(asfunc,ASFUNC)
CXXHFW    double fasfunc(double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double asfunc(double r2, int &nf, int &ierr)
CXXWRP  {
CXXWRP    return fasfunc(&r2,&nf,&ierr);
CXXWRP  }
C-----------------------------------------------------------------------
            
C     ===============================================
      double precision function asfunc(r2,nfout,ierr)
C     ===============================================

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
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
      data subnam /'ASFUNC ( R2, NFOUT, IERR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      jord = iord6
      alf0 = alfq06
      ralf = q0alf6
      nfix = abs(nfix6)

C--   Update thresholds on the renormalisation scale in the VFNS
      if(nfix.eq.0 .or. nfix.eq.1) call sqcRmass2(qthrs6,rthrs6)
C--   Point to base set
      call sparParTo5(1)
C--   Evolve alfas
      asfunc = dqcAsEvol(r2,ralf,alf0,rthrs6,jord,nfout,ierr)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double altabn(int iset, int iq, int n, int &ierr);
C-----------------------------------------------------------------------
CXXHFW  #define faltabn FC_FUNC(altabn,ALTABN)
CXXHFW    double faltabn(int*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double altabn(int iset, int iq, int n, int &ierr)
CXXWRP  {
CXXWRP    return faltabn(&iset,&iq,&n,&ierr);
CXXWRP  }
C-----------------------------------------------------------------------

C     ================================================
      double precision function AlTabN(jset,jq,n,ierr)
C     ================================================

C--   Returns (as/2pi)^n on the factorisation scale
C--
C--   jset  (in) : pdf set 0=base,1=unpol,2=pol,etc
C--   jq    (in) : mu2 grid point >0 nf = 4,5,6 at iqc,b,t
C--                               <0 nf = 3,4,5 at iqc,b,t
C--   n     (in) :  1, 2, 3 --> a^n in expansion (a, a^2, a^3)
C--                 0,-1,-2 --> a^n in expansion (1, a  , a^2)
C--   ierr (out) :  1  iq close or below Lambda^2
C--                 2  iq outside grid

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
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
      data subnam /'ALTABN ( ISET, IQ, N, IERR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Make sure jset is in range
      call sqcIlele(subnam,'ISET',0,jset,mset0,'ISET does not exist')
C--   Make sure jset is filled (fatal error if not)
      if(.not.Lfill7(jset)) call sqcSetMsg(subnam,'ISET',jset)
C--   Check range
      call sqcIlele(subnam,'N',-2,n,20,' ')
C--   Check iq inside grid 
      iq = abs(jq)
      if(iq.lt.1 .or. iq.gt.ntt2) then
        altabn = 0.D0
        ierr   = 2
        return
      endif
C--   Point to correct set
*mb      write(6,*) 'altabn call sparparto5'
*mb      write(6,*) 'iset,key =', jset, ikeyf7(jset)
      call sparParTo5(ikeyf7(jset))
C--   Not below Lambda, thank you
      if(iq.lt.itmin5) then
        altabn = 0.D0
        ierr   = 1
        return
      endif

C--   Find out which iz index to take
      iz   = izfit5(iq)
      nfiz = itfiz5(-iz)
*mb      iz   = izfitU2(iq)
*mb      nfiz = nffiz2(iz)
*mb      write(6,*) 'iz,nf =',iz1,iz,nfiz1,nfiz

C--   Take iz-1 at threshold if jq preceeded by minus sign
      if(jq.lt.0 .and. iz.ne.1) then
        nfizm1 = itfiz5(1-iz)
        if(nfizm1.eq.nfiz-1) iz=iz-1
      endif
C--   Do the work
*mb      write(6,*) 'call antab8 iz,jset  =',iz,jset
      ierr = 0
      if(n.eq.0)                      then
        altabn = 1.D0
      elseif(n.lt.0)                  then
        altabn = antab8(iz,n,jset)
      elseif(n.gt.0 .and. n.le.iord6) then
        altabn = antab8(iz,n,jset)
      else  
        altabn = (antab8(iz,0,jset))**n
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
          
C==   ==================================================================
C==   PDF evolution ====================================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void evolfg(int iset, double (*func)(int*, double*), double *def, int iq0, double &epsi);
C-----------------------------------------------------------------------
CXXHFW  #define fevolfg FC_FUNC(evolfg,EVOLFG)
CXXHFW    void fevolfg(int*, double (*)(int* ,double*), double*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void evolfg(int iset, double (*func) (int*, double*), double *def, int iq0, doubleCXXWRP &epsi)
CXXWRP  {
CXXWRP    fevolfg(&iset,func,def,&iq0,&epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     =========================================
      subroutine EvolFG(kset,func,def,iq0,epsi)
C     =========================================

C--   Evolve all pdfs.
C-- 
C--   kset      (in) :  10*jset + itype (<0 use old evolfg)
C--   func(j,x) (in) :  Function that returns input pdf_j(x) at iq0
C--   def(i,j)  (in) :  Contribution of flavour (i) to input pdf (j)
C--   iq0       (in) :  Starting scale
C--   epsi     (out) :  Largest deviation quad - lin at midpoint
C--   
C--   Flavor indices :  tb bb cb sb ub db  g  d  u  s  c  b  t
C--                     -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   q+- indices    :   g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--                      0  1  2  3  4  5  6  7  8  9 10 11 12  
C--   Si/ns  indices :   g  si 2+ 3+ 4+ 5+ 6+ va 2- 3- 4- 5- 6-
C--                      0  1  2  3  4  5  6  7  8  9 10 11 12    

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
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

      dimension def(-6:6,12)

      character*80 subnam1,subnam2,subnam
      data subnam1 /'EVOLFG ( ISET, FUNC, DEF, IQ0, EPSI )'/
      data subnam2 /'EVOLFG-OLD ( ISET, FUNC, DEF, IQ0, EPSI )'/

      character*38 wtmsg(4)
C--                12345678901234567890123456789012345678
      data wtmsg /'No unpolarised weight tables available',
     +            'No polarised weight tables available  ',
     +            'No time-like weight tables available  ',
     +            'No custom weight tables available     '/

C--   Old or new routine thats the question
      if(kset.lt.0) then
        subnam = subnam2
      else
        subnam = subnam1
      endif

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam2,ichk,iset,idel)
        first = .false.
      endif

C--   Get jset and ityp
      jset = abs(kset)/10
      ityp = abs(kset)-10*jset
      if(jset.eq.0) jset = ityp

C--   Check jset is in range
      call sqcIlele(subnam,'ISET',1,jset,mset0,
     +     'Invalid PDF set identifier')

C--   Check ityp is in range
      call sqcIlele(subnam,'ITYPE',1,ityp,3,
     + 'ITYPE must be unpolarised (1), polarised (2) or time-like (3)')

C--   Check status bits
      call sqcChkflg(jset,ichk,subnam)

C--   Check weights available
      if(isetp7(ityp).eq.0) call sqcErrMsg(subnam,wtmsg(ityp))

C--   Book pdf and alfas tables in stor7 (if they do not exist already)
      ntab  = 13
      ifrst = 0
      noalf = 0
      call sqcPdfBook(jset,ntab,ifrst,noalf,nwlast,ierr)

C--   Handle error code from sqcPdfBook
      if(ierr.ge.-3) then
        call sqcMemMsg(subnam,nwlast,ierr)
      elseif(ierr.eq.-4) then
        call sqcNtbMsg(subnam,'ISET',jset)
      elseif(ierr.eq.-5) then
        call sqcErrMsg(subnam,'ISET exists but has no pointer tables')
      else
        stop 'EvolFG: unkown error code from sqcPdfBook'
      endif

C--1  Get pdf_key
      kset0   = isetf7(0)
      ksetw   = isetf7(jset)
      keypdf  = int(dparGetPar(stor7,ksetw,idipver8))
      if(keypdf.lt.0 .or. keypdf.gt.mset0) stop 'EVOLFG: invalid key'
C--2  Get base_key
      keybase = iparGetGroupKey(pars8,1,6)
C--3  Do nothing if pdf_key = base_key
      if(keypdf.ne.keybase) then
        call sparCountDn(keypdf)
        call sparBaseToKey(keybase)
        call sparCountUp(keybase)
        call sparParAtoB(pars8,keybase,stor7,ksetw)
        call sparAlfAtoB(pars8,keybase,stor7,kset0)  
      endif

C--   Point to base set
      call sparParTo5(1)

C--   Check perturbative order
      call sqcIlele(subnam,'IORD',1,iord6,mxord7(ityp),' ')

C--   Check starting scale iq0
      if(abs(nfix6).eq.0 .and. kset.lt.0) then
C--     Old evolution code in the VFNS
        call sqcIlele(subnam,'IQ0',1,abs(iq0),itchm2-1,
     +               'IQ0 should be below the charm threshold')
      elseif(abs(nfix6).eq.1 .and. kset.lt.0) then
C--     Old evolution cannot run with intrinsic heavy flavours
        call sqcErrMsg(subnam,
     +                'Cannot evolve with intrinsic heavy flavours')
      else
C--     Just check that iq0 is within the grid boundaries
        call sqcIlele(subnam,'IQ0',1,abs(iq0),ntt2,
     +               'IQ0 outside the grid boundaries')
      endif

C--   Dummy call to func
      dummy = func(-1,0.D0)

C--   Go..
      if(kset.gt.0) then
        call sqcEvolFG(ityp,jset,func,def,iq0,epsi,nfheavy,ierr)
      else
        call sqcEvolFG_old(ityp,jset,func,def,iq0,epsi,nfheavy,ierr)
      endif

C--   Handle error from EvolFG
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,
     +       'IQ0 outside the grid boundaries or cuts')
      elseif(ierr.eq.2) then
        call sqcErrMsg(subnam,
     +       'Attempt to evolve with too large alpha-s')
      elseif(ierr.eq.3) then
        call sqcErrMsg(subnam,
     +       'Input pdfs not linearly independent')
      elseif(ierr.eq.4) then
        call sqcErrMsg(subnam,
     +  'Intrinsic heavy quark input must be mixture of h and hbar')
      endif

C--   Check max deviation
      if(dflim6.gt.0.D0 .and. epsi.gt.dflim6) call sqcErrMsg(subnam, 
     +          'Possible spline oscillation detected')

C--   Store the type of evolution and nfheavy
      call sparSetPar(stor7,ksetw,idnfheavy8,dble(nfheavy))
      call sparSetPar(stor7,ksetw,idievtyp8,dble(ityp))
      itypf7(jset) = ityp
C--   Set filled
      Lfill7(jset) = .true.
C--   Store parameter key for jset
      ikeyf7(jset) = keybase

C--   Update status bits
      call sqcSetflg(iset,idel,jset)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void evsgns(int iset, double (*func)(int*, double*), int *isns, int n, int iq0, double &epsi);
C-----------------------------------------------------------------------
CXXHFW  #define fevsgns FC_FUNC(evsgns,EVSGNS)
CXXHFW    void fevsgns(int*, double (*)(int* ,double*), int*, int*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void evsgns(int iset, double (*func)(int*, double*), int *isns, int n, int iq0, double &epsi)
CXXWRP  {
CXXWRP    fevsgns(&iset,func,isns,&n,&iq0,&epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     ============================================
      subroutine EVSGNS(kset,func,isns,n,iq0,epsi)
C     ============================================

C--   Evolve all pdfs.
C--
C--   kset      (in) :  10*jset + itype (<0 use old evolfg)
C--   func(j,x) (in) :  Function that returns input pdf_j(x) at iq0
C--   isns(i)   (in) :  si/ns specifiers 1=SG 2=NS+ -1=V -2=NS-
C--   n         (in) :  number of pdfs to evolve
C--   iq0       (in) :  Starting scale
C--   epsi     (out) :  Largest deviation quad - lin at midpoint
C--
C--   Flavor indices :  tb bb cb sb ub db  g  d  u  s  c  b  t
C--                     -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   q+- indices    :   g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--                      0  1  2  3  4  5  6  7  8  9 10 11 12
C--   Si/ns  indices :   g  si 2+ 3+ 4+ 5+ 6+ va 2- 3- 4- 5- 6-
C--                      0  1  2  3  4  5  6  7  8  9 10 11 12

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'pstor8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension isns(*)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam  /'EVSGNS ( ISET, FUNC, ISNS, N, IQ0, EPSI )'/

      character*38 wtmsg(4)
C--                12345678901234567890123456789012345678
      data wtmsg /'No unpolarised weight tables available',
     +            'No polarised weight tables available  ',
     +            'No time-like weight tables available  ',
     +            'No custom weight tables available     '/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Get jset and ityp
      jset = abs(kset)/10
      ityp = abs(kset)-10*jset
      if(jset.eq.0) jset = ityp

C--   Check jset is in range
      call sqcIlele(subnam,'ISET',1,jset,mset0,
     +     'Invalid PDF set identifier')

C--   Check ityp is in range
      call sqcIlele(subnam,'ITYPE',1,ityp,3,
     + 'ITYPE must be unpolarised (1), polarised (2) or time-like (3)')

C--   Check status bits
      call sqcChkflg(jset,ichk,subnam)

C--   Check N  in range
      call sqcIlele(subnam,'N',1,n,mpdf0-1,' ')

C--   Check weights available
      if(isetp7(ityp).eq.0) call sqcErrMsg(subnam,wtmsg(ityp))

C--   Book pdf and alfas tables in stor7 (if they do not exist already)
      ntab  = n+1
      ifrst = 0
      noalf = 0
      call sqcPdfBook(jset,ntab,ifrst,noalf,nwlast,ierr)

C--   Handle error code from sqcPdfBook
      if(ierr.ge.-3) then
        call sqcMemMsg(subnam,nwlast,ierr)
      elseif(ierr.eq.-4) then
        call sqcNtbMsg(subnam,'ISET',jset)
      elseif(ierr.eq.-5) then
        call sqcErrMsg(subnam,'ISET exists but has no pointer tables')
      else
        stop 'EVSGNS: unkown error code from sqcPdfBook'
      endif

C--1  Get pdf_key
      kset0   = isetf7(0)
      ksetw   = isetf7(jset)
      keypdf  = int(dparGetPar(stor7,ksetw,idipver8))
      if(keypdf.lt.0 .or. keypdf.gt.mset0) stop 'EVSGNS: invalid key'
C--2  Get base_key
      keybase = iparGetGroupKey(pars8,1,6)
C--3  Do nothing if pdf_key = base_key
      if(keypdf.ne.keybase) then
        call sparCountDn(keypdf)
        call sparBaseToKey(keybase)
        call sparCountUp(keybase)
        call sparParAtoB(pars8,keybase,stor7,ksetw)
        call sparAlfAtoB(pars8,keybase,stor7,kset0)
      endif

C--   Point to base set
      call sparParTo5(1)

C--   Check perturbative order
      call sqcIlele(subnam,'IORD',1,iord6,mxord7(ityp),' ')

C--   Only FFNS or MFNS
      call  sqcIlele(subnam,'NFIX',3,abs(nfix6),6,
     +   'Can only evolve in the FFNS or MFNS')

C--   Check starting scale iq0
      call sqcIlele(subnam,'IQ0',1,abs(iq0),ntt2,
     +                   'IQ0 outside the grid boundaries')

C--   Check isns
      do i = 1,n
        call sqcIlele(subnam,'ISNS(i)',1,abs(isns(i)),2,
     +   'Invalid singlet-nonsinglet specifier')
      enddo
      do i = 2,n
        if(isns(i).eq.1) then
          call sqcErrMsg(subnam,
     +   'Found Singlet/Gluon specifier not in position ISNS(1)' )
        endif
      enddo

C--   Dummy call to func
      dummy = func(-1,0.D0)

C--   Go..
      call sqcEvSGNS(ityp,jset,func,isns,n,iq0,epsi,nfheavy,ierr)

C--   Handle error from EvSGNS
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,
     +       'IQ0 outside the grid boundaries or cuts')
      elseif(ierr.eq.2) then
        call sqcErrMsg(subnam,
     +       'Attempt to evolve with too large alpha-s')
      endif

C--   Check max deviation
      if(dflim6.gt.0.D0 .and. epsi.gt.dflim6) call sqcErrMsg(subnam,
     +          'Possible spline oscillation detected')

C--   Store the type of evolution
      call sparSetPar(stor7,ksetw,idnfheavy8,dble(nfheavy))
      call sparSetPar(stor7,ksetw,idievtyp8,5.D0)
      itypf7(jset) = 5
C--   Set filled
      Lfill7(jset) = .true.
C--   Store parameter key for jset
      ikeyf7(jset) = keybase

C--   Update status bits
      call sqcSetflg(iset,idel,jset)

      return
      end

C==   ==================================================================
C==   Pdfs from elsewhere   ============================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void pdfcpy(int iset1, int iset2);
C-----------------------------------------------------------------------
CXXHFW  #define fpdfcpy FC_FUNC(pdfcpy,PDFCPY)
CXXHFW    void fpdfcpy(int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void pdfcpy(int iset1, int iset2)
CXXWRP  {
CXXWRP    fpdfcpy(&iset1,&iset2);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================
      subroutine PdfCpy(jset1,jset2)
C     ==============================

C--   Copy jset1 [1,mset0] to jset2 [1,mset0]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'
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

      character*80 subnam
      data subnam /'PDFCPY ( ISET1, ISET2 )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check jset1,2
      call sqcIlele(subnam,'ISET1',1,jset1,mset0,' ')
      call sqcIlele(subnam,'ISET2',1,jset2,mset0,' ')

C--   Check status bits  (also check if pdf set is filled)
      call sqcChkflg(jset1,ichk,subnam)

C--   Check jset1 is filled
      if(.not.Lfill7(jset1)) call sqcSetMsg(subnam,'ISET1',jset1)

C--   Do nothing if jset1 = jset2
      if(jset1.eq.jset2) return

C--   Now book jset2 if it does not already exist
      ntab  = ilast7(jset1)-ifrst7(jset1)+1
      ifrst = 0
      noalf = 0
      call sqcPdfBook(jset2,ntab,ifrst,noalf,nwlast,ierr)

C--   Handle error code from sqcPdfBook
      if(ierr.ge.-3) then
        call sqcMemMsg(subnam,nwlast,ierr)
      elseif(ierr.eq.-4) then
        call sqcNtbMsg(subnam,'ISET',jset2)
      elseif(ierr.eq.-5) then
        call sqcErrMsg(subnam,'ISET exists but has no pointer tables')
      else
        stop 'PdfCpy unkown error code from sqcPdfBook'
      endif

C--1  Get keys
      kset1   = isetf7(jset1)
      key1    = int(dparGetPar(stor7,kset1,idipver8))
      kset2   = isetf7(jset2)
      key2    = int(dparGetPar(stor7,kset2,idipver8))
C--2  Handle counters
      if(key1.eq.key2) then
        call sqcPdfCpy(kset1,kset2)                                !copy
      else
        call sparCountDn(key2)
        call sparCountUp(key1)
        call sqcPdfCpy(kset1,kset2)                                !copy
        call sparParAtoB(pars8,key1,stor7,kset2)
        evtyp  = dparGetPar(stor7,kset1,idievtyp8)
        fheavy = dparGetPar(stor7,kset1,idnfheavy8)
        call sparSetPar(stor7,kset2,idievtyp8,evtyp)
        call sparSetPar(stor7,kset2,idnfheavy8,fheavy)
      endif
C--   Set filled
      Lfill7(jset2) = .true.
C--   Set parameter key of jset2
      ikeyf7(jset2) = ikeyf7(jset1)
      itypf7(jset2) = itypf7(jset1)

C--   Update status bits
      call sqcSetflg(iset,idel,jset2)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void extpdf(double (*func)(int*, double*, double*, bool*), int iset, int n, double offset, double &epsi);
C-----------------------------------------------------------------------
CXXHFW  #define fextpdf FC_FUNC(extpdf,EXTPDF)
CXXHFW    void fextpdf(double (*)(int*, double* ,double* ,bool*), int*, int*, double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void extpdf(double (*func)(int*, double*, double*, bool*), int iset, int n, double offset, double &epsi)
CXXWRP  {
CXXWRP    fextpdf(func,&iset,&n,&offset,&epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================================
      subroutine ExtPdf(func,jset,n,offset,epsi)
C     ==========================================

C--   Read pdfs into memory
C--
C--   func    (external) user supplied function
C--   jset    (in)       pdf set [1,mset0]
C--   n       (in)       number of pdfs beyond gluon,q,qbar
C--   offset  (in)       threshold offset mu2h*(1+-offset)
C--   epsi    (out)      estimate of spline accuracy
C--
C--   subr(x,mu2,qqbar) should return qqbar(-6:6+n)= ...,ub,db,gl,d,u,..

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'
      include 'point5.inc'
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

      character*80 subnam
      data subnam /'EXTPDF ( FUNC, ISET, N, DELTA, EPSI )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check pdf set in range
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Check extra tables are in range
      call sqcIlele(subnam,'N',0,n,mpdf0-13,' ')

C--   Check status bits  (also check if pdf set is filled)
      call sqcChkflg(jset,ichk,subnam)

C--   Invalidate pdfs in jset
C--   This blocks access to jset via subr
      call sqcDelbit(ibPdfs4,istat4(1,jset),mbp0)

C--   Point to base set
      call sparParTo5(1)                                !point to params

C--   Initialize store
      if(.not.Lwtini7) call sqcIniWt

C--   Now book jset if it does not already exist
      ntab  = 13+n
      ifrst = 0
      noalf = 0
      call sqcPdfBook(jset,ntab,ifrst,noalf,nwlast,ierr)

C--   Handle error code from sqcPdfBook
      if(ierr.ge.-3) then
        call sqcMemMsg(subnam,nwlast,ierr)
      elseif(ierr.eq.-4) then
        call sqcNtbMsg(subnam,'ISET',jset)
      elseif(ierr.eq.-5) then
        call sqcErrMsg(subnam,'ISET exists but has no pointer tables')
      else
        stop 'EXTPDF: unkown error code from sqcPdfBook'
      endif

C--1  Find current key of jset (this is to be changed later)
      ksetw   = isetf7(jset)
      keypdf  = int(dparGetPar(stor7,ksetw,idipver8))
      if(keypdf.lt.0 .or. keypdf.gt.mpl0) stop 'EXTPDF: invalid key'
C--2  Get base_key
      keybase = iparGetGroupKey(pars8,1,6)
C--3  Do nothing if pdf_key = base_key
      if(keypdf.ne.keybase) then
        call sparCountDn(keypdf)
        call sparBaseToKey(keybase)
        call sparCountUp(keybase)
        call sparParAtoB(pars8,keybase,stor7,ksetw)
      endif

C--   Global id of gluon table
      ig0  = iqcIdPdfLtoG(jset, 0)

C--   Offset cut-off
      off = max(abs(offset),2.D0*aepsi6)

C--   Now import pdf set as defined by func
      call sqcExtPdf(func,ig0,n,off,nfheavy)                       !fill

C--   Check for spline oscillations
      epsi = 0.D0
      do id = ifrst7(jset),ilast7(jset)
        idg   = iqcIdPdfLtoG(jset,id)
        iq1   = itfiz5(izmic2)
        iq2   = itfiz5(izmac2)
        do iq = iq1,iq2
          eps  = dqcSplChk(idg,iq)
          epsi = max(epsi,eps)
        enddo
      enddo

C--   Check max deviation
      if(dflim6.gt.0.D0 .and. epsi.gt.dflim6) call sqcErrMsg(subnam,
     +          'Possible spline oscillation detected')

C--   Global identifiers
      ksetw  = isetf7(jset)
      igfrst = iqcIdPdfLtoG(jset,ifrst7(jset))
      iglast = iqcIdPdfLtoG(jset,ilast7(jset))
      idmax  = iqcGetNumberOfTables(stor7,ksetw,5)+ifrst7(jset)-1
      igmax  = iqcIdPdfLtoG(jset,idmax)

C--   Validate the pdfs, set iset pointer and store the cuts
      do idg = igfrst,iglast
        call sqcValidate(stor7,idg)
      enddo
C--   Invalidate pdfs that were not filled by sqcPdfExt (if any)
      do idg = iglast+1,igmax
        call sqcInvalidate(stor7,idg)
      enddo

C--   Store the type of evolution
      call sparSetPar(stor7,ksetw,idnfheavy8,dble(nfheavy))
      call sparSetPar(stor7,ksetw,idievtyp8,4.D0)
      itypf7(jset) = 4

C--   Set filled
      Lfill7(jset) = .true.
C--   Set parameter key for jset
      ikeyf7(jset) = keybase

C--   Update status bits
      call sqcSetflg(iset,idel,jset)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void usrpdf(double (*func)(int*, double*, double*, bool*), int iset, int n, double offset, double &epsi);
C-----------------------------------------------------------------------
CXXHFW  #define fusrpdf FC_FUNC(usrpdf,USRPDF)
CXXHFW    void fusrpdf(double (*)(int*, double* ,double* ,bool*), int*, int*, double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void usrpdf(double (*func)(int*, double*, double*, bool*), int iset, int n, double offset, double &epsi)
CXXWRP  {
CXXWRP    fusrpdf(func,&iset,&n,&offset,&epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================================
      subroutine UsrPdf(func,jset,n,offset,epsi)
C     ==========================================

C--   Read pdfs into memory
C--
C--   func    (external) user supplied function
C--   jset    (in)       pdf set [1,mset0]
C--   n       (in)       number of pdfs beyond gluon
C--   offset  (in)       threshold offset mu2h*(1+-offset)
C--   epsi    (out)      estimate of spline accuracy
C--
C--   subr(x,mu2,qqbar) should return qqbar(-6:6+n)= ...,ub,db,gl,d,u,..

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'
      include 'point5.inc'
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

      character*80 subnam
      data subnam /'USRPDF ( FUNC, ISET, N, DELTA, EPSI )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check pdf set in range
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Check extra tables are in range
      call sqcIlele(subnam,'N',1,n,mpdf0-1,' ')

C--   Check status bits  (also check if pdf set is filled)
      call sqcChkflg(jset,ichk,subnam)

C--   Invalidate pdfs in jset
C--   This blocks access to jset via subr
      call sqcDelbit(ibPdfs4,istat4(1,jset),mbp0)

C--   Point to base set
      call sparParTo5(1)                                !point to params

C--   Initialize store
      if(.not.Lwtini7) call sqcIniWt

C--   Now book jset if it does not already exist
      ntab  = 1+n
      ifrst = 0
      noalf = 0
      call sqcPdfBook(jset,ntab,ifrst,noalf,nwlast,ierr)

C--   Handle error code from sqcPdfBook
      if(ierr.ge.-3) then
        call sqcMemMsg(subnam,nwlast,ierr)
      elseif(ierr.eq.-4) then
        call sqcNtbMsg(subnam,'ISET',jset)
      elseif(ierr.eq.-5) then
        call sqcErrMsg(subnam,'ISET exists but has no pointer tables')
      else
        stop 'USRPDF: unkown error code from sqcPdfBook'
      endif

C--   New code
C--1  Find current key of jset (this is to be changed later)
      ksetw   = isetf7(jset)
      keypdf  = int(dparGetPar(stor7,ksetw,idipver8))
      if(keypdf.lt.0 .or. keypdf.gt.mpl0) stop 'USRPDF: invalid key'
C--2  Get base_key
      keybase = iparGetGroupKey(pars8,1,6)
C--3  Do nothing if pdf_key = base_key
      if(keypdf.ne.keybase) then
        call sparCountDn(keypdf)
        call sparBaseToKey(keybase)
        call sparCountUp(keybase)
        call sparParAtoB(pars8,keybase,stor7,ksetw)
      endif

C--   Global id of gluon table
      ig0  = iqcIdPdfLtoG(jset, 0)

C--   Offset cut-off
      off = max(abs(offset),2.D0*aepsi6)

C--   Now import pdf set as defined by func
      call sqcUsrPdf(func,ig0,n,off,nfheavy)                       !fill

C--   Check for spline oscillations
      epsi = 0.D0
      do id = ifrst7(jset),ilast7(jset)
        idg   = iqcIdPdfLtoG(jset,id)
        iq1   = itfiz5(izmic2)
        iq2   = itfiz5(izmac2)
        do iq = iq1,iq2
          eps  = dqcSplChk(idg,iq)
          epsi = max(epsi,eps)
        enddo
      enddo

C--   Check max deviation
      if(dflim6.gt.0.D0 .and. epsi.gt.dflim6) call sqcErrMsg(subnam,
     +          'Possible spline oscillation detected')

C--   Global identifiers
      ksetw  = isetf7(jset)
      igfrst = iqcIdPdfLtoG(jset,ifrst7(jset))
      iglast = iqcIdPdfLtoG(jset,ilast7(jset))
      idmax  = iqcGetNumberOfTables(stor7,ksetw,5)+ifrst7(jset)-1
      igmax  = iqcIdPdfLtoG(jset,idmax)

C--   Validate the pdfs, set iset pointer and store the cuts
      do idg = igfrst,iglast
        call sqcValidate(stor7,idg)
      enddo
C--   Invalidate pdfs that were not filled by sqcPdfExt (if any)
      do idg = iglast+1,igmax
        call sqcInvalidate(stor7,idg)
      enddo

C--   Store the type of evolution
      call sparSetPar(stor7,ksetw,idnfheavy8,dble(nfheavy))
      call sparSetPar(stor7,ksetw,idievtyp8,5.D0)
      itypf7(jset) = 5

C--   Set filled
      Lfill7(jset) = .true.
C--   Set parameter key for jset
      ikeyf7(jset) = keybase

C--   Update status bits
      call sqcSetflg(iset,idel,jset)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int nptabs(int iset);
C-----------------------------------------------------------------------
CXXHFW  #define fnptabs FC_FUNC(nptabs,NPTABS)
CXXHFW    int fnptabs(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int nptabs(int iset)
CXXWRP  {
CXXWRP    return fnptabs(&iset);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      integer function NpTabs(jset)
C     =============================

C--   Return the number of pdf tables in jset

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*80 subnam
      data subnam /'NPTABS ( ISET )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif
C--   Check jset in range
      call sqcIlele(subnam,'ISET',1,jset,mset0,'ISET does not exist')
C--   Go ...
      if(Lfill7(jset)) then
        NpTabs = ilast7(jset)-ifrst7(jset)+1
      else
        NpTabs = 0
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int ievtyp(int iset);
C-----------------------------------------------------------------------
CXXHFW  #define fievtyp FC_FUNC(ievtyp,IEVTYP)
CXXHFW    int fievtyp(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int ievtyp(int iset)
CXXWRP  {
CXXWRP    return fievtyp(&iset);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      integer function Ievtyp(jset)
C     =============================

C--   Return the evolution type of jset

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*80 subnam
      data subnam /'IEVTYP ( ISET )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif
C--   Check jset in range
      call sqcIlele(subnam,'ISET',1,jset,mset0,'ISET does not exist')
C--   Go ...
      if(Lfill7(jset)) then
        Ievtyp = int(dparGetPar(stor7,isetf7(jset),idievtyp8))
      else
        Ievtyp = 0
      endif

      return
      end
