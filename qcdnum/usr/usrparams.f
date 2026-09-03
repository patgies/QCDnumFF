
C--   This is the file usrparams.f with evolution parameter stuff

C--   subroutine setord(iord)
C--   subroutine getord(iord)
C--   subroutine setalf(as,22)
C--   subroutine getalf(as,r2)
C--   subroutine setcbt(nfix,iqc,iqb,iqt)
c--   subroutine mixfns(nfix,r2c,r2b,r2t)
C--   subroutine getcbt(nfix,iqc,iqb,iqt)
C--   subroutine setabr(ar,br)
C--   subroutine getabr(ar,br)
C--
C--   integer function nfrmiq(jset,iq,ithrs)
C--   integer function nflavs(iq,ithrs)
C--
C--   subroutine CpyPar(array,n,jset)
C--   subroutine UsePar(jset)
C--   integer function KeyPar(jset)
C--   integer function KeyGrp(jset,igroup)
C--   subroutine Pushcp
C--   subroutine Pullcp

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM parameter routines from usrparams.f              */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM parameter routines from usrparams.f                */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM parameter routines from usrparams.f                */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------


C==   ==================================================================
C==   Routines to set and get evolution parameters =====================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void setord(int iord);
C-----------------------------------------------------------------------
CXXHFW  #define fsetord FC_FUNC(setord,SETORD)
CXXHFW    void fsetord(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setord(int iord)
CXXWRP  {
CXXWRP    fsetord(&iord);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================
      subroutine setord(iord)
C     =======================

C--   Set order 1,2,3 = LO, NLO, NNLO.

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

      character*80 subnam
      data subnam /'SETORD ( IORD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Already set?
      if(iord.eq.iord6)   return
C--   Check user input
      call sqcIlele(subnam,'IORD',1,iord,3,' ')
C--   Do the work
      iord6   = iord
C--   Invalidate parameter store
      call smb_sbit1(ipbits8,ipsbit8)
C--   Invalidate alphas table
      call smb_sbit1(ipbits8,iasbit8)
C--   Construct base parameter set
      call sparMakeBase

C--   Update status bits
      call sqcSetflg(iset,idel,0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void getord(int &iord);
C-----------------------------------------------------------------------
CXXHFW  #define fgetord FC_FUNC(getord,GETORD)
CXXHFW    void fgetord(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void getord(int &iord)
CXXWRP  {
CXXWRP    fgetord(&iord);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================
      subroutine getord(iord)
C     =======================

C--   Get current value of iord = 1,2,3 for LO, NLO, NNLO.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*80 subnam
      data subnam /'GETORD ( IORD )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Do the work
      iord = iord6

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void setalf(double as0, double r20);
C-----------------------------------------------------------------------
CXXHFW  #define fsetalf FC_FUNC(setalf,SETALF)
CXXHFW    void fsetalf(double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setalf(double as0, double r20)
CXXWRP  {
CXXWRP    fsetalf(&as0,&r20);
CXXWRP  }
C-----------------------------------------------------------------------

C     ========================
      subroutine setalf(as,r2)
C     ========================

C--   Set input value of alpha_s(r2)

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

      character*80 subnam
      data subnam /'SETALF ( AS, R2 )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Already set?
      if(as.eq.alfq06 .and. r2.eq.q0alf6)   return
C--   Check user input
      call sqcDlele(subnam,'AS',1.D-5,as,aslim6,
     + 'Remark: the upper AS limit can be changed by a call to SETVAL')
      call sqcDlele(subnam,'R2',qlimd6,abs(r2),qlimu6,
     + 'Remark: these R2 limits can be changed by a call to SETVAL')
C--   Do the work
      alfq06  = as
      q0alf6  = r2
C--   Invalidate parameter store
      call smb_sbit1(ipbits8,ipsbit8)
C--   Invalidate alphas table
      call smb_sbit1(ipbits8,iasbit8)
C--   Construct base parameter set
      call sparMakeBase

C--   Update status bits
      call sqcSetflg(iset,idel,0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void getalf(double &as0, double &r20);
C-----------------------------------------------------------------------
CXXHFW  #define fgetalf FC_FUNC(getalf,GETALF)
CXXHFW    void fgetalf(double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void getalf(double &as0, double &r20)
CXXWRP  {
CXXWRP    fgetalf(&as0,&r20);
CXXWRP  }
C-----------------------------------------------------------------------

C     ========================
      subroutine getalf(as,r2)
C     ========================

C--   Get input value of alpha_s(R2)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*80 subnam
      data subnam /'GETALF ( AS, R2 )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Do the work
      as = alfq06
      r2 = q0alf6

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void setcbt(int nfix, int iqc, int iqb, int iqt);
C-----------------------------------------------------------------------
CXXHFW  #define fsetcbt FC_FUNC(setcbt,SETCBT)
CXXHFW    void fsetcbt(int*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setcbt(int nfix, int iqc, int iqb, int iqt)
CXXWRP  {
CXXWRP    fsetcbt(&nfix,&iqc,&iqb,&iqt);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine setcbt(nfix,iqc,iqb,iqt)
C     ===================================

C--   If 3 .le. nfix .le. 6 set FFNS, nfix = 0,1 set VFNS with
C--   c, b and t thresholds defined on the factorization scale.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpars6.inc'
      include 'qpdfs7.inc'
      include 'pstor8.inc'

      dimension iqh(4:6), iq1(3:6), iq2(3:6)

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam, etxt
      data subnam /'SETCBT ( NFIX, IQC, IQB, IQT )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first   = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check iqc,b,t
      iqh(4) = iqc
      iqh(5) = iqb
      iqh(6) = iqt
      call sqcChkIqh(ntt2,nfix,iqh,iq1,iq2,nfmin,nfmax,ierr)

C--   Error messages
      if(ierr.eq.1) then
        call smb_itoch(nfix,etxt,ltxt)
        call sqcErrMsg(subnam,
     +      'NFIX = '//etxt(1:ltxt)//
     +      ' must be 0,1 (VFNS) or 3,4,5,6 (FFNS)')
      elseif(ierr.eq.2) then
        call sqcErrMsg(subnam,
     +      'None of the IQC, IQB, IQT are inside the grid')
      elseif(ierr.eq.3) then
        call sqcErrMsg(subnam,
     +      'Threshold combination Charm-Top not allowed')
      elseif(ierr.eq.4) then
        call sqcErrMsg(subnam,
     +      'Found thresholds not ascending or too close together')
      endif

C--   Now correctly fill the /qpars6/ common block
      if(nfix.eq.0 .or. nfix.eq.1) then
C--     VFNS
        call sqcThrVFNS(nfix,iqh,nfmin,nfmax)
      else
C--     FFNS
        call sqcThrFFNS(nfix)
      endif

C--   Invalidate parameter store
      call smb_sbit1(ipbits8,ipsbit8)
C--   Invalidate flavour map
      call smb_sbit1(ipbits8,infbit8)
C--   Invalidate alphas map
      call smb_sbit1(ipbits8,iasbit8)
C--   Invalidate iz-cuts
      call smb_sbit1(ipbits8,izcbit8)
C--   Construct base parameter set
      call sparMakeBase

C--   Update status bits
      call sqcSetflg(iset,idel,0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void mixfns(int nfix, double r2c, double r2b, double r2t);
C-----------------------------------------------------------------------
CXXHFW  #define fmixfns FC_FUNC(mixfns,MIXFNS)
CXXHFW    void fmixfns(int*, double*, double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void mixfns(int nfix, double r2c, double r2b, double r2t)
CXXWRP  {
CXXWRP    fmixfns(&nfix,&r2c,&r2b,&r2t);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine mixfns(nfix,r2c,r2b,r2t)
C     ===================================

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      dimension thrin(4:6), throut(4:6)

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'MIXFNS ( NFIX, R2C, R2B, R2T )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first   = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check nfix
      call sqcIlele(subnam,'NFIX',3,nfix,6,' ')

C--   Check thresholds and return ordered threshold list in throut
      thrin(4) = r2c
      thrin(5) = r2b
      thrin(6) = r2t
      qmin     = exp(tgrid2(1))
      qmax     = exp(tgrid2(ntt2))
      call sqcChkRqh(qmin,qmax,thrin,throut,ierr)

C--   Error messages
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,
     +      'None of the R2C, R2B, R2T are inside the grid')
      elseif(ierr.eq.2) then
        call sqcErrMsg(subnam,
     +      'Threshold combination (R2C,xxx,R2T) not allowed')
      elseif(ierr.eq.3) then
        call sqcErrMsg(subnam,
     +      'Found thresholds not ascending or too close together')
      endif

C--   Now fill the /qpars6/ common block
      call sqcThrMFNS(nfix,throut(4),throut(5),throut(6))

C--   Invalidate parameter store
      call smb_sbit1(ipbits8,ipsbit8)
C--   Invalidate flavour map
      call smb_sbit1(ipbits8,infbit8)
C--   Invalidate alphas map
      call smb_sbit1(ipbits8,iasbit8)
C--   Invalidate iz-cuts
      call smb_sbit1(ipbits8,izcbit8)
C--   Construct base parameter set
      call sparMakeBase

C--   Update status bits
      call sqcSetflg(iset,idel,0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void getcbt(int &nfix, double &q2c, double &q2b, double &q2t);
C-----------------------------------------------------------------------
CXXHFW  #define fgetcbt FC_FUNC(getcbt,GETCBT)
CXXHFW    void fgetcbt(int*, double*, double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void getcbt(int &nfix, double &q2c, double &q2b, double &q2t)
CXXWRP  {
CXXWRP    fgetcbt(&nfix,&q2c,&q2b,&q2t);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine getcbt(nfix,q2c,q2b,q2t)
C     ===================================

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'GETCBT ( NFIX, Q2C, Q2B, Q2T )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first   = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Do the work
      nfix   = nfix6
      if(nfix6.lt.0) then
C--     MFNS
        q2c  =  rthrs6(4)
        q2b  =  rthrs6(5)
        q2t  =  rthrs6(6)
      else
C--     FFNS or VFNS
        q2c  =  qthrs6(4)
        q2b  =  qthrs6(5)
        q2t  =  qthrs6(6)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void setabr(double ar, double br);
C-----------------------------------------------------------------------
CXXHFW  #define fsetabr FC_FUNC(setabr,SETABR)
CXXHFW    void fsetabr(double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setabr(double ar, double br)
CXXWRP  {
CXXWRP    fsetabr(&ar,&br);
CXXWRP  }
C-----------------------------------------------------------------------

C     ========================
      subroutine setabr(ar,br)
C     ========================

C--   Set mu2r = ar*mu2f + br

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'SETABR ( AR, BR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first   = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Already set?
      if(ar.eq.aar6 .and. br.eq.bbr6)   return
C--   Check user input
      call sqcDlele(subnam,'AR',1.D-2,ar,1.D2,' ')
      call sqcDlele(subnam,'BR',-1.D2,br,1.D2,' ')
C--   Do the work
      aar6 = ar
      bbr6 = br
C--   Update thresholds on the renormalisation scale in the VFNS
      if(abs(nfix6).eq.0 .or. abs(nfix6).eq.1)
     +                    call sqcRmass2(qthrs6,rthrs6)

C--   Invalidate parameter store
      call smb_sbit1(ipbits8,ipsbit8)
C--   Invalidate alphas map
      call smb_sbit1(ipbits8,iasbit8)
C--   Construct base parameter set
      call sparMakeBase

C--   Update status bits
      call sqcSetflg(iset,idel,0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void getabr(double &ar, double &br);
C-----------------------------------------------------------------------
CXXHFW  #define fgetabr FC_FUNC(getabr,GETABR)
CXXHFW    void fgetabr(double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void getabr(double &ar, double &br)
CXXWRP  {
CXXWRP    fgetabr(&ar,&br);
CXXWRP  }
C-----------------------------------------------------------------------

C     ========================
      subroutine getabr(ar,br)
C     ========================

C--   Get current values of Q2charm, Q2bottom and Q2top.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*80 subnam
      data subnam /'GETABR ( AR, BR )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Do the work
      ar = aar6
      br = bbr6

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int nfrmiq(int iset, int iq, int &ithresh);
C-----------------------------------------------------------------------
CXXHFW  #define fnfrmiq FC_FUNC(nfrmiq,NFRMIQ)
CXXHFW    int fnfrmiq(int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int nfrmiq(int iset, int iq, int &ithresh)
CXXWRP  {
CXXWRP    return fnfrmiq(&iset,&iq,&ithresh);
CXXWRP  }
C-----------------------------------------------------------------------

C     ======================================
      integer function nfrmiq(jset,iq,ithrs)
C     ======================================

C--   Returns number of flavors

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
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
      data subnam  /'NFRMIQ ( ISET, IQ, ITHRESH )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check jset in range
      call sqcIlele(subnam,'ISET',0,jset,mset0,'ISET does not exist')
C--   Check jset exists and is filled
      if(.not.Lfill7(jset)) call sqcSetMsg(subnam,'ISET',jset)
C--   Point to correct set
      call sparParTo5(ikeyf7(jset))

      if(abs(iq).lt.1 .or. abs(iq).gt.ntt2) then
C--     iq out of range
        nfrmiq = 0
        ithrs  = 0
      else
C--     iq in range, set ithrs flag
        ithrs  = 0
        it     = iq
        iz     = izfit5( it)
        nfrmiq = itfiz5(-iz)
        if(it.gt.0) then
          if( it.eq.itchm2 .or. it.eq.itbot2 .or. it.eq.ittop2) ithrs=1
        elseif(it.lt.0) then
          if(-it.eq.itchm2 .or.-it.eq.itbot2 .or.-it.eq.ittop2) ithrs=-1
        else
          stop 'NFRMIQ: encounter it = 0'
        endif
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int nflavs(int iq, int &ithresh);
C-----------------------------------------------------------------------
CXXHFW  #define fnflavs FC_FUNC(nflavs,NFLAVS)
CXXHFW    int fnflavs(int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int nflavs(int iq, int &ithresh)
CXXWRP  {
CXXWRP    return fnflavs(&iq,&ithresh);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      integer function nflavs(iq,ithrs)
C     =================================

C--   Returns number of flavors (routine kept 4 backward compatibility)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
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
      data subnam  /'NFLAVS ( IQ, ITHRESH )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Point to base set
      call sparParTo5(1)

      if(abs(iq).lt.1 .or. abs(iq).gt.ntt2) then
C--     iq out of range
        nflavs = 0
        ithrs  = 0
      else
C--     iq in range, set ithrs flag
        ithrs  = 0
        it     = iq
        iz     = izfit5( it)
        nflavs = itfiz5(-iz)
        if(it.gt.0) then
          if( it.eq.itchm2 .or. it.eq.itbot2 .or. it.eq.ittop2) ithrs=1
        elseif(it.lt.0) then
          if(-it.eq.itchm2 .or.-it.eq.itbot2 .or.-it.eq.ittop2) ithrs=-1
        else
          stop 'NFLAVS: encounter it = 0'
        endif
      endif

      return
      end

C=======================================================================
C==   Routines to manage evolution parameters ==========================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    void cpypar(double *array, int n, int iset);
C-----------------------------------------------------------------------
CXXHFW  #define fcpypar FC_FUNC(cpypar,CPYPAR)
CXXHFW    void fcpypar(double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void cpypar(double *array, int n, int iset)
CXXWRP  {
CXXWRP    fcpypar(array,&n,&iset);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===============================
      subroutine CpyPar(array,n,jset)
C     ===============================

C--   Copy parameter list to a local array

C--   jset  (in) : pdf set id 0=base, 1=unpol, 2=pol, etc.
C--   pars (out) : array with parameter values
C--   n     (in) : dimension of pars declared in the calling routine

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
      data subnam /'CPYPAR ( ARRAY, N, ISET )'/

      dimension array(*)

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check jset in range
      call sqcIlele(subnam,'ISET',0,jset,mset0,' ')
C--   Check n in range
      call sqcIlele(subnam,'N',13,n,9999,' ')

C--   Do the work
      if(jset.eq.0) then
        call sparListPar(1,array,ierr)
      elseif(Lfill7(jset)) then
        call sparListPar(ikeyf7(jset),array,ierr)
        array(13) = int(dparGetPar(stor7,isetf7(jset),idievtyp8))
      else
        ierr = 1
      endif

      if(ierr.ne.0) then
C--     Jset does not exist or has no parameters
        if(ierr.eq.1) then
        write(6,*) 'slot',ikeyf7(jset),'does not exist'
        else
        write(6,*) 'slot',ikeyf7(jset),'is empty'
        endif
        call sqcSetMsg(subnam,'ISET',jset)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void usepar(int iset);
C-----------------------------------------------------------------------
CXXHFW  #define fusepar FC_FUNC(usepar,USEPAR)
CXXHFW    void fusepar(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void usepar(int iset)
CXXWRP  {
CXXWRP    fusepar(&iset);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================
      subroutine UsePar(jset)
C     =======================

C--   Copy parameters of jset back into qpars6 and remake base set

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
      data subnam /'USEPAR ( ISET )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check jset in range
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')
C--   Jset does not exist or is not filled
      if(.not.Lfill7(jset)) then
        call sqcSetMsg(subnam,'ISET',jset)
      endif

C--   Reset /qpars6/ and remake base
      key   = ikeyf7(jset)
      call sparRemakeBase(key)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int keypar(int iset);
C-----------------------------------------------------------------------
CXXHFW  #define fkeypar FC_FUNC(keypar,KEYPAR)
CXXHFW    int fkeypar(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int keypar(int iset)
CXXWRP  {
CXXWRP    return fkeypar(&iset);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      integer function KeyPar(jset)
C     =============================

C--   Returns the version number of a parameter set

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'KEYPAR ( ISET )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check jset in range
      call sqcIlele(subnam,'ISET',0,jset,mset0,' ')

      if(jset.ne.0) then
C--     Version number of pdf set
        if(.not.Lfill7(jset)) then
          call sqcSetMsg(subnam,'ISET',jset)
          KeyPar = 0
        else
          key    = ikeyf7(jset)
          KeyPar = iparGetGroupKey(pars8,key,6)
        endif
      else
        KeyPar = iparGetGroupKey(pars8,1,6)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int keygrp(int iset, int igrp);
C-----------------------------------------------------------------------
CXXHFW  #define fkeygrp FC_FUNC(keygrp,KEYGRP)
CXXHFW    int fkeygrp(int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int keygrp(int iset, int igrp)
CXXWRP  {
CXXWRP    return fkeygrp(&iset,&igrp);
CXXWRP  }
C-----------------------------------------------------------------------

C     ====================================
      integer function KeyGrp(jset,igroup)
C     ====================================

C--   Returns the version number of a parameter selection

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'KEYGRP ( ISET, IGROUP )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check jset in range
      call sqcIlele(subnam,'ISET',0,jset,mset0,' ')
C--   Check igroup in range
      call sqcIlele(subnam,'IGROUP',1,igroup,6,
     + '1=order, 2=alfa, 3=thresholds, 4=scale, 5=cuts, 6=all')

      KeyGrp = 0

      if(jset.ne.0) then
C--     Version number of pdf set
        if(.not.Lfill7(jset)) then
          call sqcSetMsg(subnam,'ISET',jset)
        else
          key    = ikeyf7(jset)
          KeyGrp = iparGetGroupKey(pars8,key,igroup)
        endif
      else
        KeyGrp = iparGetGroupKey(pars8,1,igroup)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void pushcp();
C-----------------------------------------------------------------------
CXXHFW  #define fpushcp FC_FUNC(pushcp,PUSHCP)
CXXHFW    void fpushcp();
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void pushcp()
CXXWRP  {
CXXWRP    fpushcp();
CXXWRP  }
C-----------------------------------------------------------------------

C     =================
      subroutine Pushcp
C     =================

C--   Save base to LIFO buffer

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'PUSHCP'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Push base into  buffer
      call sparBufBase(1,ierr)
C--   Error
      if(ierr.eq.1) call sqcErrMsg(subnam,
     + 'LIFO buffer full: please call PULLCP first')

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void pullcp();
C-----------------------------------------------------------------------
CXXHFW  #define fpullcp FC_FUNC(pullcp,PULLCP)
CXXHFW    void fpullcp();
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void pullcp()
CXXWRP  {
CXXWRP    fpullcp();
CXXWRP  }
C-----------------------------------------------------------------------

C     =================
      subroutine Pullcp
C     =================

C--   Restore base from LIFO buffer

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
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
      data subnam /'PULLCP'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Pull base from buffer
      call sparBufBase(-1,ierr)
C--   Error
      if(ierr.eq.-1 .or. ierr.eq.2) call sqcErrMsg(subnam,
     + 'LIFO buffer empty: please call PUSHCP first')
C--   Restore /qpars6/ and remake base
      call sparRemakeBase(1)

      return
      end
