
C--   This is the file usrini.f containing initialization routines

C--   subroutine qcinitCPP(lun,fname,ls)
C--   subroutine qcinit(lun,fname)
C--   subroutine setlunCPP(lun,fname,ls)
C--   subroutine setlun(lun,fname)
C--   integer function nxtLun(lmin)
C--
C--   subroutine setvalCPP(chopt,ls,dval)
C--   subroutine setval(chopt,dval)
C--   subroutine getvalCPP(chopt,ls,dval)
C--   subroutine getval(chopt,dval)
C--   subroutine setintCPP(chopt,ls,ival)
C--   subroutine setint(chopt,ival)
C--   subroutine sqcSetNopt(ival)
C--   subroutine getintCPP(chopt,ls,ival)
C--   subroutine getint(chopt,ival)
C--
C--   subroutine qstoreCPP(chopt,ls,ival,dval)
C--   subroutine qstore(chopt,i,val)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM init routines from usrini.f                      */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM init routines from usrini.f                        */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM init routines from usrini.f                        */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C==   ===============================================================
C==   Initialization ================================================
C==   ===============================================================

C-----------------------------------------------------------------------
CXXHDR    void qcinit(int lun, string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fqcinitcpp FC_FUNC(qcinitcpp,QCINITCPP)
CXXHFW    void fqcinitcpp(int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void qcinit(int lun, string fname)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fqcinitcpp(&lun,cfname,&ls);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================
      subroutine qcinitCPP(lun,fname,ls)
C     ==================================

C--   Initialize qcdnum Cxx interface

      implicit double precision (a-h,o-z)

      character*(100) fname

      if(ls.gt.100) stop 'qcinitCPP: input file name > 100 characters'

      call qcinit(lun,fname(1:ls))

      return
      end

C     ============================
      subroutine qcinit(lun,fname)
C     ============================

C--   Initialize qcdnum

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qvers1.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'
      include 'qpars6.inc'

      character*(*) fname

C--   Error; -6 is allowed for std output without banner
      if(lun.le.0.and.lun.ne.-6) goto 500

C--   Version   12345678 
      ivers1 =  180000
      cvers1 = '18-00-00'
      cdate1 = '08-03-22'
      
C--   Initialize qcdnum status
      do j = 1,mset0
        do i = 1,mbp0
          istat4(i,j) = 0
        enddo
      enddo  

C--   Set initialization flag
      iniflg4 = 123456
      
C--   Assign status bits
      call sqcBitIni

C--   Initialize constants
      call sqcIniCns

C--   Setup transformation from udscbt basis to si/ns basis
      call sqcPdfMat
      call sqcPdfMatn
      
C--   Initialize weight tables and pdf sets      
      call sqcIniWt

C--   Set a few status bits
      do i = 1,mset0
        call sqcSetbit(ibinit4,istat4(1,i),mbp0)   !initialization done
      enddo  

C--   Set output stream
      call sqcSetLun(abs(lun),fname)
C--   Print banner; no banner if lun = -6
      if(lun.ne.-6) call sqcBanner(lunerr1)
C--   Print please refer to ...
      call sqcReftoo(lunerr1)

      return

 500  continue

C--   Error message
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in QCINIT ( LUN, FNAME ) ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) 'LUN = ',lun,' should be positive'

      stop
      end

C-----------------------------------------------------------------------
CXXHDR    void setlun(int lun, string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fsetluncpp FC_FUNC(setluncpp,SETLUNCPP)
CXXHFW    void fsetluncpp(int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setlun(int lun, string fname)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fsetluncpp(&lun,cfname,&ls);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================
      subroutine setlunCPP(lun,fname,ls)
C     ==================================

      implicit double precision (a-h,o-z)

      character*(100) fname

      if(ls.gt.100) stop 'setlunCPP: input file name > 100 characters'

      call setlun(lun,fname(1:ls))

      return
      end
      
C     ============================
      subroutine setlun(lun,fname)
C     ============================

C--   (Re)Set logical unit number and output file name.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*(*) fname

      character*80 subnam
      data subnam /'SETLUN ( LUN, FNAME )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif
C--   Check if input is in allowed range
      call sqcIlele(subnam,'LUN',1,lun,99,
     +             'LUN should be between 1 and 99')
C--   Check that fname nonempty for lun .ne. 6
      if(lun.ne.6 .and. imb_lenoc(fname).eq.0) then
        call sqcErrMsg(subnam,'FNAME is empty')
      endif

C--   Do the work
      call sqcSetLun(lun,fname)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int nxtlun(int lmin);
C-----------------------------------------------------------------------
CXXHFW  #define fnxtlun FC_FUNC(nxtlun,NXTLUN)
CXXHFW    int fnxtlun(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int nxtlun(int lmin)
CXXWRP  {
CXXWRP    return fnxtlun(&lmin);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      integer function nxtLun(lmin)
C     =============================

C--   Find first free lun in the range [max(lmin,10)-99]
C--   Returns zero if no free lun found

      nxtLun = iqcLunFree(lmin)

      return
      end

C==   ==================================================================
C==   Routines to set and get QCDNUM parameters ========================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void setval(string opt, double val);
C-----------------------------------------------------------------------
CXXHFW  #define fsetvalcpp FC_FUNC(setvalcpp,SETVALCPP)
CXXHFW    void fsetvalcpp(char*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setval(string opt, double val)
CXXWRP  {
CXXWRP    int ls = opt.size();
CXXWRP    char *copt = new char[ls+1];
CXXWRP    strcpy(copt,opt.c_str());
CXXWRP    fsetvalcpp(copt,&ls,&val);
CXXWRP    delete[] copt;
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine setvalCPP(chopt,ls,dval)
C     ===================================

      implicit double precision (a-h,o-z)

      character*(100) chopt

      if(ls.gt.100) stop 'setvalCPP: input CHOPT size > 100 characters'

      call setval(chopt(1:ls),dval)

      return
      end

C     =============================
      subroutine setval(chopt,dval)
C     =============================

C--   Set double precision value in /qpars6/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*(*) chopt
      character*4   opt

      character*80 subnam
      data subnam /'SETVAL ( CHOPT, DVAL )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Use mbutil for character string manipulations
      len        = min(imb_lenoc(chopt),4)
C--   Avoid changing the argument so make a copy 
      opt(1:len) = chopt(1:len)
C--   OK now convert to upper case
      call smb_cltou(opt)
C--   Do the work
      if    (opt(1:len).eq.'EPSI') then
        call sqcDlele(subnam,'EPSI',1.D-10,dval,1.D-4,' ')
        aepsi6 = dval
      elseif(opt(1:len).eq.'EPSG') then
        call sqcDlele(subnam,'EPSG',1.D-9,dval,1.D-1,' ')
        gepsi6 = dval
      elseif(opt(1:len).eq.'ELIM') then
        call sqcDlele(subnam,'ELIM',-1.D10,dval,1.D10,' ')
        dflim6 = dval  
      elseif(opt(1:len).eq.'ALIM') then
        call sqcDlele(subnam,'ALIM',1.D-10,dval,1.D10,' ')
        aslim6 = dval
      elseif(opt(1:len).eq.'QMIN') then
        call sqcDlele(subnam,'QMIN',1.D-1,dval,qlimu6,' ')
        qlimd6 = dval
      elseif(opt(1:len).eq.'QMAX') then
        call sqcDlele(subnam,'QMAX',qlimd6,dval,1.D11,' ')
        qlimu6 = dval
      elseif(opt(1:len).eq.'NULL') then
        qnull6 = dval
      else
        call sqcErrMsg(subnam,'CHOPT = '//chopt//' : unknown option')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void getval(string opt, double &val);
C-----------------------------------------------------------------------
CXXHFW  #define fgetvalcpp FC_FUNC(getvalcpp,GETVALCPP)
CXXHFW    void fgetvalcpp(char*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void getval(string opt, double &val)
CXXWRP  {
CXXWRP    int ls = opt.size();
CXXWRP    char *copt = new char[ls+1];
CXXWRP    strcpy(copt,opt.c_str());
CXXWRP    fgetvalcpp(copt,&ls,&val);
CXXWRP    delete[] copt;
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine getvalCPP(chopt,ls,dval)
C     ===================================

      implicit double precision (a-h,o-z)

      character*(100) chopt

      if(ls.gt.100) stop 'getvalCPP: input CHOPT size > 100 characters'

      call getval(chopt(1:ls),dval)

      return
      end

C     =============================
      subroutine getval(chopt,dval)
C     =============================

C--   Get double precision value from /qpars6/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*(*) chopt
      character*4   opt

      character*80 subnam
      data subnam /'GETVAL ( CHOPT, DVAL )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Use mbutil for character string manipulations
      len        = min(imb_lenoc(chopt),4)
C--   Avoid changing the argument so make a copy 
      opt(1:len) = chopt(1:len)
C--   OK now convert to upper case
      call smb_cltou(opt)
C--   Do the work
      if    (opt(1:len).eq.'EPSI') then
        dval = aepsi6
      elseif(opt(1:len).eq.'EPSG') then
        dval = gepsi6
      elseif(opt(1:len).eq.'ELIM') then
        dval = dflim6  
      elseif(opt(1:len).eq.'ALIM') then
        dval = aslim6
      elseif(opt(1:len).eq.'QMIN') then
        dval = qlimd6
      elseif(opt(1:len).eq.'QMAX') then
        dval = qlimu6
      elseif(opt(1:len).eq.'NULL') then
        dval = qnull6
      else
        call sqcErrMsg(subnam,'CHOPT = '//chopt//' : unknown option')
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void setint(string opt, int ival);
C-----------------------------------------------------------------------
CXXHFW  #define fsetintcpp FC_FUNC(setintcpp,SETINTCPP)
CXXHFW    void fsetintcpp(char*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setint(string opt, int ival)
CXXWRP  {
CXXWRP    int ls = opt.size();
CXXWRP    char *copt = new char[ls+1];
CXXWRP    strcpy(copt,opt.c_str());
CXXWRP    fsetintcpp(copt,&ls,&ival);
CXXWRP    delete[] copt;
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine setintCPP(chopt,ls,ival)
C     ===================================

      implicit double precision (a-h,o-z)

      character*(100) chopt

      if(ls.gt.100) stop 'setintCPP: input CHOPT size > 100 characters'

      call setint(chopt(1:ls),ival)

      return
      end

C     =============================
      subroutine setint(chopt,ival)
C     =============================

C--   Set integer value in /qpars6/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*(*) chopt
      character*4   opt

      character*80 subnam
      data subnam /'SETINT ( CHOPT, IVAL )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Use mbutil for character string manipulations
      len        = min(imb_lenoc(chopt),4)
C--   Avoid changing the argument so make a copy 
      opt(1:len) = chopt(1:len)
C--   OK now convert to upper case
      call smb_cltou(opt)
C--   Do the work
      if    (opt(1:len).eq.'ITER') then
        call sqcIlele(subnam,'ITER',-9999,ival,5,' ')
        niter6 = ival
      elseif(opt(1:len).eq.'TLMC') then
        itlmc6 = ival
      elseif(opt(1:len).eq.'NOPT') then
        call sqcSetNopt(ival)
      elseif(opt(1:len).eq.'EDBG') then
        idbug6 = ival
      else
        call sqcErrMsg(subnam,'CHOPT = '//chopt//' : unknown option')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     ===========================
      subroutine sqcSetNopt(ival)
C     ===========================

C--   Sets number of perturbative terms for evdglap
C--   Produces error messages if ival not ok

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      character*10 cval
      character*80 subr
      character*1  quot
      data quot /''''/
      character*35 msg(3)
C--                       1         2         3
C--              12345678901234567890123456789012345
      data msg /'NOPT cannot be .le. zero           ',
     +          'NOPT cannot have more than 3 digits',
     +          'All digits in NOPT must be nonzero '/

      character*1 digit(9)
      data digit /'1','2','3','4','5','6','7','8','9'/

C--   Convert to character string
      call smb_itoch(ival,cval,mord)

C--   Check ival not negative
      if(ival.le.0) then
        imsg = 1
        goto 500
      endif

C--   Get number of digits in ival
      if(mord.gt.mord0) then
        imsg = 2
        goto 500
      endif

C--   Get each digit in ival; reject zero digits
      do i = 1,mord
        if(cval(i:i) .eq. '0') then
          imsg = 3
          goto 500
        else
          do j = 1,9
            if(cval(i:i) .eq. digit(j)) nnopt6(i) = j
          enddo
        endif
      enddo

      nnopt6(0) = mord
      inopt6    = ival

      return

C--   Error message
 500  continue
      write(subr,'(''GETINT ( '',A,''NOPT'',A,'' , '',A,'' )'')')
     + quot,quot,cval(1:mord)
      call sqcErrMsg(subr,msg(imsg))
      return

      end

C-----------------------------------------------------------------------
CXXHDR    void getint(string opt, int &ival);
C-----------------------------------------------------------------------
CXXHFW  #define fgetintcpp FC_FUNC(getintcpp,GETINTCPP)
CXXHFW    void fgetintcpp(char*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void getint(string opt, int &ival)
CXXWRP  {
CXXWRP    int ls = opt.size();
CXXWRP    char *copt = new char[ls+1];
CXXWRP    strcpy(copt,opt.c_str());
CXXWRP    fgetintcpp(copt,&ls,&ival);
CXXWRP    delete[] copt;
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine getintCPP(chopt,ls,ival)
C     ===================================

      implicit double precision (a-h,o-z)

      character*(100) chopt

      if(ls.gt.100) stop 'getintCPP: input CHOPT size > 100 characters'

      call getint(chopt(1:ls),ival)

      return
      end

C     =============================
      subroutine getint(chopt,ival)
C     =============================

C--   Get integer value from /qpars6/, /qluns1/ or /qcdnum/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qvers1.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*(*) chopt
      character*4   opt

      character*80 subnam
      data subnam /'GETINT ( CHOPT, IVAL )'/
      
C--   Use mbutil for character string manipulations
      len        = min(imb_lenoc(chopt),4)
C--   Avoid changing the argument so make a copy 
      opt(1:len) = chopt(1:len)
C--   OK now convert to upper case
      call smb_cltou(opt)

C--   Set version number to zero if qcdnum is not initialised
      if(opt(1:len).eq.'VERS' .and. iniflg4.ne.123456) then
        ival = 0
        return
      endif  

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Do the work
      if    (opt(1:len).eq.'ITER') then
        ival = niter6
      elseif(opt(1:len).eq.'TLMC') then
        ival = itlmc6
      elseif(opt(1:len).eq.'NOPT') then
        ival = inopt6
      elseif(opt(1:len).eq.'EDBG') then
        ival = idbug6
      elseif(opt(1:len).eq.'LUNQ') then
        ival = lunerr1
      elseif(opt(1:len).eq.'MSET') then
        ival = mset0
      elseif(opt(1:len).eq.'MQS0') then
        ival = mqs0
      elseif(opt(1:len).eq.'MXG0') then
        ival = mxg0
      elseif(opt(1:len).eq.'MXX0') then   
        ival = mxx0
      elseif(opt(1:len).eq.'MQQ0') then   
        ival = mqq0
      elseif(opt(1:len).eq.'MKY0') then
        ival = mky0
      elseif(opt(1:len).eq.'MPT0') then   
        ival = mpt0
      elseif(opt(1:len).eq.'MST0') then
        ival = mst0
      elseif(opt(1:len).eq.'MBF0') then   
        ival = mbf0
      elseif(opt(1:len).eq.'MCE0') then
        ival = mce0
      elseif(opt(1:len).eq.'NWF0') then
        ival = nwf0
      elseif(opt(1:len).eq.'VERS') then   
        ival = ivers1
      else
        call sqcErrMsg(subnam,'CHOPT = '//chopt//' : unknown option')
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void qstore(string opt, int ival, double &val);
C-----------------------------------------------------------------------
CXXHFW  #define fqstorecpp FC_FUNC(qstorecpp,QSTORECPP)
CXXHFW    void fqstorecpp(char*, int*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void qstore(string opt, int ival, double &val)
CXXWRP  {
CXXWRP    int ls = opt.size();
CXXWRP    char *copt = new char[ls+1];
CXXWRP    strcpy(copt,opt.c_str());
CXXWRP    fqstorecpp(copt,&ls,&ival,&val);
CXXWRP    delete[] copt;
CXXWRP  }
C-----------------------------------------------------------------------

C     ========================================
      subroutine qstoreCPP(chopt,ls,ival,dval)
C     ========================================

      implicit double precision (a-h,o-z)

      character*(100) chopt

      if(ls.gt.100) stop 'qstoreCPP: input CHOPT size > 100 characters'

      call qstore(chopt(1:ls),ival,dval)

      return
      end

C     ==============================
      subroutine qstore(chopt,i,val)
C     ==============================

C--   Set and get values from qstore
C--
C--   chopt   (in)   R/W/L/U = read/write/lock/unlock
C--   i       (in)   index [1,msq0]
C--   val  (inout)   value
C--
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*(*) chopt

      character*80 subnam
      data subnam /'QSTORE ( ACTION, I, VAL )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Check user input
      call sqcIlele(subnam,'I',1,i,mqs0,' ')

C--   Do the work
      if    (chopt(1:1).eq.'R' .or. chopt(1:1).eq.'r') then
        val = qstore6(i)
      elseif(chopt(1:1).eq.'W' .or. chopt(1:1).eq.'w') then
        if(.not.Lqswrite6) then
          call sqcErrMsg(subnam,'QSTORE is locked, please unlock')
        endif
        qstore6(i) = val
      elseif(chopt(1:1).eq.'L' .or. chopt(1:1).eq.'l') then
        Lqswrite6 = .false.
      elseif(chopt(1:1).eq.'U' .or. chopt(1:1).eq.'u') then
        Lqswrite6 = .true.
      else
        call sqcErrMsg(subnam,
     +  'ACTION = '//chopt//' : first character should be R,W,L,U')
      endif

      return
      end

