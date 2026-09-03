
C--  This is the file istore.f with istore routines
C--  This file also contains C++ wrappers for all the Fortran routines

C--   sws_IwInit(iw,nw)             Create istore
C--   swsIwEbuf(iw,txt,opt)         Handle message text buffer
C--   swsIwEmsg(iw,n,srname)        Out-of-space error message
C--   sws_SetIwn(iw,nw)             Update istore size info
C--   iws_Iarray(iw,imi,ima)        Create new array(imi:ima)
C--   iws_IAread(iw,iarr,n)         Create array(1:n) and read iarr
C--   iws_DAread(iw,darr,n)         Create array(1:n) and read int(darr)
C--   sws_IwWipe(iw,ia)             Wipe istore
C--
C--   iws_IhSize()                  Return header size
C--   iws_IsaIstore(iw)             Yes/no iw is an istore
C--   iws_IwSize(iw)                Total size of iw
C--   iws_IwNused(iw)               Words used in iw
C--   iws_IwNarrays(iw)             Number of arrays in iw
C--   iws_IwNheader(iw)             Header size
C--   iws_IwObjectType(iw,ia)       Object type
C--   iws_IwAsize(iw,ia)            Array size
C--   iws_IwAnumber(iw,ia)          Array number
C--   iws_IwAfprint(iw,ia)          Fingerprint of root or array
C--   iws_IwAdim(iw,ia)             Array dimension
C--   iws_IwAimin(iw,ia)            Imin of array
C--   iws_IwAimax(iw,ia)            Imax of array
C--   iws_IaAbegin(iw,ia)           Ia begin body
C--   iws_IaAend(iw,ia)             Ia end body
C--   iws_IwKnul(iw,ia)             K0 pointer coefficient
C--   iws_ArrayI(iw,ia,i)           Ia of element A(i)

C--   sws_IwTree(iw)                Print istore tree
C--   swsIWtree(iw,iroot)           Print istore tree
C--   swsiwprnt(iw,ia,iroot)        Print istore summary line
C--   swsaprint(iw,ia,iroot)        Print array  summary line
C--   sws_IwHead(iw,ia)             Print header


C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  WSTORE istore routines from istore.f                    */
CXXHDR    /************************************************************/
CXXHDR
CXXHDR    // alrady defined in file wstore.f
CXXHDR    // inline int iaFtoC(int ia) { return ia-1; };
CXXHDR    // inline int iaCtoF(int ia) { return ia+1; };
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  WSTORE istore routines from istore.f                      */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  WSTORE istore routines from istore.f                      */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    void sws_IwInit(int *iw, int nw, string txt);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_iwinitcpp FC_FUNC(sws_iwinitcpp,SWS_IWINITCPP)
CXXHFW    void fsws_iwinitcpp(int*, int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_IwInit(int *iw, int nw, string txt)
CXXWRP    {
CXXWRP      int ls = txt.size();
CXXWRP      char *ctxt = new char[ls+1];
CXXWRP      strcpy(ctxt,txt.c_str());
CXXWRP      fsws_iwinitcpp(iw, &nw, ctxt, &ls);
CXXWRP      delete[] ctxt;
CXXWRP    }
C-----------------------------------------------------------------------

C     ======================================
      subroutine sws_IwInitCPP(iw,nw,txt,ls)
C     ======================================

      implicit double precision (a-h,o-z)

      dimension iw(*)
      character*(100) txt

      if(ls.gt.100) stop
     +             'WSTORE::SWS_IWINIT: input text > 100 characters'

      call sws_IwInit(iw,nw,txt(1:ls))

      return
      end

C     ================================
      subroutine sws_IwInit(iw,nw,txt)
C     ================================

C--   Create new istore
C--
C--   iw  (in): integer array iw(nw)
C--   nw  (in): dimension of iw declared in calling routine

C--   Author: Michiel Botje h24@nikhef.nl   30-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      character*(*) txt
      character*10  date, time, zone
      dimension     ival(8)
      character*20  ntxt

      save icnt
      data icnt /0/

C--   Check input
      if(nw.le.0) stop
     + 'WSTORE:SWS_IWINIT: cannot have istore size NW <= 0'
C--   Check size
      nhead  = nwHeader1
      nsize  = nhead
      nneed  = nsize + 1
      if(nw.lt.nneed) then
        call smb_itoch(nneed,ntxt,leng)
        write(6,*)
     + 'WSTORE:IWS_IWINIT: workspace size must be at least ',
     +  ntxt(1:leng),' words'
        if(imb_lastc(txt).ne.0) write(6,*) txt
        stop
      endif
C--   Istore fingerprint
      call date_and_time(date,time,zone,ival)      !millisec resolution
      icnt  = icnt+1                               !count calls
      iseed = 0                                    !call-dependent seed
      do i = 1,4
        call smb_cbyte(mod(icnt+i,256),1,iseed,i)  !set 4 bytes of iseed
      enddo
      ihash = imb_ihash(iseed,ival,8)              !hash year-date-time
C--   Initialise istore
      call smb_Ifill(iw,nw,0)
C--   Workspace header
      iw(ICword1+iroot0) = iWsVersion0             !Control word
      iw(IwAddr1+iroot0) = 0                       !IW address
      iw(NFTabl1+iroot0) = 0                       !Fskip table
      iw(NBTabl1+iroot0) = 0                       !Bskip table
      iw(Ifprnt1+iroot0) = ihash                   !Fingerprint
      iw(INwUse1+iroot0) = nsize                   !Words used
      iw(NObjec1+iroot0) = 0                       !Number of objects
      iw(IWlstT1+iroot0) = nhead                   !IW last table in w
      iw(INwMax1+iroot0) = nw                      !Istore total size
      iw(INhead1+iroot0) = nhead                   !Header size
C--   Store message
      call swsIwEbuf(iw,txt,'in')

      return
      end

C     ================================
      subroutine swsIwEbuf(iw,txt,opt)
C     ================================

C--   Store (opt = 'in') or retrieve (opt = 'out') message text

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)
      character*(*) txt, opt

      save ebuf, ifp, nebuf, first
      character*80 ebuf(mebuf0)
      dimension ifp(mebuf0)
      logical first
      data first/.true./

C--   Initialise
      if(first) then
        nebuf = 0
        do i = 1,mebuf0
          call smb_cfill(' ',ebuf(i))
          ifp(i) = 0
        enddo
        first = .false.
      endif

C--   Workspace fingerprint
      ifpw = int(iw(IFprnt1+iroot0))
C--   Find workspace entry
      iws  = 0
      do i = 1,nebuf
        if(ifp(i).eq.ifpw) iws = i
      enddo
      if(opt(1:1).eq.'i' .or. opt(1:1).eq.'I') then
C--     Store new message
C--     Empty message
        if(imb_lastc(txt).eq.0) return
        if(iws.ne.0) then
C--       Overwrite
          ebuf(iws) = txt
        else
C--       New entry
          nebuf = nebuf+1
          if(nebuf.gt.mebuf0) then
            write(6,*) 'WSTORE:SWS_IWINIT: message buffer size exceeded'
            write(6,*)
     +                'Please increase MEBUF0 in mbutil/inc/wstore0.inc'
            stop
          endif
          ebuf(nebuf) = txt
          ifp(nebuf)  = ifpw
        endif
      elseif(opt(1:1).eq.'o' .or. opt(1:1).eq.'O') then
C--     Retrieve message
        if(iws.eq.0) then
C--       Workspace not found
          call smb_cfill(' ',txt)
        else
          txt = ebuf(iws)
        endif
      else
        stop 'WSTORE:swsIwEbuf: unknown option'
      endif

      return
      end

C     =================================
      subroutine swsIwEmsg(iw,n,srname)
C     =================================

C--   Out-of-space message

      implicit double precision(a-h,o-z)

      dimension iw(*)
      character*(*) srname
      character*20  ntxt
      character*80  txt

      i1 = imb_frstc(srname)
      i2 = imb_lastc(srname)
      call smb_itoch(n,ntxt,leng)
      call swsIwEbuf(iw,txt,'out')
      write(6,*) srname(i1:i2),': workspace size must be at least ',
     +  ntxt(1:leng),' words'
      if(imb_lastc(txt).ne.0) write(6,*) txt
      stop

      end

C-----------------------------------------------------------------------
CXXHDR    void sws_SetIwN(int *iw, int nw);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_setiwn FC_FUNC(sws_setiwn,SWS_SETIWN)
CXXHFW    void fsws_setiwn(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_SetIwN(int *iw, int nw)
CXXWRP    {
CXXWRP     fsws_setiwn(iw, &nw);
CXXWRP    }
C-----------------------------------------------------------------------

C     ============================
      subroutine sws_SetIwN(iw,nw)
C     ============================

C--   Update size information in the istore header
C--
C--   iw              (in): istore
C--   nw              (in): actual istore size

C--   Author: Michiel Botje h24@nikhef.nl   31-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

C--   Check if workspace
      if(iw(iCword1+iroot0).ne.iWsVersion0) then
        stop 'WSTORE:SWS_SETIWN: IW is not an istore'
      endif
C--   Check input
      nwold = iw(INwMax1+iroot0)
      if(nw.lt.nwold) stop
     +               'WSTORE:SWS_SETIWN: cannot decrease istore size'
C--   Update header
      iw(INwMax1+iroot0) = nw                    !Istore total size

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_Iarray(int* iw, int imi, int ima);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iarray FC_FUNC(iws_iarray,IWS_IARRAY)
CXXHFW    int fiws_iarray(int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_Iarray(int* iw, int imi, int ima)
CXXWRP    {
CXXWRP     int ia = fiws_iarray(iw, &imi, &ima);
CXXWRP     return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     =======================================
      integer function iws_Iarray(iw,imi,ima)
C     =======================================

C--   Book new array
C--
C--   iw          (in): istore
C--   imi,ima     (in): index limits
C--   iws_Iarray (out): IA address of the new array

C--   Author: Michiel Botje h24@nikhef.nl   31-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*), karr(0:1), imin(1), imax(1), ival(1)

C--   Check if workspace
      if(iw(iCword1+iroot0).ne.iWsVersion0) then
        stop 'WSTORE:IWS_IARRAY: IW is not an istore'
      endif
C--   Check input and get array size
      if(imi.gt.ima) stop 'WSTORE:IWS_IARRAY: imin > imax'
      nw = ima-imi+1
C--   Addresses
      iaR    = iroot0                             !root address
      iaL    = iw(IWlstT1+iaR) + iaR              !IA current table
      iaT    = iw(INwUse1+iaR) + iaR              !IA new table
      NBskip = iaL-iaT                            !Bskip table
C--   Check workspace size
      NHskip = iw(INhead1+iaR)                    !header size
      NWused = iw(INwUse1+iaR)                    !words used
      NTsize = nw + NHskip                        !table size
      NWneed = NWused + NTsize + 1                !words needed
      Nspace = iw(INwMax1+iaR)                    !total size
C--   Run out-of-space
      if(NWneed.gt.Nspace) call swsIwEmsg(iw,NWneed,'WSTORE:IWS_IARRAY')
C--   Initialise
      do i = iaT,NWneed
        iw(i) = 0
      enddo
C--   Array definition
      imin(1) = imi
      imax(1) = ima
      ittb1   = NHskip
      call smb_bkmat(imin,imax,karr,1,ittb1,ittb2)
      if(ittb2.ne.NTsize-1) stop
     +                     'WSTORE:IWS_IARRAY: problem with table size'
C--   Array fingerprint
      ival(1) = 1
      ihash   = imb_ihash(0    ,ival,1)
      ihash   = imb_ihash(ihash,karr,2)
      ihash   = imb_ihash(ihash,imin,1)
      ihash   = imb_ihash(ihash,imax,1)
C--   Number of arrays in iw
      NObjec  = iw(NObjec1+iaR) + 1               !number of arrays
C--   Fill array header
      iw(ICword1+iaT) = iCWTable0                 !control word
      iw(IwAddr1+iaT) = NWused                    !IW
      iw(NFTabl1+iaT) = 0                         !Fskip array
      iw(NBTabl1+iaT) = NBskip                    !Bskip array
      iw(Ifprnt1+iaT) = ihash                     !fingerprint
      iw(INwUse1+iaT) = NTsize                    !array size
      iw(IObjec1+iaT) = NObjec                    !array number
      iw(KnulTb1+iaT) = karr(0)                   !k0 value
      iw(IminTb1+iaT) = imi                       !Imin
      iw(ImaxTb1+iaT) = ima                       !Imax
      iw(ItBody1+iaT) = ittb1                     !IT begin body
      iw(ItBend1+iaT) = ittb2                     !IT end body
C--   Update istore header
      iw(NFTabl1+iaR) = NHskip                    !Fskip
      iw(INwUse1+iaR) = NWused+NTsize             !Words used
      iw(NObjec1+iaR) = NObjec                    !number of arrays
      iw(IWlstT1+iaR) = NWused                    !IW this table
C--   Return address
      iws_Iarray = iaT
C--   Done for first array
      if(NBskip.eq.0) return
C--   Update Header previous array
      iw(NFTabl1+iaL) = iaT - iaL                 !Fskip table

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IAread(int* iw, int* iarr, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iaread FC_FUNC(iws_iaread,IWS_IAREAD)
CXXHFW    int fiws_iaread(int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IAread(int* iw, int* iarr, int n)
CXXWRP    {
CXXWRP     int ia = fiws_iaread(iw, iarr, &n);
CXXWRP     return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ======================================
      integer function iws_IAread(iw,iarr,n)
C     ======================================

C--   Create array(1:n) and read iarr
C--
C--   iw          (in): istore
C--   iarr        (in): integer array to read in
C--   iws_IAread (out): IA address of the new array(1:n)

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

C--   Check if workspace
      if(iw(iCword1+iroot0).ne.iWsVersion0) then
        stop 'WSTORE:IWS_IAREAD: IW is not an istore'
      endif
C--   Check input
      if(n.le.0) stop 'WSTORE:IWS_IAREAD: n must be > 0'
C--   Check workspace size
      iaR    = iroot0                             !root address
      NHskip = iw(INhead1+iaR)                    !header size
      NWused = iw(INwUse1+iaR)                    !words used
      NTsize = n + NHskip                         !table size
      NWneed = NWused + NTsize + 1                !words needed
      Nspace = iw(INwMax1+iaR)                    !total size
C--   Run out-of-space
      if(NWneed.gt.Nspace) call swsIwEmsg(iw,NWneed,'WSTORE:IWS_IAREAD')
C--   Book new array
      iaT = iws_Iarray(iw,1,n)
C--   Copy
      ia1 = iws_IaAbegin(iw,iaT)
      call smb_Icopy(iarr,iw(ia1),n)
C--   Done ...
      iws_IAread = iaT

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_DAread(int* iw, double* darr, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_daread FC_FUNC(iws_daread,IWS_DAREAD)
CXXHFW    int fiws_daread(int*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_DAread(int* iw, double* darr, int n)
CXXWRP    {
CXXWRP     int ia = fiws_daread(iw, darr, &n);
CXXWRP     return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ======================================
      integer function iws_DAread(iw,darr,n)
C     ======================================

C--   Create array(1:n) and read int(darr)
C--
C--   iw          (in): istore
C--   darr        (in): double precision array to read in
C--   iws_DAread (out): IA address of the new array(1:n)

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

C--   Check if workspace
      if(iw(iCword1+iroot0).ne.iWsVersion0) then
        stop 'WSTORE:IWS_DAREAD: IW is not an istore'
      endif
C--   Check input
      if(n.le.0) stop 'WSTORE:IWS_DAREAD: n must be > 0'
C--   Check workspace size
      iaR    = iroot0                             !root address
      NHskip = iw(INhead1+iaR)                    !header size
      NWused = iw(INwUse1+iaR)                    !words used
      NTsize = n + NHskip                         !table size
      NWneed = NWused + NTsize + 1                !words needed
      Nspace = iw(INwMax1+iaR)                    !total size
C--   Run out-of-space
      if(NWneed.gt.Nspace) call swsIwEmsg(iw,NWneed,'WSTORE:IWS_DAREAD')
C--   Book new array
      iaT = iws_Iarray(iw,1,n)
C--   Copy
      ia1 = iws_IaAbegin(iw,iaT)
      call smb_VDtoI(darr,iw(ia1),n)
C--   Done ...
      iws_DAread = iaT

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void sws_IwWipe(int* obj);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_iwwipe FC_FUNC(sws_iwwipe,SWS_IWWIPE)
CXXHFW    void fsws_iwwipe(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_IwWipe(int* obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     fsws_iwwipe(iw, &ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ============================
      subroutine sws_IwWipe(iw,ia)
C     ============================

C--   Wipe iw starting at ia

C--   Author: Michiel Botje h24@nikhef.nl   02-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

C--   Check input
      if(iw(iroot0).ne.iWsVersion0) stop
     +                        'WSTORE:SWS_IWWIPE: IW is not an istore'
      if(ia.le.0 .or. ia.gt.iw(INwUse1+iroot0)) stop
     +                        'WSTORE:SWS_IWWIPE: IA out of range'

      iaR    = iroot0                              !Root address
      iaT    = iw(NFTabl1+iaR) + iaR               !First array address
      narrs  = iw(NObjec1+iaR)                     !Number of arrays
      nwuse  = iw(INwUse1+iaR)                     !Words used
      nwmax  = iw(INwMax1+iaR)                     !Total size
      nhead  = iw(INhead1+iaR)                     !Header size

      if(narrs.eq.0) then
C--     empty istore
        return
      elseif(ia.eq.iaR .or. ia.eq.iaT) then
C--     ia = root or first array
        nsize           = nhead                    !Size after wipe
        iw(NFTabl1+iaR) = 0                        !Fskip table
        iw(INwUse1+iaR) = nsize                    !Words used
        iw(NObjec1+iaR) = 0                        !Number of objects
        iw(IWlstT1+iaR) = nhead                    !IW last table in w
        call smb_Ifill(iw(nsize+1),nwmax-nsize,0)  !Wipe

      elseif(iw(iCword1+ia).ne.iCWTable0) then
C--     ia is not an array address
        stop 'WSTORE:SWS_IWWIPE: IA is not an array address'

      else
C--     ia is array address but not of first array
        iaT   = iw(NBTabl1+ia) + ia                !Last arr after wipe
        nsize = ia-1                               !Size after wipe
        narrs = iw(IObjec1+iaT)                    !Nb arrays after wipe
        iw(NFTabl1+iaT) = 0                        !Forward skip
        iw(INwUse1+iaR) = nsize                    !Words used
        iw(NObjec1+iaR) = narrs                    !Number of objects
        iw(IWlstT1+iaR) = iaT-iaR                  !IW last table in w
        call smb_Ifill(iw(nsize+1),nwmax-nsize,0)  !Wipe

      endif

      return
      end

C=======================================================================
C==   Access to information ============================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    int iws_IhSize();
C-----------------------------------------------------------------------
CXXHFW  #define fiws_ihsize FC_FUNC(iws_ihsize,IWS_IHSIZE)
CXXHFW    int fiws_ihsize();
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IhSize()
CXXWRP    {
CXXWRP     return fiws_ihsize();
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================
      integer function iws_IhSize()
C     =============================

C--   Return header size

C--   Author: Michiel Botje h24@nikhef.nl   31-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      iws_IhSize = nwHeader0

      return
      end

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    int iws_IsaIstore(int *obj);
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  #define fiws_isaistore FC_FUNC(iws_isaistore,IWS_ISAISTORE)
CXXHFW    int fiws_isaistore(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IsaIstore(int *obj)
CXXWRP    {
CXXWRP     int  ia = *(obj+1);
CXXWRP     int *iw = obj-ia;
CXXWRP     return fiws_isaistore(iw);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function iws_IsaIstore(iw)
C     ==================================

C--   Return 1 (0) if iw is (not) an istore

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      if(iw(iCword0+iroot0).eq.iWsVersion0) then
        iws_IsaIstore = 1
      else
        iws_IsaIstore = 0
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwSize(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwsize FC_FUNC(iws_iwsize,IWS_IWSIZE)
CXXHFW    int fiws_iwsize(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwSize(int *obj)
CXXWRP    {
CXXWRP     int  ia = *(obj+1);
CXXWRP     int *iw = obj-ia;
CXXWRP     return fiws_iwsize(iw);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================
      integer function iws_IwSize(iw)
C     ===============================

C--   Total size of iw

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwSize = iw(iNwMax1+iroot0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwNused(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwnused FC_FUNC(iws_iwnused,IWS_IWNUSED)
CXXHFW    int fiws_iwnused(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwNused(int *obj)
CXXWRP    {
CXXWRP     int  ia = *(obj+1);
CXXWRP     int *iw = obj-ia;
CXXWRP     return fiws_iwnused(iw);
CXXWRP    }
C-----------------------------------------------------------------------

C     ================================
      integer function iws_IwNused(iw)
C     ================================

C--   Words used in iw (without trailer)

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwNused = iw(INwUse1+iroot0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwNarrays(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwnarrays FC_FUNC(iws_iwnarrays,IWS_IWNARRAYS)
CXXHFW    int fiws_iwnarrays(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwNarrays(int *obj)
CXXWRP    {
CXXWRP     int  ia = *(obj+1);
CXXWRP     int *iw = obj-ia;
CXXWRP     return fiws_iwnarrays(iw);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function iws_IwNarrays(iw)
C     ==================================

C--   Number of arrays in iw

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwNarrays = iw(NObjec1+iroot0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaLastObj(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_ialastobj FC_FUNC(iws_ialastobj,IWS_IALASTOBJ)
CXXHFW    int fiws_ialastobj(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaLastObj(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = fiws_ialastobj(iw);
CXXWRP     return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function iws_IaLastObj(iw)
C     ==================================

C--   Address of last object in iw

C--   Author: Michiel Botje h24@nikhef.nl   02-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iaR  = iroot0                           !Root address
      nobj = iw(NObjec1+iaR)                  !Number of arrays
      if(nobj.eq.0) then
        iws_IaLastObj = iaR                   !Return root when iw empty
      else
        iws_IaLastObj = iw(IWlstT1+iaR) + iaR !Return ia of last array
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwNheader(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwnheader FC_FUNC(iws_iwnheader,IWS_IWNHEADER)
CXXHFW    int fiws_iwnheader(int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwNheader(int *obj)
CXXWRP    {
CXXWRP     int  ia = *(obj+1);
CXXWRP     int *iw = obj-ia;
CXXWRP     return fiws_iwnheader(iw);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function iws_IwNheader(iw)
C     ==================================

C--   Header size

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwNheader = iw(INhead1+iroot0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwObjectType(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwobjecttype FC_FUNC(iws_iwobjecttype,IWS_IWOBJECTTYPE)
CXXHFW    int fiws_iwobjecttype(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwObjectType(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     return fiws_iwobjecttype(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ========================================
      integer function iws_IwObjectType(iw,ia)
C     ========================================

C--   Return object type

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      if(iw(ICword1+ia).eq.iWsVersion0)     then
        iws_IwObjectType = 1
      elseif(iw(ICword1+ia).eq.iCWTable0)   then
        iws_IwObjectType = 2
      else
        iws_IwObjectType = 0
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwAsize(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwasize FC_FUNC(iws_iwasize,IWS_IWASIZE)
CXXHFW    int fiws_iwasize(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwAsize(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     return fiws_iwasize(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      integer function iws_IwAsize(iw,ia)
C     ===================================

C--   Return array size

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwAsize = iw(INwUse1+ia)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwAnumber(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwanumber FC_FUNC(iws_iwanumber,IWS_IWANUMBER)
CXXHFW    int fiws_iwanumber(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwAnumber(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     return fiws_iwanumber(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     =====================================
      integer function iws_IwAnumber(iw,ia)
C     =====================================

C--   Return array number

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwAnumber = iw(IObjec1+ia)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwAfprint(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwafprint FC_FUNC(iws_iwafprint,IWS_IWAFPRINT)
CXXHFW    int fiws_iwafprint(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwAfprint(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     return fiws_iwafprint(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     =====================================
      integer function iws_IwAfprint(iw,ia)
C     =====================================

C--   Return fingerprint of root or array

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwAfprint = iw(IFprnt1+ia)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwAdim(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwadim FC_FUNC(iws_iwadim,IWS_IWADIM)
CXXHFW    int fiws_iwadim(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwAdim(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     return fiws_iwadim(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function iws_IwAdim(iw,ia)
C     ==================================

C--   Return K0 pointer coefficient

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      imi        = iw(IminTb1+ia)
      ima        = iw(ImaxTb1+ia)
      iws_IwAdim = ima-imi+1

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwAimin(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwaimin FC_FUNC(iws_iwaimin,IWS_IWAIMIN)
CXXHFW    int fiws_iwaimin(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwAimin(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     return fiws_iwaimin(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      integer function iws_IwAimin(iw,ia)
C     ===================================

C--   Return imin of array

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwAimin = iw(IminTb1+ia)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwAimax(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwaimax FC_FUNC(iws_iwaimax,IWS_IWAIMAX)
CXXHFW    int fiws_iwaimax(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwAimax(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     return fiws_iwaimax(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      integer function iws_IwAimax(iw,ia)
C     ===================================

C--   Return imax of array

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwAimax = iw(ImaxTb1+ia)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaAbegin(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iaabegin FC_FUNC(iws_iaabegin ,IWS_IAABEGIN)
CXXHFW    int fiws_iaabegin(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaAbegin(int *obj)
CXXWRP    {
CXXWRP     int  ja   = *(obj+1);
CXXWRP     int *iw   = obj-ja;
CXXWRP     int  ia   = iaCtoF(ja);
CXXWRP     int iak   = fiws_iaabegin(iw, &ia);
CXXWRP     return iaFtoC(iak);
CXXWRP    }
C-----------------------------------------------------------------------

C     ====================================
      integer function iws_IaAbegin(iw,ia)
C     ====================================

C--   First word of array body

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IaAbegin = iw(ItBody1+ia)+ia

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaAend(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iaaend FC_FUNC(iws_iaaend ,IWS_IAAEND)
CXXHFW    int fiws_iaaend(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaAend(int *obj)
CXXWRP    {
CXXWRP     int  ja   = *(obj+1);
CXXWRP     int *iw   = obj-ja;
CXXWRP     int  ia   = iaCtoF(ja);
CXXWRP     int iak   = fiws_iaaend(iw, &ia);
CXXWRP     return iaFtoC(iak);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function iws_IaAend(iw,ia)
C     ==================================

C--   Last word of array body

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IaAend = iw(ItBend1+ia)+ia

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IwKnul(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iwknul FC_FUNC(iws_iwknul,IWS_IWKNUL)
CXXHFW    int fiws_iwknul(int*,int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IwKnul(int *obj)
CXXWRP    {
CXXWRP     int  ja = *(obj+1);
CXXWRP     int *iw = obj-ja;
CXXWRP     int  ia = iaCtoF(ja);
CXXWRP     return fiws_iwknul(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==================================
      integer function iws_IwKnul(iw,ia)
C     ==================================

C--   Return K0 pointer coefficient

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      iws_IwKnul = iw(KnulTb1+ia)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_ArrayI(int *obj, int i);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_arrayi FC_FUNC(iws_arrayi ,IWS_ARRAYI)
CXXHFW    int fiws_arrayi(int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_ArrayI(int *obj, int i)
CXXWRP    {
CXXWRP     int  ja   = *(obj+1);
CXXWRP     int *iw   = obj-ja;
CXXWRP     int  ia   = iaCtoF(ja);
CXXWRP     int iak   = fiws_arrayi(iw, &ia, &i);
CXXWRP     return iaFtoC(iak);
CXXWRP    }
C-----------------------------------------------------------------------

C     ====================================
      integer function iws_ArrayI(iw,ia,i)
C     ====================================

C--   Address of array(i)

C--   Author: Michiel Botje h24@nikhef.nl   01-02-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      imi        = iw(IminTb1+ia)
      ima        = iw(ImaxTb1+ia)
      if(i.ge.imi .and. i.le.ima) then
        iws_ArrayI = iw(KnulTb1+ia)+ia+i
      else
        stop 'WSTORE:IWS_ARRAYI: index I not in range'
      endif

      return
      end

C=======================================================================
C==  Verbose dumps  ====================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    void sws_IwTree(int *iw);
C-----------------------------------------------------------------------
CXXHFW  #define fswsiwtree FC_FUNC(swsiwtree ,SWSIWTREE)
CXXHFW    int fswsiwtree(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_IwTree(int *iw)
CXXWRP    {
CXXWRP     int iroot    = 0;
CXXWRP     fswsiwtree(iw,&iroot);
CXXWRP    }
C-----------------------------------------------------------------------

C     =========================
      subroutine sws_IwTree(iw)
C     =========================

C--   Print istore tree

C--   Author: Michiel Botje h24@nikhef.nl   31-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

      call swsIWtree(iw,1)

      return
      end

C     ===============================
      subroutine  swsIWtree(iw,iroot)
C     ===============================

C--   The iroot argumant takes care of the address printouts
C--   Fortran iroot = 1 --> print address as is
C--   C++     iroot = 0 --> substract 1 from all addresses

C--   Author: Michiel Botje h24@nikhef.nl   31-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

C--   Check input
      if(iw(iroot0).ne.iWsVersion0) stop
     +                        'WSTORE:SWS_IWTREE: IW is not an istore'
      if(iroot.ne.0 .and. iroot.ne.1) stop
     +                        'WSTORE:SWS_IWTREE: iroot must be 0 or 1'

      ia = 1                                 !address of istore
      call swsiwprnt(iw,ia,iroot)            !print istore
      jd = iw(NFTabl1+ia)                    !distance to first array
      do while(jd.ne.0)                      !loop over arrays
        ia = ia+jd                           !address of array
        call swsaprint(iw,ia,iroot)          !print array
        jd = iw(NFTabl1+ia)                  !distance to next array
      enddo                                  !end loop over arrays

      return
      end

C     =================================
      subroutine swsiwprnt(iw,ia,iroot)
C     =================================

C--   The iroot argumant takes care of the address printouts
C--   Fortran iroot = 1 --> print address as is
C--   C++     iroot = 0 --> substract 1 from all addresses

C--   Author: Michiel Botje h24@nikhef.nl   31-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)
      character*15 number, hcode
      character*80 text

C--   Check input address
      if(iw(ia).ne.iWsVersion0) stop
     +                  'WSTORE:SWS_IWTREE: IA is not an istore address'

      iadr  = ia-iroot0+iroot
      nwds  = iw(INwUse1+ia)
      ifpt  = iw(IFprnt1+ia)
      nobj  = iw(NObjec1+ia)
      call smb_itoch(nobj,number,leng)
      call smb_hcode(ifpt,hcode)
      call smb_cfill(' ',text)
      text = 'istore with '//number(1:leng)//' arrays'
      ltxt = imb_lastc(text)

      write(6,'(/1X,''ADDRESS'',4X,''SIZE'',8X,''FINGERPRINT'',
     +           4X,''OBJECT'')')
      write(6,'(2I8,4X,A15,4X,A)') iadr, nwds, hcode, text(1:ltxt)

      return
      end

C     =================================
      subroutine swsaprint(iw,ia,iroot)
C     =================================

C--   The iroot argumant takes care of the address printouts
C--   Fortran iroot = 1 --> print address as is
C--   C++     iroot = 0 --> substract 1 from all addresses

C--   Author: Michiel Botje h24@nikhef.nl   31-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)
      character*15 number, hcode
      character*80 text

C--   Check input address
      if(iw(ia).ne.iCWTable0) stop
     +                   'WSTORE:SWS_IWTREE: IA is not an array address'

      iadr  = ia-iroot0+iroot
      nwds  = iw(INwUse1+ia)
      ifpt  = iw(IFprnt1+ia)
      iab1  = iw(ItBody1+ia)
      iab2  = iw(ItBend1+ia)
      nelm  = iab2 - iab1 + 1
      call smb_itoch(nelm,number,leng)
      call smb_hcode(ifpt,hcode)
      call smb_cfill(' ',text)
      text = ' array with '//number(1:leng)//' elements'
      ltxt = imb_lastc(text)

      write(6,'(2I8,4X,A15,4X,A)') iadr, nwds, hcode, text(1:ltxt)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void sws_IwHead(int *obj);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_iwhead FC_FUNC(sws_iwhead ,SWS_IWHEAD)
CXXHFW    int fsws_iwhead(int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_IwHead(int *obj)
CXXWRP    {
CXXWRP     int ja    = *(obj+1);
CXXWRP     int *iw   = obj-ja;
CXXWRP     int ia    = iaCtoF(ja);
CXXWRP     fsws_iwhead(iw,&ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ============================
      subroutine sws_IwHead(iw,ia)
C     ============================

C--   Print header of istore, table-set or table

C--   Author: Michiel Botje h24@nikhef.nl   31-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension iw(*)

C--   Check input
      if(iw(iroot0).ne.iWsVersion0) stop
     +                        'WSTORE:SWS_IWHEAD: IW is not an istore'
      if(ia.le.0 .or. ia.gt.iw(INwUse1+iroot0)) stop
     +                        'WSTORE:SWS_IWHEAD: IA out of range'

      if(iw(ia).eq.iWsVersion0) then
         write(6,'(/'' Istore Header'')')
         write(6,'( '' 0 Cword    '',I15  )') iw(ICword1+ia)
         write(6,'( '' 1 IW       '',I15  )') iw(IwAddr1+ia)
         write(6,'( '' 2 TFskip   '',I15  )') iw(NFTabl1+ia)
         write(6,'( '' 3 TBskip   '',I15  )') iw(NBTabl1+ia)
         write(6,'( '' 4 Fprint   '',I15  )') iw(IFprnt1+ia)
         write(6,'( '' 5 NWused   '',I15  )') iw(INwUse1+ia)
         write(6,'( '' 6 Nobj     '',I15  )') iw(NObjec1+ia)
         write(6,'( '' 7 IW Ltab  '',I15  )') iw(IWlstT1+ia)
         write(6,'( '' 8 NWtotal  '',I15  )') iw(iNwMax1+ia)
         write(6,'( '' 9 Nheader  '',I15  )') iw(INhead1+ia)
       elseif(iw(ia).eq.iCWTable0) then
         write(6,'(/'' Array Header'')')
         write(6,'( '' 0 Cword    '',I15  )') iw(ICword1+ia)
         write(6,'( '' 1 IW       '',I15  )') iw(IwAddr1+ia)
         write(6,'( '' 2 TFskip   '',I15  )') iw(NFTabl1+ia)
         write(6,'( '' 3 TBskip   '',I15  )') iw(NBTabl1+ia)
         write(6,'( '' 4 Fprint   '',I15  )') iw(IFprnt1+ia)
         write(6,'( '' 5 NWused   '',I15  )') iw(INwUse1+ia)
         write(6,'( '' 6 Iobj     '',I15  )') iw(IObjec1+ia)
         write(6,'( '' 7 K0       '',I15  )') iw(KnulTb1+ia)
         write(6,'( '' 8 Imin     '',I15  )') iw(IminTb1+ia)
         write(6,'( '' 9 Imax     '',I15  )') iw(ImaxTb1+ia)
         write(6,'( ''10 IT Bbody '',I15  )') iw(ItBody1+ia)
         write(6,'( ''11 IT Ebody '',I15  )') iw(ItBend1+ia)
      else
         stop 'WSTORE:SWS_IWHEAD: IA is not a root or array address'
      endif

      return
      end



