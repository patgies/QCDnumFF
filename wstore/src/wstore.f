
C--  This is the file wstore.f with workspace routines
C--  This file also contains C++ wrappers for all the Fortran routines
C--
C--   iws_Version()                       Return wstore version number
C--
c--   iwsWsInitCPP(w,nw,nt,txt,ls)        Wrapper for WsInit in C++
C--   iws_WsInit(w,nw,nt,txt)             Create workspace
C--   swsWsEbuf(w,txt,opt)                Handle message text buffer
C--   swsWsEmsg(w,n,srname)               Out-of-space error message
C--   sws_SetWsN(w,nw)                    Update workspace size info
C--   sws_Stampit(w)                      Set new time stamp
C--   iws_NewSet(w)                       Create new table-set
C--   iwsEtrailer(w)                      Yes/no w has empty trailer set
C--   iws_WTable(w,imi,ima,n)             Create new table
C--   iws_WClone(w1,ia1,w2)               Clone a table-set or table
C--   iwsTClone(w1,ia1,w2)                Clone a table
C--   sws_TbCopy(w1,ia1,w2,ia2,itag)      Copy table contents
C--   sws_WsWipe(w,ia)                    Wipe workspace
C--   swsTsDumpCPP(fnam,ls,key,w,ia,ierr) Wrapper for sws_TsDump in C++
C--   sws_TsDump(fname,key,w,ia,ierr)     Dump table-set to disk
C--   iwsTsReadCPP(fnam,ls,key,w,ierr)    Wrapper for sws_TsRead in C++
C--   iws_TsRead(fname,key,w,ierr)        Append table set from disk
C--   sws_WsMark(mws,mset,mtab)           Get markers
C--   iws_HdSize()                        Header size
C--   iws_TbSize(imi,ima,ndim)            Table size

C--   swsGetMeta(w,ia,nd,kk,imi,ima)      Extract metadata
C--   iws_Tpoint(w,iw,index,n)            Pointer function

C--   iws_IsaWorkspace(w)                 True if w is a workspace
C--   iws_SizeOfW(w)                      Total size of w
C--   iws_WordsUsed(w)                    Words used
C--   iws_Nheader(w)                      Header size
C--   iws_Ntags(w)                        Tag field size
C--   iws_HeadSkip(w)                     Header + tag field size
C--   iws_IaRoot()                        Root address
C--   iws_IaDrain(w)                      Drain-word address
C--   iws_IaNull(w)                       Null-word address
C--   iws_ObjectType(w,ia)                Object type
C--   iws_ObjectSize(w,ia)                Object size
C--   iws_Nobjects(w,ia)                  Number of objects
C--   iws_ObjectNumber(w,ia)              Objects number
C--   iws_FingerPrint(w,ia)               Fingerprint
C--   iws_IaFirstTag(w,ia)                Tagfield address
C--   iws_TableDim(w,ia)                  Number of table dimensions
C--   iws_IaKARRAY(w,ia)                  KARRAY first word
C--   iws_IaIMIN(w,ia)                    IMIN first word
C--   iws_IaIMAX(w,ia)                    IMAX first word
C--   iws_BeginTbody(w,ia)                Begin table body
C--   iws_EndTbody(w,ia)                  End table body

C--   iws_TFskip(w,ia)                    Forward skip table
C--   iws_TBskip(w,ia)                    Reverse skip table
C--   iws_SFskip(w,ia)                    Forward skip tset
C--   iws_SBskip(w,ia)                    Reverse skip tset

C--   sws_WsTree(w)                       Dump workspace tree
C--   sws_WsHead(w,ia)                    Dump object header


C--   Block of words: frst, last, next, size
C--
C--   size  =  last - frst + 1  =  next - frst   -->   size >= 1
C--   next  =  frst + size      =  last + 1
C--   last  =  frst + size - 1  =  next - 1
C--   frst  =  last - size + 1  =  next - size

C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  WSTORE workspace routines from wstore.f                 */
CXXHDR    /************************************************************/
CXXHDR
CXXHDR    inline int iaFtoC(int ia) { return ia-1; };
CXXHDR    inline int iaCtoF(int ia) { return ia+1; };
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  WSTORE workspace routines from wstore.f                   */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  WSTORE workspace routines from wstore.f                   */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    int iws_Version();
C-----------------------------------------------------------------------
CXXHFW  #define fiws_version FC_FUNC(iws_version,IWS_VERSION)
CXXHFW    int fiws_version();
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_Version()
CXXWRP    {
CXXWRP      return fiws_version();
CXXWRP    }
C-----------------------------------------------------------------------

C     ==============================
      integer function iws_Version()
C     ==============================

      implicit double precision(a-h,o-z)

      include 'wstore0.inc'

      iws_Version = iWSversion0

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_WsInit(double *w, int nw, int nt, string txt);
C-----------------------------------------------------------------------
CXXHFW  #define fiwswsinitcpp FC_FUNC(iwswsinitcpp,IWS_WSINITCPP)
CXXHFW    int fiwswsinitcpp(double*, int*, int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_WsInit(double *w, int nw, int nt, string txt)
CXXWRP    {
CXXWRP      int ls = txt.size();
CXXWRP      char *ctxt = new char[ls+1];
CXXWRP      strcpy(ctxt,txt.c_str());
CXXWRP      int ia = fiwswsinitcpp(w, &nw, &nt, ctxt, &ls);
CXXWRP      delete[] ctxt;
CXXWRP      return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================================
      integer function iwsWsInitCPP(w,nw,nt,txt,ls)
C     =============================================

      implicit double precision (a-h,o-z)

      dimension w(*)
      character*(100) txt

      if(ls.gt.100) stop
     +             'WSTORE:IWS_WSINIT: input text > 100 characters'

      iwsWsInitCpp = iws_WsInit(w,nw,nt,txt(1:ls))

      return
      end

C     ========================================
      integer function iws_WsInit(w,nw,nt,txt)
C     ========================================

C--   Create new workspace with first (empty) table-set
C--
C--   w           (in): double precision array w(nw)
C--   nw          (in): dimension of w declared in calling routine
C--   nt          (in): requested tag field size
C--   txt         (in): comment line when w runs out of space
C--   iws_WsInit (out): IA address of the first table-set

C--   Author: Michiel Botje h24@nikhef.nl   17-11-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      character*(*) txt
      character*10  date, time, zone
      dimension     ival(8)
      character*20  ntxt

      save icnt
      data icnt /0/

C--   Check input
      if(nw.le.0) stop
     + 'WSTORE:IWS_WSINIT: cannot have workspace size NW <= 0'
      if(nt.lt.0) stop
     + 'WSTORE:IWS_WSINIT: cannot have tag field size NT < 0'
C--   Check size (for Wspace and first table set)
      nhskip = nwHeader0 + nt
      nsize  = nhskip + nhskip
      nneed  = nsize + 1
      if(nw.lt.nneed) then
        call smb_itoch(nneed,ntxt,leng)
        write(6,*)
     + 'WSTORE:IWS_WSINIT: workspace size must be at least ',
     +  ntxt(1:leng),' words'
        if(imb_lastc(txt).ne.0) write(6,*) txt
        stop
      endif
C--   Workspace fingerprint
      call date_and_time(date,time,zone,ival)      !millisec resolution
      icnt  = icnt+1                               !count calls
      iseed = 0                                    !call-dependent seed
      do i = 1,4
        call smb_cbyte(mod(icnt+i,256),1,iseed,i)  !set 4 bytes of iseed
      enddo
      ihash = imb_ihash(iseed,ival,8)              !hash year-date-time
C--   Initialise store
      call smb_Vfill(w,nw,0.D0)
C--   Workspace header
      w(ICword0+iroot0) = dble(iCWorkspace0)    !Control word
      w(IwAddr0+iroot0) = 0.D0                  !IW address
      w(NFTabl0+iroot0) = 0.D0                  !Fskip table
      w(NBTabl0+iroot0) = 0.D0                  !Bskip table
      w(NFTset0+iroot0) = dble(nhskip)          !Fskip tset
      w(NBTset0+iroot0) = 0.D0                  !Bskip tset
      w(Ifprnt0+iroot0) = dble(ihash)           !Fingerprint
      w(NObjec0+iroot0) = 1.D0                  !Number of objects
      w(IObjec0+iroot0) = 1.D0                  !Object number
      w(IWlstS0+iroot0) = dble(nhskip)          !IW last table-set in w
      w(IWlstT0+iroot0) = dble(nhskip+nhskip)   !IW last table in w
      w(INwMax0+iroot0) = dble(nw)              !Wspace total size
      w(INwUse0+iroot0) = dble(nsize)           !Words used
      w(IDrain0+iroot0) = 0.D0                  !Drain word
      w(iZnull0+iroot0) = 1.D20                 !Null word
C--   Initial table-set fingerprint
      ival(1) = nwHeader0
      ival(2) = nt
      ihash   = 0
      ihash   = imb_ihash(ihash,ival,2)
C--   First table-set address
      iaS     = nhskip + iroot0
C--   Fill Header first table set
      w(ICword0+iaS) = dble(iCWTableSet0)       !control word
      w(IwAddr0+iaS) = dble(nhskip)             !IW
      w(NFTabl0+iaS) = 0.D0                     !Fskip table
      w(NBTabl0+iaS) = 0.D0                     !Bskip table
      w(NFTset0+iaS) = 0.D0                     !Fskip tset
      w(NBTset0+iaS) = 0.D0                     !Bskip tset
      w(Ifprnt0+iaS) = dble(ihash)              !fingerprint
      w(NObjec0+iaS) = 0.D0                     !number of tables
      w(IObjec0+iaS) = 1.D0                     !Object number
      w(INwUse0+iaS) = dble(nhskip)             !Oject  size
      w(INhead0+iaS) = dble(nwHeader0)          !Header size
      w(INtags0+iaS) = dble(nt)                 !Tag size
      w(IHskip0+iaS) = dble(nhskip)             !Hskip
      w(ISlstT0+iaS) = dble(nhskip)             !IS first table in ts
C--   Store message
      call swsWsEbuf(w,txt,'in')
C--   Return address first table-set
      iws_WsInit = iaS

      return
      end

C     ===============================
      subroutine swsWsEbuf(w,txt,opt)
C     ===============================

C--   Store (opt = 'in') or retrieve (opt = 'out') message text

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)
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
      ifpw = int(w(IFprnt0+iroot0))
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
            write(6,*) 'WSTORE:IWS_WSINIT: message buffer size exceeded'
            write(6,*) '-- Increase MEBUF0 in wstore/inc/wstore0.inc'
            write(6,*) '-- and recompile WSTORE'
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
        stop 'WSTORE:swsWsEbuf: unknown option'
      endif

      return

      end

C     ================================
      subroutine swsWsEmsg(w,n,srname)
C     ================================

C--   Out-of-space message

      implicit double precision(a-h,o-z)

      dimension w(*)
      character*(*) srname
      character*20  ntxt
      character*80  txt

      i1 = imb_frstc(srname)
      i2 = imb_lastc(srname)
      call smb_itoch(n,ntxt,leng)
      call swsWsEbuf(w,txt,'out')
      write(6,*) srname(i1:i2),': workspace size must be at least ',
     +  ntxt(1:leng),' words'
      if(imb_lastc(txt).ne.0) write(6,*) txt
      stop

      end

C-----------------------------------------------------------------------
CXXHDR    void sws_SetWsN(double *w, int nw);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_setwsn FC_FUNC(sws_setwsn,SWS_SETWSN)
CXXHFW    void fsws_setwsn(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_SetWsN(double *w, int nw)
CXXWRP    {
CXXWRP     fsws_setwsn(w, &nw);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===========================
      subroutine sws_SetWsN(w,nw)
C     ===========================

C--   Update size information in the workspace header
C--
C--   w               (in): workspace
C--   nw              (in): actual workspace size

C--   Author: Michiel Botje h24@nikhef.nl   25-02-20

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

C--   Check if workspace
      if(int(w(iCword0+iroot0)).ne.iCWorkspace0) then
        stop 'WSTORE:SWS_SETWSN: W is not a workspace'
      endif
C--   Check input
      nwuse = int(w(INwUse0+iroot0))+1
      if(nw.lt.nwuse) stop
     +    'WSTORE:SWS_SETWSN: cannot set NW < number of words used'
C--   Update header
      w(INwMax0+iroot0) = dble(nw)              !Wspace total size

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void sws_Stampit(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_stampit FC_FUNC(sws_stampit,SWS_STAMPIT)
CXXHFW    void fsws_stampit(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_Stampit(double *w)
CXXWRP    {
CXXWRP     fsws_stampit(w);
CXXWRP    }
C-----------------------------------------------------------------------

C     =========================
      subroutine sws_Stampit(w)
C     =========================

C--   Set new timetamp in workspace w

C--   Author: Michiel Botje h24@nikhef.nl   05-11-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      character*10  date, time, zone
      dimension     ival(8),ifpold(1)

      save icnt
      data icnt /0/

C--   Check if workspace
      if(int(w(iCword0+iroot0)).ne.iCWorkspace0) then
        stop 'WSTORE:SWS_STAMPIT: W is not a workspace'
      endif
C--   Old timestamp
      ifpold = int(w(IFprnt0+iroot0))
C--   New timestamp
      call date_and_time(date,time,zone,ival)      !millisec resolution
      icnt  = icnt+1                               !count calls
      iseed = 0                                    !call-dependent seed
      do i = 1,4
        call smb_cbyte(mod(icnt+i,256),1,iseed,i)  !set 4 bytes of iseed
      enddo
      ihash = imb_ihash(iseed,ival,8)              !hash year-date-time
      ihash = imb_ihash(ihash,ifpold,1)            !add old timetamp

C--   Update header
      w(IFprnt0+iroot0) = dble(ihash)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_NewSet(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_newset FC_FUNC(iws_newset,IWS_NEWSET)
CXXHFW    int fiws_newset(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_NewSet(double *w)
CXXWRP    {
CXXWRP     int ia = fiws_newset(w);
CXXWRP     return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==============================
      integer function iws_NewSet(w)
C     ==============================

C--   Book new table-set (if current one not empty)
C--
C--   w               (in): workspace
C--   iws_NewSet     (out): IA address of new table-set

C--   Author: Michiel Botje h24@nikhef.nl   17-11-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*), ival(2)

C--   Check if workspace
      if(int(w(iCword0+iroot0)).ne.iCWorkspace0) then
        stop 'WSTORE:IWS_NEWSET: W is not a workspace'
      endif
C--   Get info
      iaR    = iroot0                             !root address
      iaL    = int(w(IWlstS0+iaR)) + iaR          !IA current table-set
      iaS    = int(w(INwUse0+iaR)) + iaR          !IA new     table-set
      iaT    = int(w(IWlstT0+iaR)) + iaR          !IA last table
      NBskip = iaL-iaS                            !Bskip tset
      NTskip = iaT-iaS                            !Bskip table
      NObjec = int(w(NObjec0+iaR))                !# table-sets in w
C--   Do nothing if current set is empty
      if( iwsEtrailer(w) .eq.1 ) then
        iws_NewSet = iaL
        return
      endif
C--   Check workspace size
      NHskip = int(w(IHskip0+iaL))                !header+tag size
      NWused = int(w(INwUse0+iaR))                !words used
      NWneed = NWused + NHskip + 1                !words needed
      Nspace = int(w(INwMax0+iaR))                !total size
C--   Run out-of-space
      if(NWneed.gt.Nspace) call swsWsEmsg(w,NWneed,'WSTORE:IWS_NEWSET')
C--   Initial fingerprint
      Nhead   = int(w(INhead0+iaL))
      Ntags   = int(w(INtags0+iaL))
      ival(1) = Nhead
      ival(2) = Ntags
      ihash   = 0
      ihash   = imb_ihash(ihash,ival,2)
      NObjec  = NObjec + 1                        !update #table sets
C--   Fill Header
      w(ICword0+iaS) = dble(iCWTableSet0)         !control word
      w(IwAddr0+iaS) = dble(NWused)               !IW
      w(NFTabl0+iaS) = 0.D0                       !Fskip table
      w(NBTabl0+iaS) = dble(NTskip)               !Bskip table
      w(NFTset0+iaS) = 0.D0                       !Fskip tset
      w(NBTset0+iaS) = dble(NBskip)               !Bskip tset
      w(Ifprnt0+iaS) = dble(ihash)                !fingerprint
      w(NObjec0+iaS) = 0.D0                       !number of tables
      w(IObjec0+iaS) = dble(NObjec)               !Object number
      w(INwUse0+iaS) = dble(NHskip)               !Oject  size
      w(INhead0+iaS) = dble(Nhead)                !Header size
      w(INtags0+iaS) = dble(Ntags)                !Tag size
      w(IHskip0+iaS) = dble(NHskip)               !Hskip
      w(ISlstT0+iaS) = dble(NHskip)               !IS first table in ts
C--   Update Wspace Header
      w(NObjec0+iaR) = dble(NObjec)               !Number of objects
      w(IWlstS0+iaR) = dble(NWused)               !IW this table set
      w(IWlstT0+IaR) = dble(NWused+NHskip)        !IW first table in ts
      w(INwUse0+iaR) = dble(NWused+NHskip)        !Words used
C--   Return address
      iws_NewSet = iaS
C--   Done for first set
      if(NBskip.eq.0) return
C--   Update Header previous set
      w(NFTset0+iaL) = dble(iaS - iaL)            !SFskip tset
C--   Update Header previous tables
      nTab  = int(w(iaL+NObjec0))
      iaTab = int(w(iaL+NFTabl0)) + iaL
      do i = 1,nTab
        w(NFTset0+iaTab) = dble(iaS - iaTab)      !SFskip table in tset
        iaTab = int(w(iaTab+NFTabl0)) + iaTab
      enddo

      return
      end

C     ===============================
      integer function iwsEtrailer(w)
C     ===============================

C--   Return 1 (0) if w has (not) a trailing empty table-set
C--
C--   Author: Michiel Botje h24@nikhef.nl   05-03-20

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iaR    = iroot0                             !root address
      iaS    = int(w(IWlstS0+iaR)) + iaR          !IA trailing table-set
      NTabs  = int(w(NObjec0+iaS))                !# tables trailing set

      if(NTabs.ne.0) then
        iwsEtrailer = 0
      else
        iwsEtrailer = 1
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_WTable(double* w, int* imi, int* ima, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_wtable FC_FUNC(iws_wtable,IWS_WTABLE)
CXXHFW    int fiws_wtable(double*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_WTable(double* w, int* imi, int* ima, int n)
CXXWRP    {
CXXWRP     int ia = fiws_wtable(w, imi, ima, &n);
CXXWRP     return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ========================================
      integer function iws_WTable(w,imi,ima,n)
C     ========================================

C--   Book new table
C--
C--   w           (in): workspace
C--   imi,ima     (in): index limits
C--   n           (in): number of dimensions [1,mdim0]
C--   iws_WTable (out): IA address of the new table

C--   Author: Michiel Botje h24@nikhef.nl   17-11-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*), imi(*), ima(*), ival(2)

C--   Check if workspace
      if(int(w(iCword0+iroot0)).ne.iCWorkspace0) then
        stop 'WSTORE:IWS_WTABLE: W is not a workspace'
      endif
C--   Check input and get table size
      if(n.lt.1 .or. n.gt.mdim0) then
        stop 'WSTORE:IWS_WTABLE: Ndim not in range [1,25]'
      endif
      nw = 1
      do i = 1,n
        if(imi(i).gt.ima(i)) stop 'WSTORE:IWS_WTABLE: imin > imax'
        nw = nw*(ima(i)-imi(i)+1)
      enddo
C--   Addresses
      iaR    = iroot0                             !root address
      iaS    = int(w(IWlstS0+iaR)) + iaR          !IA current table-set
      iaL    = int(w(IWlstT0+iaR)) + iaR          !IA current table
      iaT    = int(w(INwUse0+iaR)) + iaR          !IA new table
      NBskip = iaL-iaT                            !Bskip table
      NSskip = iaS-iaT                            !Bskip table set
C--   Check workspace size
      NSsize = int(w(INwUse0+iaS))                !table-set size
      NHskip = int(w(IHskip0+iaS))                !header+tag size
      NWused = int(w(INwUse0+iaR))                !words used
      NTsize = nw + 3*n + 2 + NHskip              !table size
      NWneed = NWused + NTsize + 1                !words needed
      Nspace = int(w(INwMax0+iaR))                !total size
C--   Run out-of-space
      if(NWneed.gt.Nspace) call swsWsEmsg(w,NWneed,'WSTORE:IWS_WTABLE')
C--   Initialise
      do i = iaT,NWneed
        w(i) = 0.D0
      enddo
C--   Local addresses of first word of metadata items
      itnd1          = NHskip
      itka1          = NHskip+1
      itmi1          = NHskip+n+2
      itma1          = NHskip+n+n+2
      ittb1          = NHskip+n+n+n+2
C--   Fill metadata
      w(iaT+itnd1)   = dble(n)
      call smb_dkmat(imi,ima,w(iaT+itka1),n,ittb1,ittb2)
      if(ittb2.ne.NTsize-1) stop
     +                     'WSTORE:IWS_WTABLE: problem with table size'
      call smb_Vitod(imi,w(iaT+itmi1),n)
      call smb_Vitod(ima,w(iaT+itma1),n)
C--   Table fingerprint
      ival(1)        = n
      jh             = imb_ihash(0,ival,1)
      jh             = imb_jhash(jh,w(iaT+itka1),n+1)
      jh             = imb_ihash(jh,imi,n)
      jh             = imb_ihash(jh,ima,n)
      ihash          = imb_jhash(0,w(iaT+itnd1),n+n+n+2)
      if(ihash.ne.jh) stop
     +               'WSTORE:IWS_WTABLE: problem with fingerprint'
      NObjec = int(w(NObjec0+iaS)) + 1            !# tables
C--   Fill Header
      w(ICword0+iaT) = dble(iCWTable0)            !control word
      w(IwAddr0+iaT) = dble(NWused)               !IW
      w(NFTabl0+iaT) = 0.D0                       !Fskip table
      w(NBTabl0+iaT) = dble(NBskip)               !Bskip table
      w(NFTset0+iaT) = 0.D0                       !Fskip tset
      w(NBTset0+iaT) = dble(NSskip)               !Bskip tset
      w(Ifprnt0+iaT) = dble(ihash)                !fingerprint
      w(NObjec0+iaT) = 0.D0                       !number of objects
      w(IObjec0+iaT) = dble(NObjec)               !Object number
      w(INwUse0+iaT) = dble(NTsize)               !Table size
      w(ItMeta0+iaT) = dble(itnd1)                !IT Metadata
      w(ItImin0+iaT) = dble(itmi1)                !IT IMIN
      w(ItImax0+iaT) = dble(itma1)                !IT IMAX
      w(ItBody0+iaT) = dble(ittb1)                !IT begin body
      w(ItBend0+iaT) = dble(ittb2)                !IT end body
C--   Update Wspace header
      w(IWlstT0+iaR) = dble(NWused)               !IW this table
      w(INwUse0+iaR) = dble(NWused+NTsize)        !Words used
      w(NFTabl0+iaR) = dble(NHskip+NHskip)        !Fskip table
C--   Update table-set fingerprint
      jhash   = int(w(Ifprnt0+iaS))               !old fingerprint
      ival(1) = ihash                             !table fingerprint
      jhash   = imb_ihash(jhash,ival,1)           !new fingerprint
C--   Update table set header
      w(NFTabl0+iaS) = dble(NHskip)               !Fskip table
      w(Ifprnt0+iaS) = dble(jhash)                !fingerprint
      w(NObjec0+iaS) = dble(NObjec)               !Number of objects
      w(INwUse0+iaS) = dble(NSsize+NTsize)        !Table-set size
      w(ISlstT0+iaS) = dble(NSsize)               !IS this table
C--   Return address
      iws_WTable = iaT
C--   Done for first table
      if(NBskip.eq.0) return
C--   Update Header previous table
      w(NFTabl0+iaL) = iaT - iaL                  !Fskip table

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_WClone(double *w1, int ia1, double *w2);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_wclone FC_FUNC(iws_wclone,IWS_WCLONE)
CXXHFW    int fiws_wclone(double*, int*, double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_WClone(double *w1, int ia1, double *w2)
CXXWRP    {
CXXWRP     int ja1    = iaCtoF(ia1);
CXXWRP     int ia     = fiws_wclone(w1, &ja1, w2);
CXXWRP     return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ========================================
      integer function iws_WClone(w1, ia1, w2)
C     ========================================

C--   Clone a table-set or table

C--   Author: Michiel Botje h24@nikhef.nl   09-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w1(*), w2(*)

C--   Check input
      if(int(w1(iCword0+iroot0)).ne.iCWorkspace0) stop
     +                        'WSTORE:IWS_WCLONE: W1 is not a workspace'
      if(ia1.le.0 .or. ia1.gt.int(w1(INwUse0+iroot0))) stop
     +                        'WSTORE:IWS_WCLONE: IA1 out of range'
      if(int(w2(iCword0+iroot0)).ne.iCWorkspace0) stop
     +                        'WSTORE:IWS_WCLONE: W2 is not a workspace'
      if(int(w1(iCword0+ia1)).ne.iCWTableSet0 .and.
     +   int(w1(iCword0+ia1)).ne.iCWTable0) stop
     +  'WSTORE:IWS_WCLONE: object to clone is not a table-set or table'

C--   Check workspace size
      iaR    = iroot0                             !workspace base
      iaS    = int(w1(IWlstS0+iaR)) + iaR         !IA current table-set
      NHskip = int(w1(IHskip0+iaS))               !header+tag size
      NCsize = int(w1(INwUse0+ia1))               !object size
      NWused = int(w2(INwUse0+iaR))               !words used
      Nspace = int(w2(INwMax0+iaR))               !total size
      if(iwsEtrailer(w2).eq.1) then
        NWneed = NWused + NCsize - NHskip + 1     !words needed
      else
        NWneed = NWused + NCsize + 1              !words needed
      endif
C--   Run out-of-space
      if(NWneed.gt.Nspace) call swsWsEmsg(w2,NWneed,'WSTORE:IWS_WCLONE')

      if(int(w1(ICword0+ia1)).eq.iCWTableSet0) then
C--     Clone table set
C--     New table set in w2
        iaS   = iws_NewSet(w2)
C--     Clone w1 tables to w2
        iaT1  = ia1 + int(w1(ia1+NFTabl0))
        iaT2  = iwsTClone(w1,iaT1,w2)
        iskip = int(w1(iaT1+NFTabl0))
        do while(iskip.ne.0)
          iaT1  = iaT1 + iskip
          iaT2  = iwsTClone(w1,iaT1,w2)
          iskip = int(w1(iaT1+NFTabl0))
        enddo
        iws_WClone = iaS
      else
C--     Clone table
        iws_WClone = iwsTClone(w1,ia1,w2)
      endif

      return
      end

C     =======================================
      integer function iwsTClone(w1, ia1, w2)
C     =======================================

C--   Clone table

C--   Author: Michiel Botje h24@nikhef.nl   02-03-20

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w1(*), w2(*), ival(2)

C--   Workspace size
      iaR    = iroot0                             !workspace base
      NCsize = int(w1(INwUse0+ia1))               !object size
      NWused = int(w2(INwUse0+iaR))               !words used

C--   Clone table
      iaS    = int(w2(IWlstS0+iaR)) + iaR         !IA this table-set
      iaL    = int(w2(IWlstT0+iaR)) + iaR         !IA current table
      iaT    = int(w2(INwUse0+iaR)) + iaR         !IA new     table
      NBskip = iaL-iaT                            !Bskip table
      NSskip = iaS-iaT                            !Bskip table set
C--   Clone
      call smb_vcopy(w1(ia1),w2(iaT),NCsize)
      NSsize = int(w2(INwUse0+iaS))               !table-set size
      NHskip = int(w2(IHskip0+iaS))               !header+tag size
      NObjec = int(w2(NObjec0+iaS)) + 1           !update #tables
C--   Update table header
      w2(IwAddr0+iaT) = dble(NWused)              !IW
      w2(NFTabl0+iaT) = 0.D0                      !Fskip table
      w2(NBTabl0+iaT) = dble(NBskip)              !Bskip table
      w2(NFTset0+iaT) = 0.D0                      !Fskip tset
      w2(NBTset0+iaT) = dble(NSskip)              !Bskip tset
      w2(IObjec0+iaT) = dble(NObjec)              !Object number
C--   Update Wspace Header
      w2(IWlstT0+iaR) = dble(NWused)              !IW this table
      w2(INwUse0+iaR) = dble(NWused+NCsize)       !Words used
C--   Update table-set fingerprint
      jhash   = int(w2(Ifprnt0+iaS))              !old fingerprint
      ival(1) = int(w2(IFprnt0+iaT))              !table fingerprint
      jhash   = imb_ihash(jhash,ival,1)           !new fingerprint
C--   Update table-set header
      w2(NFTabl0+iaS) = dble(NHskip)              !Fskip table
      w2(Ifprnt0+iaS) = dble(jhash)               !fingerprint
      w2(NObjec0+iaS) = dble(NObjec)              !Number of objects
      w2(INwUse0+iaS) = dble(NSsize+NCsize)       !Table-set size
      w2(ISlstT0+iaS) = dble(NSsize)              !IS this table
C--   Return address
      iwsTClone = iaT
C--   Update header previous table
      if(NBskip.eq.0) return
      w2(NFTabl0+iaL) = iaT - iaL                 !Fskip table

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void sws_TbCopy(double *w1, int ia1,
CXXHDR                    double *w2, int ia2, int itag);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_tbcopy FC_FUNC(sws_tbcopy,SWS_TBCOPY)
CXXHFW    void fsws_tbcopy(double*, int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_TbCopy(double *w1, int ia1,
CXXWRP                    double *w2, int ia2, int itag)
CXXWRP    {
CXXWRP     int ja1    = iaCtoF(ia1);
CXXWRP     int ja2    = iaCtoF(ia2);
CXXWRP     fsws_tbcopy(w1, &ja1, w2, &ja2, &itag);
CXXWRP    }
C-----------------------------------------------------------------------

C     =========================================
      subroutine sws_TbCopy(w1,ia1,w2,ia2,itag)
C     =========================================

C--   Copy table content and also yes(no) tags for itag = 1(0)

C--   Author: Michiel Botje h24@nikhef.nl   09-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w1(*), w2(*)

C--   Check input
      if(int(w1(iCword0+iroot0)).ne.iCWorkspace0) stop
     +                   'WSTORE:SWS_TBCOPY: W1 is not a workspace'
      if(int(w2(iCword0+iroot0)).ne.iCWorkspace0) stop
     +                   'WSTORE:SWS_TBCOPY: W2 is not a workspace'
      if(ia1.le.0 .or. ia1.gt.int(w1(INwUse0+iroot0))) stop
     +                   'WSTORE:SWS_TBCOPY: IA1 out of range'
      if(ia2.le.0 .or. ia2.gt.int(w2(INwUse0+iroot0))) stop
     +                   'WSTORE:SWS_TBCOPY: IA2 out of range'
      if(int(w1(iCword0+ia1)).ne.iCWTable0) stop
     +                   'WSTORE:SWS_TBCOPY: source object is not table'
      if(int(w2(iCword0+ia2)).ne.iCWTable0) stop
     +                   'WSTORE:SWS_TBCOPY: target object is not table'
      if(itag.ne.0 .and. itag.ne.1) stop
     +                   'WSTORE:SWS_TBCOPY: itag should be 0 or 1'
C--   Check target = source
      ifw1 = int(w1(IFprnt0+iroot0))
      ifw2 = int(w2(IFprnt0+iroot0))
      if(ifw1.eq.ifw2 .and. ia1.eq.ia2) return
C--   Compare fingerprints
      ift1 = int(w1(IFprnt0+ia1))
      ift2 = int(w2(IFprnt0+ia2))
      if(ift1.ne.ift2) stop
     +  'WSTORE:SWS_TBCOPY: source and target fingerprints do not match'
C--   Copy table body
      it1 = int(w1(ItBody0+ia1))
      it2 = int(w1(ItBend0+ia1))
      do it = it1,it2
        w2(ia2+it) = w1(ia1+it)
      enddo
      if(itag.eq.0) return
C--   Copy tag field
      is1 = ia1 + int(w1(NBTset0+ia1))
      is2 = ia2 + int(w2(NBTset0+ia2))
      nh1 = int(w1(INhead0+is1))
      nh2 = int(w2(INhead0+is2))
      if(nh1.ne.nh2) stop 'WSTORE:SWS_TBCOPY: different header size'
      nt1 = int(w1(INtags0+is1))
      nt2 = int(w2(INtags0+is2))
      if(nt1.ne.nt2) stop 'WSTORE:SWS_TBCOPY: different tag-field size'
      ja1 = ia1+nh1
      ja2 = ia2+nh2
      do j = 0,nt1
        w2(ja2+j) = w1(ja1+j)
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void sws_WsWipe(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_wswipe FC_FUNC(sws_wswipe,SWS_WSWIPE)
CXXHFW    void fsws_wswipe(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_WsWipe(double *w, int ia)
CXXWRP    {
CXXWRP     int ja     = iaCtoF(ia);
CXXWRP     fsws_wswipe(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===========================
      subroutine sws_WsWipe(w,ia)
C     ===========================

C--   Wipe workspace starting at object ia

C--   Author: Michiel Botje h24@nikhef.nl   29-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*), ival(2)

      iaR = iroot0                                         !root address
      iaN = int(w(INwUse0+iaR))                          !last word used

C--   Check input
      if(int(w(iCword0+iaR)).ne.iCWorkspace0) stop
     +                        'WSTORE:SWS_WSWIPE: W is not a workspace'
      if(ia.le.0 .or. ia.gt.iaN) stop
     +                        'WSTORE:SWS_WSWIPE: IA out of range'

      iaS    = int(w(NFTset0+iaR)) + iaR          !first tbset
      iaT    = int(w(NFTabl0+iaR)) + iaR          !first table
      nwmax  = int(w(INwMax0+iaR))                !total size
      nhskip = int(w(IHskip0+iaS))                !hskip
      nhead  = int(w(INhead0+iaS))                !header size
      ntags  = int(w(INtags0+iaS))                !tag size

C-(1) if ia = root, first set or first table then restore WSINIT status
      if(ia.eq.iaR .or. ia.eq.iaS .or. ia.eq.iaT) then
        ival(1)        = nhead
        ival(2)        = ntags
        ihash          = 0
        ihash          = imb_ihash(ihash,ival,2)  !tbset initial fprint
        nsize          = nhskip + nhskip          !nused after wipe
C--     Update root header
        w(NFTabl0+iaR) = 0.D0                     !Fskip table
        w(NFTset0+iaR) = dble(nhskip)             !Fskip tbset
        w(NObjec0+iaR) = 1.D0                     !Number of objects
        w(IWlstS0+iaR) = dble(nhskip)             !IW last tbset in w
        w(IWlstT0+iaR) = dble(nhskip+nhskip)      !IW last table in w
        w(INwUse0+iaR) = dble(nsize)              !Words used
C--     Update 1st table-set header
        w(NFTabl0+iaS) = 0.D0                     !Fskip table
        w(NBTabl0+iaS) = 0.D0                     !Bskip table
        w(NFTset0+iaS) = 0.D0                     !Fskip tset
        w(NBTset0+iaS) = 0.D0                     !Bskip tset
        w(IFprnt0+iaS) = dble(ihash)              !Fingerprint
        w(NObjec0+iaS) = 0.D0                     !number of tables
        w(INwUse0+iaS) = dble(nhskip)             !Oject  size
        w(ISlstT0+iaS) = dble(nhskip)             !IS first table in ts
C--     Wipe
        call smb_Vfill(w(nsize+1),nwmax-nsize,0.D0)
        return
      endif

C-(2) Wipe starting with table-set
      if(int(w(ICword0+ia)).eq.iCWTableSet0)    then
        iaS    = ia                             !IA current  tbset
        iaL    = int(w(NBTset0+iaS)) + iaS      !IA previous tbset
C--     Update header of tables in previous table-set
        nskip  = int(w(NFTabl0+iaL))            !Dist to table1 in iaL
        iaT    = iaL
        do while(nskip.ne.0)                    !Loop over tables in iaL
          iaT            = iaT + nskip          !IA of table
          w(NFTset0+iaT) = 0.D0                 !Set SFskip to zero
          nskip          = int(w(NFTabl0+iaT))  !Distance to next table
        enddo
C--     Update header of previous table-set
        w(NFTset0+iaL)   = 0.D0                 !Set SFskip to zero
C--     Update root header
        ntsets = int(w(IObjec0+iaS)) - 1        !Nb of sets after wipe
        nwsize = iaS - 1                        !Nused after wipe
        w(NObjec0+iaR) = dble(ntsets)           !Number of objects
        w(IWlstS0+iaR) = dble(iaL-iaR)          !IW last tbset in w
        w(IWlstT0+iaR) = dble(iaT-iaR)          !IW last table in w
        w(INwUse0+iaR) = dble(nwsize)           !Words used
C--     Wipe
        call smb_Vfill(w(nwsize+1),nwmax-nwsize,0.D0)
        return
      endif

C-(3) Wipe starting with first table in set
      if(int(w(ICword0+ia)).eq.iCWTable0 .and.
     +   int(w(IObjec0+ia)).eq.1)               then
        iaS    = int(w(NBTset0+ia)) + ia        !IA current tbset
C--     Update current table-set header
        ival(1) = nhead
        ival(2) = ntags
        ihash   = 0
        ihash   = imb_ihash(ihash,ival,2)
        w(NFTabl0+iaS) = 0.D0                   !Fskip table
        w(NFTset0+iaS) = 0.D0                   !Fskip tset
        w(IFprnt0+iaS) = dble(ihash)            !Fingerprint
        w(NObjec0+iaS) = 0.D0                   !number of tables
        w(INwUse0+iaS) = dble(nhskip)           !Oject  size
        w(ISlstT0+iaS) = dble(nhskip)           !IS first table in ts
C--     Update root header
        ntsets = int(w(IObjec0+iaS))            !Nb of sets after wipe
        nwsize = ia - 1                         !Nused after wipe
        iaT    = int(w(NBTabl0+iaS)) + iaS      !IA last table in w
        w(NObjec0+iaR) = dble(ntsets)           !Number of objects
        w(IWlstS0+iaR) = dble(iaS-iaR)          !IW last tbset in w
        w(IWlstT0+iaR) = dble(iaS+nhskip-iaR)   !IW last table in w
        w(INwUse0+iaR) = dble(nwsize)           !Words used
C--     Wipe
        call smb_Vfill(w(nwsize+1),nwmax-nwsize,0.D0)
        return
      endif

C-(4) Wipe starting with table (not first table in set)
      if(int(w(ICword0+ia)).eq.iCWTable0)       then
        iaT    = ia                             !IA current  table
        iaS    = int(w(NBTset0+iaT)) + iaT      !IA current tbset
C--     Update header of previous tables in this set
        nskip            = int(w(NBTabl0+iaT))  !Dist to previous table
        iaL              = iaT + nskip          !IA previous table
        w(NFTabl0+iaL)   = 0.D0                 !TFskip previous table
        iaL              = iaT
        do while(nskip.ne.0)                    !Loop over tables in iaS
          iaL            = iaL + nskip          !IA of table
          w(NFTset0+iaL) = 0.D0                 !Set SFskip to zero
          nskip          = int(w(NBTabl0+iaL))  !Dist to previous table
        enddo
C--     Update header of current table-set
        ival(1) = nhead
        ival(2) = ntags
        ihash   = 0
        ihash   = imb_ihash(ihash,ival,2)
        nskip = int(w(NFTabl0+iaS))             !Dist first table in iaS
        iaT   = iaS
        ntabs = 0
        do while(nskip.ne.0)                    !Loop over tables in iaS
          ntabs   = ntabs + 1                   !Number of tables
          iaT     = iaT + nskip                 !IA of table
          ifp     = int(w(IFprnt0+iaT))         !Fingerprint of table
          ival(1) = ifp
          ihash   = imb_ihash(ihash,ival,1)     !Update fprint of iaS
          nskip   = int(w(NFTabl0+iaT))         !Dist to next table
        enddo
        w(NFTset0+iaS) = 0.D0                   !Set SFskip to zero
        w(IFprnt0+iaS) = dble(ihash)            !Fingerprint of iaS
        w(NObjec0+iaS) = dble(ntabs)            !Number of tables in iaS
        w(INwUse0+iaS) = dble(ia-iaS)           !Size of iaS
        w(ISlstT0+iaS) = dble(iaT-iaS)          !IS last table
C--     Update root header
        ntsets = int(w(IObjec0+iaS))            !Nb of sets after wipe
        nwsize = ia - 1                         !Nused after wipe
        w(NObjec0+iaR) = dble(ntsets)           !Number of objects
        w(IWlstS0+iaR) = dble(iaS-iaR)          !IW last tbset in w
        w(IWlstT0+iaR) = dble(iaT-iaR)          !IW last table in w
        w(INwUse0+iaR) = dble(nwsize)           !Words used
C--     Wipe
        call smb_Vfill(w(nwsize+1),nwmax-nwsize,0.D0)
        return
      endif

C-(5) No valid input address
      stop 'WSTORE:SWS_WSWIPE: IA not root, table-set or table address'

      end

C-----------------------------------------------------------------------
CXXHDR    void sws_TsDump(string fname, int key,
CXXHDR                    double *w, int ia, int &ierr);
C-----------------------------------------------------------------------
CXXHFW  #define fswstsdumpcpp FC_FUNC(swstsdumpcpp,SWSTSDUMPCPP)
CXXHFW    void fswstsdumpcpp(char*, int*, int*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_TsDump(string fname, int key,
CXXWRP                    double *w, int ia, int &ierr)
CXXWRP    {
CXXWRP     int ls = fname.size();
CXXWRP     char *cfname = new char[ls+1];
CXXWRP     strcpy(cfname,fname.c_str());
CXXWRP     int ja     = iaCtoF(ia);
CXXWRP     fswstsdumpcpp(cfname, &ls, &key, w, &ja, &ierr);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================================
      subroutine swsTsDumpCPP(fname,ls,key,w,ia,ierr)
C     ===============================================

      implicit double precision (a-h,o-z)

      character*(100) fname
      dimension w(*)

      if(ls.gt.100) stop
     +             'WSTORE:SWS_TSDUMP: input file name > 100 characters'

      call sws_TsDump(fname(1:ls),key,w,ia,ierr)

      return
      end

C     ==========================================
      subroutine sws_TsDump(fname,key,w,ia,ierr)
C     ==========================================

C--   Author: Michiel Botje h24@nikhef.nl   04-03-20

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      character*(*) fname
      dimension w(*)

C--   Check input
      if(int(w(ICword0+iroot0)).ne.iCWorkspace0) stop
     +            'WSTORE:SWS_TSDUMP: W is not a workspace'
      if(ia.le.0 .or. ia.gt.int(w(INwUse0+iroot0))) stop
     +            'WSTORE:SWS_TSDUMP: IA out of range'
      if(int(w(ICword0+ia)).ne.iCWTableSet0) stop
     +            'WSTORE:SWS_TSDUMP: object to dump is not a table-set'
      if(int(w(NObjec0+ia)).eq.0) stop
     +            'WSTORE:SWS_TSDUMP: cannot dump empty table-set'

C--   Open file
      lun = imb_nextL(0)
      if(lun.eq.0) stop
     +            'WSTORE:SWS_TSDUMP: no logical unit number available'
      open(unit=lun,file=fname,form='unformatted',
     +     status='unknown',err=500)
      ierr = 0
      ivers  = int(w(ICword0+iroot0))
      nhead  = int(w(INhead0+ia))
      ntags  = int(w(INtags0+ia))
      ifprnt = int(w(IFprnt0+ia))
      nwords = int(w(INwUse0+ia))
      ntabs  = int(w(NObjec0+ia))
      ifskip = int(w(NFTabl0+ia))
      i1     = ia + nhead
      i2     = ia + nwords-1
      nw     = i2-i1+1
C--   Write header info
      write(lun,err=500) key,ivers,nhead,ntags,ifprnt,ntabs,nw
C--   Write the table-set but without the header
      write(lun,err=500) (w(i),i=i1,i2)
      close(lun)
C--   Debug print
*      write(6,'(/'' SWS_TSDUMP: key   = '',I10)') key
*      write(6,'( ''             ivers = '',I10)') ivers
*      write(6,'( ''             nhead = '',I10)') nhead
*      write(6,'( ''             ntags = '',I10)') ntags
*      write(6,'( ''             fprnt = '',I10)') ifprnt
*      write(6,'( ''             nwtab = '',I10)') nw
*      write(6,'( ''             ntabs = '',I10)') ntabs
*      write(6,'( ''             ia1   = '',I10)') i1
*      write(6,'( ''             ia2   = '',I10)') i2

      return

  500 continue
C--   Open or write error
      ierr = -1
      return

      end

C-----------------------------------------------------------------------
CXXHDR    int iws_TsRead(string fname, int key, double *w, int &ierr);
C-----------------------------------------------------------------------
CXXHFW  #define fiwstsreadcpp FC_FUNC(iwstsreadcpp,IWSTSREADCPP)
CXXHFW    int fiwstsreadcpp(char*, int*, int*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_TsRead(string fname, int key, double *w, int &ierr)
CXXWRP    {
CXXWRP     int ls = fname.size();
CXXWRP     char *cfname = new char[ls+1];
CXXWRP     strcpy(cfname,fname.c_str());
CXXWRP     int ia = fiwstsreadcpp(cfname, &ls, &key, w, &ierr);
CXXWRP     return iaFtoC(ia);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================================
      integer function iwsTsReadCPP(fname,ls,key,w,ierr)
C     ===================================================

      implicit double precision (a-h,o-z)

      character*(100) fname
      dimension w(*)

      if(ls.gt.100) stop
     +             'WSTORE:IWS_TSREAD: input file name > 100 characters'

      iwsTsReadCPP = iws_TsRead(fname(1:ls),key,w,ierr)

      return
      end

C     =============================================
      integer function iws_TsRead(fname,key,w,ierr)
C     =============================================

C--   Author: Michiel Botje h24@nikhef.nl   04-03-20

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      character*(*) fname
      dimension w(*)

C--   Init
      ierr = 0
      iws_TsRead = 0
C--   Check input
      if(int(w(ICword0+iroot0)).ne.iCWorkspace0) stop
     +             'WSTORE:IWS_TSREAD: W is not a workspace'
C--   Open file
      lun = imb_nextL(0)
      if(lun.eq.0) stop
     +             'WSTORE:IWS_TSREAD: no logical unit number available'
      open(unit=lun,file=fname,form='unformatted',status='old',err=500)
C--   Read header info
      read(lun,err=500,end=500) keyr,ivers,nhead,ntags,ifprnt,ntabs,nw
      iaR    = iroot0                             !root address
      iaL    = int(w(IWlstS0+iaR)) + iaR          !IA current table-set
C--   Check key,ivers,nhead,ntags
      if(key.ne.0 .and. key.ne.keyr)       then
        ierr = -2
        return
      elseif(int(w(ICword0+iaR)).ne.ivers) then
        ierr = -2
        return
      elseif(int(w(INhead0+iaL)).ne.nhead) then
        ierr = -2
        return
      elseif(int(w(INtags0+iaL)).ne.ntags) then
        ierr = -2
        return
      endif
C--   Check workspace size
      NHskip = nhead + ntags
      NWused = int(w(INwUse0+iaR))                !words used
      Nspace = int(w(INwMax0+iaR))                !total size
      if( iwsEtrailer(w) .eq. 1 ) then
        NWneed = NWused + nw + 1                  !words needed
      else
        NWneed = NWused + NHskip + nw + 1         !words needed
      endif
C--   Run out-of-space
      if(NWneed.gt.Nspace) call swsWsEmsg(w,NWneed,'WSTORE:IWS_TSREAD')
C--   New table-set in w
      iaS    = iws_NewSet(w)
C--   Get nw-used but subtract ntags since tagfield will be overwritten
      NWused = int(w(INwUse0+iaR))-ntags
C--   Read table-set without header
      i1     = iaS + nhead
      i2     = i1  + nw-1
      read(lun,err=500,end=500) (w(i),i=i1,i2)
C--   Set distance to root in all table headers
      nfskip = NHskip
      iaT = iaS + nfskip
      do while(nfskip.ne.0)
        w(IwAddr0+iaT) = dble(iaT-iaR)
        nfskip         = int(w(NFTabl0+iaT))
        iaT            = iaT + nfskip
      enddo
C--   Update Wspace Header
      w(IWlstS0+iaR) = dble(iaS-iaR)              !IW last table-set
      w(IWlstT0+iaR) = dble(iaT-iaR)              !IW last table
      w(INwUse0+iaR) = dble(NWused+nw)            !Words used
      w(NFTabl0+iaR) = dble(NHskip+NHskip)        !Fskip table
C--   Update Table-set Header
      w(NFTabl0+iaS) = dble(NHskip)               !Fskip table
      w(Ifprnt0+iaS) = dble(ifprnt)               !fingerprint
      w(NObjec0+iaS) = dble(ntabs)                !Number of objects
      w(INwUse0+iaS) = dble(nhead+nw)             !Table-set size
      w(ISlstT0+iaS) = dble(iaT-iaS)              !IS last table
C--   Return table set address
      iws_TsRead = iaS
C--   Close file
      close(lun)

      return

  500 continue
C--   Open or read error
      ierr = -1
      return

      end

C-----------------------------------------------------------------------
CXXHDR    void sws_WsMark(int &mws, int &mset, int &mtab );
C-----------------------------------------------------------------------
CXXHFW  #define fsws_wsmark FC_FUNC(sws_wsmark,SWS_WSMARK)
CXXHFW    void fsws_wsmark(int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_WsMark(int &mws, int &mset, int &mtab)
CXXWRP    {
CXXWRP     fsws_wsmark(&mws, &mset, &mtab);
CXXWRP    }
C-----------------------------------------------------------------------

C     ======================================
      subroutine sws_WsMark(mws, mset, mtab)
C     ======================================

C--   Return current markers of root, table-set and table

C--   Author: Michiel Botje h24@nikhef.nl   30-01-21

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      mws  = iCWorkspace0
      mset = iCWTableSet0
      mtab = iCWTable0

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_HdSize();
C-----------------------------------------------------------------------
CXXHFW  #define fiws_hdsize FC_FUNC(iws_hdsize,IWS_HDSIZE)
CXXHFW    int fiws_hdsize();
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_HdSize()
CXXWRP    {
CXXWRP     return fiws_hdsize();
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================
      integer function iws_HdSize()
C     =============================

C--   Return header size

C--   Author: Michiel Botje h24@nikhef.nl   02-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      iws_HdSize = nwHeader0

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_TbSize(int *imi, int *ima, int ndim);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_tbsize FC_FUNC(iws_tbsize,IWS_TBSIZE)
CXXHFW    int fiws_tbsize(int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_TbSize(int *imi, int *ima, int ndim)
CXXWRP    {
CXXWRP     return fiws_tbsize(imi, ima, &ndim);
CXXWRP    }
C-----------------------------------------------------------------------

C     =========================================
      integer function iws_TbSize(imi,ima,ndim)
C     =========================================

C--   Table size (body + metadata)
C--
C--   imi, ima  (in): index limits
C--   ndim      (in): number of dimensions

C--   Author: Michiel Botje h24@nikhef.nl   02-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension imi(*), ima(*)

      if(ndim.lt.1 .or. ndim.gt.mdim0) stop
     +                      'WSTORE:IWS_TBSIZE: ndim out of range'
      nw = 1
      do i = 1,ndim
        if(imi(i).gt.ima(i)) stop
     +                      'WSTORE:IWS_TBSIZE: imin > imax encountered'
        nw = nw*(ima(i)-imi(i)+1)
      enddo

      iws_TbSize = nw + 3*ndim + 2

      return
      end

C=======================================================================
C==   Pointer functions  ===============================================
C=======================================================================

C     =========================================
      subroutine swsGetMeta(w,ia,nd,kk,imi,ima)
C     =========================================

C--   Extract metadata from a table
C--
C--   w            (in): workspace
C--   ia           (in): table address
C--   nd          (out): number of dimensions
C--   kk(nd+2)    (out): fingerprint + pointer coefficients
C--   imi,ima(nd) (out): index limits

C--   Author: Michiel Botje h24@nikhef.nl   27-11-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*), kk(*), imi(*), ima(*)

      mm     = int(w(ItMeta0+ia))+ia
      nd     = int(w(mm))
      mi     = mm+nd+1
      ma     = mi+nd
      kk(1)  = int(w(Ifprnt0+ia))
      kk(2)  = int(w(mm+1))
      do i = 1,nd
        kk(i+2) = int(w(mm+1+i))
        imi(i)  = int(w(mi+i))
        ima(i)  = int(w(ma+i))
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_Tpoint(double *w, int ia, int* index, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_tpoint FC_FUNC(iws_tpoint,IWS_TPOINT)
CXXHFW    int fiws_tpoint(double*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_Tpoint(double *w, int ia, int *index, int n)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     int iadr  = fiws_tpoint(w, &ja, index, &n);
CXXWRP     return iaFtoC(iadr);
CXXWRP    }
C-----------------------------------------------------------------------

C     =========================================
      integer function iws_Tpoint(w,ia,index,n)
C     =========================================

C--   Pointer function
C--
C--   w, ia     (in): Table ia in workspace w
C--   index     (in): Array with ndim index values
C--   n         (in): Dimension of index as set in the calling routine

C--   Author: Michiel Botje h24@nikhef.nl   02-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*), index(*)
      save nd, kk, imi, ima
      dimension kk(mdim0+2), imi(mdim0), ima(mdim0)

C--   Check input
      if(n.le.0) stop 'WSTORE:IWS_TPOINT: n <= 0'
      if(int(w(iroot0)).ne.iCWorkspace0) stop
     +                'WSTORE:IWS_TPOINT: W is not a workspace'
      if(ia.le.0 .or. ia.gt.int(w(INwUse0+iroot0))) stop
     +                'WSTORE:IWS_TPOINT: IA out of range'
      if(int(w(ia)).ne.iCWTable0) stop
     +                'WSTORE:IWS_TPOINT: IA is not a table address'
C--   Fingerprint
      ifp = int(w(Ifprnt0+ia))
C--   Load metadata, if needed
      if(kk(1).ne.ifp) call swsGetMeta(w,ia,nd,kk,imi,ima)
      if(n.lt.nd) stop        'WSTORE:IWS_TPOINT: n < ndim of table'
C--   Check index range and build pointer
      ipoint = kk(2) + ia
      do i = 1,nd
        if(index(i).lt.imi(i) .or. index(i).gt.ima(i)) then
          write(6,
     +        '(''WSTORE:IWS_TPOINT: index '',I3,'' out of range'')') i
          stop
        endif
        ipoint = ipoint + kk(i+2)*index(i)
      enddo
C--   Check pointer
      ibot = int(w(ItBody0+ia))+ia
      ieot = int(w(ItBend0+ia))+ia
      if(ipoint.lt.ibot .or. ipoint.gt.ieot)
     +   stop 'WSTORE:IWS_TPOINT: calculated pointer outside table body'

      iws_Tpoint = ipoint

      return
      end

C=======================================================================
C==   Access to information ============================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    int iws_IsaWorkspace(double *w);
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  #define fiws_isaworkspace FC_FUNC(iws_isaworkspace,IWS_ISAWORKSPACE)
CXXHFW    int fiws_isaworkspace(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IsaWorkspace(double *w)
CXXWRP    {
CXXWRP     return fiws_isaworkspace(w);
CXXWRP    }
C-----------------------------------------------------------------------

C     ====================================
      integer function iws_IsaWorkspace(w)
C     ====================================

C--   Return 1 (0) if w is (not) a workspace

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      if(int(w(iCword0+iroot0)).eq.iCWorkspace0) then
        iws_IsaWorkspace = 1
      else
        iws_IsaWorkspace = 0
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_SizeOfW(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_sizeofw FC_FUNC(iws_sizeofw,IWS_SIZEOFW)
CXXHFW    int fiws_sizeofw(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_SizeOfW(double *w)
CXXWRP    {
CXXWRP     return fiws_sizeofw(w);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================
      integer function iws_SizeOfW(w)
C     ===============================

C--   Total size of w

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_SizeOfW = int(w(iNwMax0+iroot0))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_WordsUsed(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_wordsused FC_FUNC(iws_wordsused,IWS_WORDSUSED)
CXXHFW    int fiws_wordsused(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_WordsUsed(double *w)
CXXWRP    {
CXXWRP     return fiws_wordsused(w);
CXXWRP    }
C-----------------------------------------------------------------------

C     =================================
      integer function iws_WordsUsed(w)
C     =================================

C--   Words used

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_WordsUsed = int(w(INwUse0+iroot0))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_Nheader(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_nheader FC_FUNC(iws_nheader,IWS_NHEADER)
CXXHFW    int fiws_nheader(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_Nheader(double *w)
CXXWRP    {
CXXWRP     return fiws_nheader(w);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================
      integer function iws_Nheader(w)
C     ===============================

C--   Header size

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iaS         = int(w(IWlstS0+iroot0))+iroot0
      iws_Nheader = int(w(INhead0+iaS))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_Ntags(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_ntags FC_FUNC(iws_ntags,IWS_NTAGS)
CXXHFW    int fiws_ntags(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_Ntags(double *w)
CXXWRP    {
CXXWRP     return fiws_ntags(w);
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================
      integer function iws_Ntags(w)
C     =============================

C--   Tag field size

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iaS       = int(w(IWlstS0+iroot0))+iroot0
      iws_Ntags = int(w(INtags0+iaS))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_HeadSkip(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_headskip FC_FUNC(iws_headskip,IWS_HEADSKIP)
CXXHFW    int fiws_headskip(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_HeadSkip(double *w)
CXXWRP    {
CXXWRP     return fiws_headskip(w);
CXXWRP    }
C-----------------------------------------------------------------------

C     ================================
      integer function iws_HeadSkip(w)
C     ================================

C--   Header + tag field size

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iaS          = int(w(IWlstS0+iroot0))+iroot0
      iws_HeadSkip = int(w(IHskip0+iaS))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaRoot();
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iaroot FC_FUNC(iws_iaroot,IWS_IAROOT)
CXXHFW    int fiws_iaroot();
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaRoot()
CXXWRP    {
CXXWRP     int ja    = fiws_iaroot();
CXXWRP     return iaFtoC(ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     =============================
      integer function iws_IaRoot()
C     =============================

C--   Root address

C--   Author: Michiel Botje h24@nikhef.nl   02-03-20

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      iws_IaRoot = iroot0

      return
      end


C-----------------------------------------------------------------------
CXXHDR    int iws_IaDrain(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iadrain FC_FUNC(iws_iadrain,IWS_IADRAIN)
CXXHFW    int fiws_iadrain(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaDrain(double *w)
CXXWRP    {
CXXWRP     int ja    = fiws_iadrain(w);
CXXWRP     return iaFtoC(ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===============================
      integer function iws_IaDrain(w)
C     ===============================

C--   Drain-word address

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      dum         = w(1)                         !avoid compiler warning
      iws_IaDrain = IDrain0 + iroot0

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaNull(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_ianull FC_FUNC(iws_ianull,IWS_IANULL)
CXXHFW    int fiws_ianull(double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaNull(double *w)
CXXWRP    {
CXXWRP     int ja    = fiws_ianull(w);
CXXWRP     return iaFtoC(ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     ==============================
      integer function iws_IaNull(w)
C     ==============================

C--   Null-word address

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      dum        = w(1)                          !avoid compiler warning
      iws_IaNull = IZnull0+iroot0

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_ObjectType(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_objecttype FC_FUNC(iws_objecttype,IWS_OBJECTTYPE)
CXXHFW    int fiws_objecttype(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_ObjectType(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_objecttype(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     =====================================
      integer function iws_ObjectType(w,ia)
C     =====================================

C--   Return object type

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_ObjectType = 0

      if(int(w(ia)).eq.iCWorkspace0)      then
        iws_ObjectType = 1
      elseif(int(w(ia)).eq.iCWTableSet0) then
        iws_ObjectType = 2
      elseif(int(w(ia)).eq.iCWTable0)    then
        iws_ObjectType = 3
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_ObjectSize(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_objectsize FC_FUNC(iws_objectsize,IWS_OBJECTSIZE)
CXXHFW    int fiws_objectsize(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_ObjectSize(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_objectsize(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     =====================================
      integer function iws_ObjectSize(w,ia)
C     =====================================

C--   Return object size

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_ObjectSize = int(w(INwUse0+ia))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_Nobjects(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_nobjects FC_FUNC(iws_nobjects,IWS_NOBJECTS)
CXXHFW    int fiws_nobjects(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_Nobjects(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_nobjects(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      integer function iws_Nobjects(w,ia)
C     ===================================

C--   Return object size

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_Nobjects = int(w(NObjec0+ia))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_ObjectNumber(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_objectnumber FC_FUNC(iws_objectnumber,IWS_OBJECTNUMBER)
CXXHFW    int fiws_objectnumber(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_ObjectNumber(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_objectnumber(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     =======================================
      integer function iws_ObjectNumber(w,ia)
C     =======================================

C--   Return object size

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_ObjectNumber = int(w(IObjec0+ia))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_FingerPrint(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_fingerprint FC_FUNC(iws_fingerprint,IWS_FINGERPRINT)
CXXHFW    int fiws_fingerprint(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_FingerPrint(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_fingerprint(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     ======================================
      integer function iws_FingerPrint(w,ia)
C     ======================================

C--   Return object size

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_FingerPrint = int(w(IFprnt0+ia))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaFirstTag(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iafirsttag FC_FUNC(iws_iafirsttag,IWS_IAFIRSTTAG)
CXXHFW    int fiws_iafirsttag(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaFirstTag(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     int ift   = fiws_iafirsttag(w, &ja);
CXXWRP     return iaFtoC(ift);
CXXWRP    }
C-----------------------------------------------------------------------

C     =====================================
      integer function iws_IaFirstTag(w,ia)
C     =====================================

C--   Null-word address

C--   Author: Michiel Botje h24@nikhef.nl   03-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      dum            = w(1)                      !avoid compiler warning
      iws_IaFirstTag = ia + nwHeader0

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_TableDim(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_tabledim FC_FUNC(iws_tabledim,IWS_TABLEDIM)
CXXHFW    int fiws_tabledim(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_TableDim(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_tabledim(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      integer function iws_TableDim(w,ia)
C     ===================================

C--   Number of dimensions of table ia

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      ja           = int(w(itMeta0+ia))+ia
      iws_TableDim = int(w(ja))

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaKARRAY(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iakarray FC_FUNC(iws_iakarray,IWS_IAKARRAY)
CXXHFW    int fiws_iakarray(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaKARRAY(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     int iak   = fiws_iakarray(w, &ja);
CXXWRP     return iaFtoC(iak);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      integer function iws_IaKARRAY(w,ia)
C     ===================================

C--   First word of KARRAY

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_IaKARRAY = int(w(itMeta0+ia))+ia+1

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaIMIN(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iaimin FC_FUNC(iws_iaimin ,IWS_IAIMIN)
CXXHFW    int fiws_iaimin(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaIMIN(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     int iak   = fiws_iaimin(w, &ja);
CXXWRP     return iaFtoC(iak);
CXXWRP    }
C-----------------------------------------------------------------------

C     =================================
      integer function iws_IaIMIN(w,ia)
C     =================================

C--   First word of IMIN

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_IaIMIN = int(w(itImin0+ia))+ia

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_IaIMAX(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_iaimax FC_FUNC(iws_iaimax ,IWS_IAIMAX)
CXXHFW    int fiws_iaimax(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_IaIMAX(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     int iak   = fiws_iaimax(w, &ja);
CXXWRP     return iaFtoC(iak);
CXXWRP    }
C-----------------------------------------------------------------------

C     =================================
      integer function iws_IaIMAX(w,ia)
C     =================================

C--   First word of IMAX

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_IaIMAX = int(w(itImax0+ia))+ia

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_BeginTbody(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_begintbody FC_FUNC(iws_begintbody ,IWS_BEGINTBODY)
CXXHFW    int fiws_begintbody(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_BeginTbody(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     int iak   = fiws_begintbody(w, &ja);
CXXWRP     return iaFtoC(iak);
CXXWRP    }
C-----------------------------------------------------------------------

C     =====================================
      integer function iws_BeginTbody(w,ia)
C     =====================================

C--   First word of table body

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_BeginTbody = int(w(itBody0+ia))+ia

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_EndTbody(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_endtbody FC_FUNC(iws_endtbody ,IWS_ENDTBODY)
CXXHFW    int fiws_endtbody(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_EndTbody(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     int iak   = fiws_endtbody(w, &ja);
CXXWRP     return iaFtoC(iak);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===================================
      integer function iws_EndTbody(w,ia)
C     ===================================

C--   Last word of table body

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      iws_EndTbody = int(w(itBend0+ia))+ia

      return
      end

C=======================================================================
C==   Navigation =======================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    int iws_TFskip(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_tfskip FC_FUNC(iws_tfskip ,IWS_TFSKIP)
CXXHFW    int fiws_tfskip(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_TFskip(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_tfskip(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     =================================
      integer function iws_TFskip(w,ia)
C     =================================

C--   Distance (signed) to next table marker

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

C--   Check input
      if(int(w(iroot0)).ne.iCWorkspace0) stop
     +                        'WSTORE:IWS_TFSKIP: W is not a workspace'
      if(ia.le.0 .or. ia.gt.int(w(INwUse0+iroot0))) stop
     +                        'WSTORE:IWS_TFSKIP: IA out of range'
C--   Go
      mark = int(w(ia))
      if(mark.eq.iCWorkspace0  .or.
     +   mark.eq.iCWTableSet0 .or. mark.eq.iCWTable0) then
        iws_TFskip = int(w(NFTabl0+ia))
      else
        iws_TFskip = 0
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_TBskip(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_tbskip FC_FUNC(iws_tbskip ,IWS_TBSKIP)
CXXHFW    int fiws_tbskip(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_TBskip(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_tbskip(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     =================================
      integer function iws_TBskip(w,ia)
C     =================================

C--   Distance (signed) to previous table marker

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

C--   Check input
      if(int(w(iroot0)).ne.iCWorkspace0) stop
     +                  'WSTORE:IWS_TBSKIP: W is not a workspace'
      if(ia.le.0) stop  'WSTORE:IWS_TBSKIP: IA out of range'

C--   Special case ia > last word used
      nwused = int(w(INwUse0+iroot0))
      if(ia.gt.nwused) then
        iaLset = int(w(IWlstS0+iroot0)) + iroot0
        ntabs  = int(w(NObjec0+iaLset))
        if(ntabs.eq.0) then
          iws_TBskip = 0                          !no table in last tset
        else
          iaLtab     = int(w(IWlstT0+iroot0)) + iroot0
          iws_TBskip = iaLtab - ia
        endif
        return
      endif

C--   Go for ia in range
      mark = int(w(ia))
      if(mark.eq.iCWorkspace0 .or.
     +   mark.eq.iCWTableSet0 .or. mark.eq.iCWTable0) then
        iws_TBskip = int(w(NBTabl0+ia))
      else
        iws_TBskip = 0
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_SFskip(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_sfskip FC_FUNC(iws_sfskip ,IWS_SFSKIP)
CXXHFW    int fiws_sfskip(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_SFskip(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_sfskip(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     =================================
      integer function iws_SFskip(w,ia)
C     =================================

C--   Distance (signed) to next table-set marker

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

C--   Check input
      if(int(w(iroot0)).ne.iCWorkspace0) stop
     +                        'WSTORE:IWS_SFSKIP: W is not a workspace'
      if(ia.le.0 .or. ia.gt.int(w(INwUse0+iroot0))) stop
     +                        'WSTORE:IWS_SFSKIP: IA out of range'
C--   Go
      mark = int(w(ia))
      if(mark.eq.iCWorkspace0 .or.
     +   mark.eq.iCWTableSet0 .or. mark.eq.iCWTable0) then
        iws_SFskip = int(w(NFTset0+ia))
      else
        iws_SFskip = 0
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int iws_SBskip(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fiws_sbskip FC_FUNC(iws_sbskip ,IWS_SBSKIP)
CXXHFW    int fiws_sbskip(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int iws_SBskip(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     return fiws_sbskip(w, &ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     =================================
      integer function iws_SBskip(w,ia)
C     =================================

C--   Distance (signed) to previous table-set marker

C--   Author: Michiel Botje h24@nikhef.nl   04-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

C--   Check input
      if(int(w(iroot0)).ne.iCWorkspace0) stop
     +                  'WSTORE:IWS_SBSKIP: W is not a workspace'
      if(ia.le.0) stop  'WSTORE:IWS_SBSKIP: IA out of range'

C--   Special case ia > last word used
      nwused = int(w(INwUse0+iroot0))
      if(ia.gt.nwused) then
        iaLset = int(w(IWlstS0+iroot0)) + iroot0
        iws_SBskip = iaLset - ia
        return
      endif

C--   Go for ia in range
      mark = int(w(ia))
      if(mark.eq.iCWorkspace0 .or.
     +   mark.eq.iCWTableSet0 .or. mark.eq.iCWTable0) then
        iws_SBskip = int(w(NBTset0+ia))
      else
        iws_SBskip = 0
      endif

      return
      end

C=======================================================================
C==  Verbose dumps  ====================================================
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    void sws_WsTree(double *w);
C-----------------------------------------------------------------------
CXXHFW  #define fswswstree FC_FUNC(swswstree ,SMBWSTREE)
CXXHFW    int fswswstree(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_WsTree(double *w)
CXXWRP    {
CXXWRP     int iroot    = 0;
CXXWRP     fswswstree(w,&iroot);
CXXWRP    }
C-----------------------------------------------------------------------

C     ========================
      subroutine sws_WsTree(w)
C     ========================

C--   Print workspace tree

C--   Author: Michiel Botje h24@nikhef.nl   12-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

      call swsWStree(w,1)

      return
      end

C     ==============================
      subroutine  swsWStree(w,iroot)
C     ==============================

C--   The iroot argumant takes care of the address printouts
C--   Fortran iroot = 1 --> print address as is
C--   C++     iroot = 0 --> substract 1 from all addresses

C--   Author: Michiel Botje h24@nikhef.nl   12-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

C--   Check input
      if(int(w(iroot0)).ne.iCWorkspace0) stop
     +                        'WSTORE:SWS_WSTREE: W is not a workspace'
      if(iroot.ne.0 .and. iroot.ne.1) stop
     +                        'WSTORE:SWS_WSTREE: iroot must be 0 or 1'

      ia = 1                                 !address of workspace
      call swswprint(w,ia,iroot)             !print workspace
      id = int(w(NFTset0+ia))                !distance to 1st table-set
      do while(id.ne.0)                      !loop over table-sets
        ia = ia+id                           !address of table-set
        call swssprint(w,ia,iroot)           !print table-set
        jd = int(w(NFTabl0+ia))              !distance to first table
        do while(jd.ne.0)                    !loop over tables
          ia = ia+jd                         !address of table
          call swstprint(w,ia,iroot)         !print table
          jd = int(w(NFTabl0+ia))            !distance to next table
        enddo                                !end loop over tables
        id = int(w(NFTset0+ia))              !distance to next table-set
      enddo                                  !end loop over table-sets

      return
      end

C     ================================
      subroutine swswprint(w,ia,iroot)
C     ================================

C--   The iroot argumant takes care of the address printouts
C--   Fortran iroot = 1 --> print address as is
C--   C++     iroot = 0 --> substract 1 from all addresses

C--   Author: Michiel Botje h24@nikhef.nl   12-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)
      character*15 number, hcode
      character*80 text

C--   Check input address
      if(int(w(ia)).ne.iCWorkspace0) stop
     +                'WSTORE:SMBWPRINT: IA is not a workspace address'

      iadr  = ia-iroot0+iroot
      nwds  = int(w(INwUse0+ia))
      ifpt  = int(w(IFprnt0+ia))
      nobj  = int(w(NObjec0+ia))
      call smb_itoch(nobj,number,leng)
      call smb_hcode(ifpt,hcode)
      call smb_cfill(' ',text)
      text = 'workspace with '//number(1:leng)//' table-sets'
      ltxt = imb_lastc(text)

      write(6,'(/1X,''ADDRESS'',4X,''SIZE'',8X,''FINGERPRINT'',
     +           4X,''OBJECT'')')
      write(6,'(2I8,4X,A15,4X,A)') iadr, nwds, hcode, text(1:ltxt)

      return
      end

C     ================================
      subroutine swssprint(w,ia,iroot)
C     ================================

C--   The iroot argumant takes care of the address printouts
C--   Fortran iroot = 1 --> print address as is
C--   C++     iroot = 0 --> substract 1 from all addresses

C--   Author: Michiel Botje h24@nikhef.nl   12-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)
      character*15 number, hcode
      character*80 text

C--   Check input address
      if(int(w(ia)).ne.iCWTableSet0) stop
     +                'WSTORE:SMBSPRINT: IA is not a table-set address'

      iadr  = ia-iroot0+iroot
      nwds  = int(w(INwUse0+ia))
      ifpt  = int(w(IFprnt0+ia))
      nobj  = int(w(NObjec0+ia))
      call smb_itoch(nobj,number,leng)
      call smb_hcode(ifpt,hcode)
      call smb_cfill(' ',text)
      text = '--- table-set with '//number(1:leng)//' tables'
      ltxt = imb_lastc(text)

      write(6,'(2I8,4X,A15,4X,A)') iadr, nwds, hcode, text(1:ltxt)

      return
      end

C     ================================
      subroutine swstprint(w,ia,iroot)
C     ================================

C--   The iroot argumant takes care of the address printouts
C--   Fortran iroot = 1 --> print address as is
C--   C++     iroot = 0 --> substract 1 from all addresses

C--   Author: Michiel Botje h24@nikhef.nl   12-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)
      character*15 number, hcode
      character*80 text

C--   Check input address
      if(int(w(ia)).ne.iCWTable0) stop
     +                    'WSTORE:SMBTPRINT: IA is not a table address'

      iadr  = ia-iroot0+iroot
      nwds  = int(w(INwUse0+ia))
      ifpt  = int(w(IFprnt0+ia))
      idim  = int(w(ItMeta0+ia))+ia
      ndim  = int(w(idim))
      call smb_itoch(ndim,number,leng)
      call smb_hcode(ifpt,hcode)
      call smb_cfill(' ',text)
      text = '------- table with '//number(1:leng)//' dimensions'
      ltxt = imb_lastc(text)

      write(6,'(2I8,4X,A15,4X,A)') iadr, nwds, hcode, text(1:ltxt)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void sws_WsHead(double *w, int ia);
C-----------------------------------------------------------------------
CXXHFW  #define fsws_wshead FC_FUNC(sws_wshead ,SWS_WSHEAD)
CXXHFW    int fsws_wshead(double*, int*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void sws_WsHead(double *w, int ia)
CXXWRP    {
CXXWRP     int ja    = iaCtoF(ia);
CXXWRP     fsws_wshead(w,&ja);
CXXWRP    }
C-----------------------------------------------------------------------

C     ===========================
      subroutine sws_WsHead(w,ia)
C     ===========================

C--   Print header of workspace, table-set or table

C--   Author: Michiel Botje h24@nikhef.nl   12-12-19

      implicit double precision (a-h,o-z)

      include 'wstore0.inc'

      dimension w(*)

C--   Check input
      if(int(w(iroot0)).ne.iCWorkspace0) stop
     +                        'WSTORE:SWS_WSHEAD: W is not a workspace'
      if(ia.le.0 .or. ia.gt.int(w(INwUse0+iroot0))) stop
     +                        'WSTORE:SWS_WSHEAD: IA out of range'

      if(int(w(ia)).eq.iCWorkspace0) then
         write(6,'(/'' Workspace Header'')')
         write(6,'( '' 0 Cword    '',I15  )') int(w(ICword0+ia))
         write(6,'( '' 1 IW       '',I15  )') int(w(IwAddr0+ia))
         write(6,'( '' 2 TFskip   '',I15  )') int(w(NFTabl0+ia))
         write(6,'( '' 3 TBskip   '',I15  )') int(w(NBTabl0+ia))
         write(6,'( '' 4 SFskip   '',I15  )') int(w(NFTset0+ia))
         write(6,'( '' 5 SBskip   '',I15  )') int(w(NBTset0+ia))
         write(6,'( '' 6 Fprint   '',I15  )') int(w(IFprnt0+ia))
         write(6,'( '' 7 Nobj     '',I15  )') int(w(NObjec0+ia))
         write(6,'( '' 8 Iobj     '',I15  )') int(w(IObjec0+ia))
         write(6,'( '' 9 NWused   '',I15  )') int(w(INwUse0+ia))
         write(6,'( ''10 IW Lset  '',I15  )') int(w(IWlstS0+ia))
         write(6,'( ''11 IW Ltab  '',I15  )') int(w(IWlstT0+ia))
         write(6,'( ''12 NWtotal  '',I15  )') int(w(iNwMax0+ia))
         write(6,'( ''13 Drain    '',E15.5)')     w(IDrain0+ia)
         write(6,'( ''14 Null     '',E15.5)')     w(IZnull0+ia)
      elseif(int(w(ia)).eq.iCWTableSet0) then
         write(6,'(/'' Table-set Header'')')
         write(6,'( '' 0 Cword    '',I15  )') int(w(ICword0+ia))
         write(6,'( '' 1 IW       '',I15  )') int(w(IwAddr0+ia))
         write(6,'( '' 2 TFskip   '',I15  )') int(w(NFTabl0+ia))
         write(6,'( '' 3 TBskip   '',I15  )') int(w(NBTabl0+ia))
         write(6,'( '' 4 SFskip   '',I15  )') int(w(NFTset0+ia))
         write(6,'( '' 5 SBskip   '',I15  )') int(w(NBTset0+ia))
         write(6,'( '' 6 Fprint   '',I15  )') int(w(IFprnt0+ia))
         write(6,'( '' 7 Nobj     '',I15  )') int(w(NObjec0+ia))
         write(6,'( '' 8 Iobj     '',I15  )') int(w(IObjec0+ia))
         write(6,'( '' 9 NWused   '',I15  )') int(w(INwUse0+ia))
         write(6,'( ''10 Hsize    '',I15  )') int(w(INhead0+ia))
         write(6,'( ''11 TagSize  '',I15  )') int(w(INtags0+ia))
         write(6,'( ''12 Hskip    '',I15  )') int(w(IHskip0+ia))
         write(6,'( ''13 IS Ltab  '',I15  )') int(w(ISlstT0+ia))
      elseif(int(w(ia)).eq.iCWTable0) then
         write(6,'(/'' Table Header'')')
         write(6,'( '' 0 Cword    '',I15  )') int(w(ICword0+ia))
         write(6,'( '' 1 IW       '',I15  )') int(w(IwAddr0+ia))
         write(6,'( '' 2 TFskip   '',I15  )') int(w(NFTabl0+ia))
         write(6,'( '' 3 TBskip   '',I15  )') int(w(NBTabl0+ia))
         write(6,'( '' 4 SFskip   '',I15  )') int(w(NFTset0+ia))
         write(6,'( '' 5 SBskip   '',I15  )') int(w(NBTset0+ia))
         write(6,'( '' 6 Fprint   '',I15  )') int(w(IFprnt0+ia))
         write(6,'( '' 7 Nobj     '',I15  )') int(w(NObjec0+ia))
         write(6,'( '' 8 Iobj     '',I15  )') int(w(IObjec0+ia))
         write(6,'( '' 9 NWused   '',I15  )') int(w(INwUse0+ia))
         write(6,'( ''10 IT Meta  '',I15  )') int(w(ItMeta0+ia))
         write(6,'( ''11 IT IMIN  '',I15  )') int(w(ItImin0+ia))
         write(6,'( ''12 IT IMAX  '',I15  )') int(w(ItImax0+ia))
         write(6,'( ''13 IT Bbody '',I15  )') int(w(ItBody0+ia))
         write(6,'( ''14 IT Ebody '',I15  )') int(w(ItBend0+ia))
      else
         stop
     +  'WSTORE:SWS_WSHEAD: IA is not a workspace, table-set or table'
      endif

      return
      end

