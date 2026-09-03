
C--   Package (par) to manage the evolution parameter repository pars8.
C--   All routines are called from the QCDNUM user interface routines.
C--
C--   Arguments: in--lower case; OUT--upper case; InOut--mixed case
C--
C--   ------------------------------------------------------------------
C--   Routines called from user interface
C--   ------------------------------------------------------------------
C--   sparInit(NUSED)
C--   - Initialise the pars8 parameter store
C--   - To be called directly after the x-q grid is defined
C--
C--   sparMakeBase
C--   - Setup base parameters and base tables in slot 1 of pars8
C--   - Called by each user parameter input routine
C--
C--   sparRemakeBase(islot)
C--   - Copy islot to /qpars6/ and remake the base in slot 1 of pars8
C--   - Called by USEPAR and PULLCP
C--
C--   sparParsForFill(w,kset)
C--   - Manage parameters after creation but before filling a pdf set
C--   - Called by EVOLFG, EVSGNS, EXTPDF, USRPDF and EVDGLAP
C--
C--   sparParsForCopy(w1,kset1,w2,kset2)
C--   - Manage parameters after creation but before copying a pdf set
C--   - Called by PDFCPY and EVPCOPY
C--
C--   sparBufBase(iadd,ierr)
C--   - Slot 1 <--> LIFO buffer
C--   - Called by PUSHCP and PULLCP
C--
C--   subroutine sparParTo5(islot)
C--   - Copy parameters and pointers of islot to /point5/
C--   - Called by all user routines that need access to a pdf set
C--
C--   sparListPar(islot,ARRAY,IERR)
C--   - Return list of parameters in an array
C--   - Called by CPYPAR
C--
C--   ------------------------------------------------------------------
C--   Package routines
C--   ------------------------------------------------------------------
C--   sparCountUp(key)                     - Increase key counters
C--   sparCountDn(key)                     - Decrease key counters
C--   ------------------------------------------------------------------
C--   sparParAtoB(w1,kset1,w2,kset2)       - Copy parameter list
C--   sparAlfAtoB(w1,kset1,w2,kset2)       - Copy alfas tables
C--   sparPntAtoB(w1,kset1,w2,kset2)       - Copy pointer tables
C--   sparPar6toA(w,kset)                  - Import pars from /qpars6/
C--   sparAtoPar6(w,kset)                  - Export pars to   /qpars6/
C--   sparBaseToKey(islot)                 - Copy slot 1 to islot
C--   sparKeyToBase(islot)                 - Copy islot to slot 1
C--   ------------------------------------------------------------------
C--   iparMakeGroupKey(ipar,np,icntr)      - Make group key
C--   sparMakeBaseKeys                     - Assign and store base keys
C--   iparGetGroupKey(w,kset,igrp)         - Get group key
C--   ------------------------------------------------------------------
C--   sparSetPar(w,kset,ipar,val)          - Set parameter value
C--   dparGetPar(w,kset,ipar)              - Get parameter value
C--   ------------------------------------------------------------------

C--   Parameter management by the par package in QCDNUM
C--   -------------------------------------------------
C--
C--   User input routines write directly in the common block /qpars6/.
C--   The package generates from /qpars6/ the base tables used by the
C--   evolution routines. The base is stored in slot 1 of the pars8
C--   repository and is continuously kept up-to-date. When a pdf set is
C--   filled a snapshot of the base tables is stored in a slot of pars8.
C--   This slot is called the key of the pdf set. Pars8 has mpl0 slots:
C--
C--          1      2                       mset0    mpl0
C--          |      |                         |        |
C--          +------+-------------------------+--------+
C--          | base |        snapshots        |  LIFO  |
C--          +------+-------------------------+--------+
C--
C--   Parameter keys
C--   --------------
C--   All snapshots in pars8 are unique: identical snapshots are never
C--   stored twice. Thus the slot number (key) of the snapshot is a
C--   unique identifier: (un)equal keys mean (un)equal parameters.
C--   Each pdf set in memory has the key of the snapshot stored with
C--   it (pdf_key). For empty pdf sets the pdf_key is zero.
C--
C--   The base set has also an associated key which is the slot where
C--   the next snapshot will be stored (base_key).
C--
C--   A pdf_key always points to a snapshot but the base_key may point
C--   to an empty slot.
C--
C--   Also maintained are key counters with the number of pdf sets
C--   in memory that share that key. A free slot has keycount = 0.
C--
C-    The management of the keycounters and the pars8 repository is
C--   distributed over several routines as follows:
C--
C--   A. Assign a key to the base set (base_key) after a user input:
C--      1. If the base parameter list is found in a slot (>= 2) in
C--         pars8 then set the base_key to the number of that slot.
C--      2. Otherwise set the base_key to the number of the first
C--         free slot (>=2) in pars8.
C--
C--   B. Evolfg, Extpdf and Evdglap: These routines fill a pdf set using
C--      the base parameters and tables as identified by the base_key.
C--      1. If pdf_key = base_key do nothing (except evolve of course).
C--      2. If not then
C--          a. Decrease the key counter[pdf_key] (if pdf_key not 0).
C--          b. Snapshot: copy the base parameters to slot [base_key].
C--          c. Set the pdf_key to the base_key.
C--          d. Increase the key counter[pdf_key].
C--
C--   C. Copy routines: here the parameters come from the source set
C--      instead of the base set so that there is no snapshot taken.
C--      1. If pdf_key = source_key do nothing (except copy the pdf).
C--      2. If not then
C--          a. Decrease the key counter[pdf_key] (if pdf_key not 0).
C--          b. Set the pdf_key to the source_key.
C--          c. Increase the key counter[pdf_key].
C--
C--   The whole mechanism flies if the key counters are initialised to
C--   zero and the pdf_keys are initialized to zero for new pdf sets.
C--
C--   Parameter keys are very handy to quickly decide if parameter lists
C--   are equal or not. There also exist keys for parameter sub-groups.

C==   ==================================================================
C==   Package interface routines  ======================================
C==   ==================================================================

C     ==========================
      subroutine sparInit(nused)
C     ==========================

C--   Initialise parameter repository and keycounters
C--   To be called just after x-q grid is defined

C--   nused  (out):  # words used (< 0 not enough space)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'pstor8.inc'

      dimension itypes(mtyp0)

C--   Initialise
      call smb_VFill(pars8,nwp0,0.D0)
      call smb_IFill(itypes,mtyp0,0)

C--   Go ...
      new       = 0
      npar      = mpar0                                 !parameter store
      nusr      = 0                                          !user store
      itypes(6) = 2*mord0+1                                !alfas tables
      itypes(7) = nsubt0                                 !pointer tables

      ksetrem = 0
      do i = 1,mpl0
        call sqcMakeTab(pars8,nwp0,itypes,npar,nusr,new,kset,nused)
        if(kset.eq.-1) then
          stop 'sparInit: try to create pars8 with no tables'
        elseif(kset.eq.-2) then
          write(lunerr1,'(''STOP sparInit: not enough space'')')
          write(lunerr1,'(''     nwp0 = '',I10)')  nwp0
          write(lunerr1,'('' required = '',I10)') -nused
          write(lunerr1,'(''last slot = '',I10)')  ksetrem
          write(lunerr1,'('' max slot = '',I10)')  mpl0
          stop
        elseif(kset.eq.-3) then
          write(lunerr1,'(''STOP sparInit: max set exceeded'')')
          write(lunerr1,'(''last slot = '',I10)')  ksetrem
          write(lunerr1,'('' max slot = '',I10)')  mpl0
          write(lunerr1,'('' max  set = '',I10)')  mst0
          stop
        elseif(kset.ne.i)  then
          write(lunerr1,'(''STOP sparInit: problem with kset'')')
          write(lunerr1,'(''this slot = '',I10)')  i
          write(lunerr1,'(''kset slot = '',I10)')  kset
          stop
        endif
        ksetrem    = kset
        keyadr8(i) = iqcFirstWordOfParams(pars8,i)-1   !address in pars8
        do j = 1,mgrp0
          keycnt8(i,j) = 0                               !group counters
        enddo
      enddo

C--   Set parameter bits
      ipbits8 = 0
      call smb_sbit1(ipbits8,infbit8)                          !no nfmap
      call smb_sbit1(ipbits8,iasbit8)                   !no alfas tables
      call smb_sbit1(ipbits8,izcbit8)                        !no iz cuts
      call smb_sbit1(ipbits8,ipsbit8)                    !no param store

      return
      end

C     =======================
      subroutine sparMakeBase
C     =======================

C--   Called by every parameter input routine

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

C--   Store not yet initialised, will update later
      if(.not.Lygrid2 .or. .not.Ltgrid2) return
C--   Nothing to do
*mb      write(6,*) 'MAKEBASE ipbits8 =', ipbits8
      if(ipbits8.eq.0)   return
      write(6,*) ' '
*mb      write(6,*) 'MakeBase --------------'
C--   Setup flavour map and pointer tables
      if(imb_gbitn(ipbits8,infbit8).eq.1) then
*mb        write(6,*) 'Hurray call NfTab'
        call sqcNfTab(pars8,1,0)
        call smb_sbit0(ipbits8,infbit8)
      endif
C--   Fill alfas table
      if(imb_gbitn(ipbits8,iasbit8).eq.1) then
*mb        write(6,*) 'Hurray call AlfTab'
        call sqcAlfTab(pars8,1,iord6)
        call smb_sbit0(ipbits8,iasbit8)
      endif
C--   Set iz cuts
      if(imb_gbitn(ipbits8,izcbit8).eq.1) then
*mb        write(6,*) 'Hurray call IzCut'
        call sqcSetIzCut
        call smb_sbit0(ipbits8,izcbit8)
      endif
C--   Fill base parameters
      if(imb_gbitn(ipbits8,ipsbit8).eq.1) then
*mb        write(6,*) 'Hurray store Params'
        call sparPar6toA( pars8, 1 )
        call smb_sbit0(ipbits6,ipsbit6)
      endif
C--   Get parameter group keys
      call sparMakeBaseKeys
C--   Put keys in /qpars6/
      ipver6 = iparGetGroupKey(pars8,1,6)
      itver6 = iparGetGroupKey(pars8,1,3)
*mb      write(6,*) 'MAKEBASE nfix,key =', nfix6,ipver6,itver6

      return
      end

C     ================================
      subroutine sparRemakeBase(islot)
C     ================================

C--   Copy parameters from islot to /qpars6/ and remake the base

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

C--   Check
      if(islot.eq.1) then
        return                         !Base slot so nothing to do
      elseif(islot.le.mset0) then
        if(keycnt8(islot,6).eq.0) stop 'sparRemakeBase: empty slot'
      else
        stop                           'sparRemakeBase: wrong slot'
      endif

C--   Set /qpars6/ parameters
      call sparAtoPar6(pars8,islot)
C--   Update everything
      ipbits8 = 0
      call smb_sbit1(ipbits8,infbit8)                          !no nfmap
      call smb_sbit1(ipbits8,iasbit8)                   !no alfas tables
      call smb_sbit1(ipbits8,izcbit8)                        !no iz cuts
      call smb_sbit1(ipbits8,ipsbit8)                    !no param store
C--   Restore base
      call sparMakeBase

      return
      end

C     =================================
      subroutine sparBufBase(iadd,ierr)
C     =================================

C--   Slot 1 <--> LIFO buffer  (slot mset0+1,....)
C--
C--   iadd  (in) : +1  slot 1 -->  buffer
C--                -1  buffer -->  slot 1
C--   ierr (out) : -1  buffer empty
C--                 0  OK
C--                 1  buffer full
C--                 2  wrong iadd --> do nothing

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

      save ipos
      data ipos /0/

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(ipos)

      if(iadd.eq.+1) then
        ipos = ipos+1
        if(ipos.le.npbuf0) then
          call sparBaseToKey( mset0+ipos )
          ierr   =  0
        else
          ierr = 1
        endif
      elseif(iadd.eq.-1) then
        if(ipos.ge.1) then
          call sparKeyToBase( mset0+ipos )
          ierr   =  0
          ipos   =  ipos-1
        else
          ierr   = -1
        endif
      else
        ierr = 2
      endif

      return
      end

C--   ============================
      subroutine sparParTo5(islot)
C--   ============================

C--   Fill the entire common block /point5/

      implicit double precision(a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'pstor8.inc'

C--   Check
      if(islot.eq.1) then
        ipver = int(dparGetPar(pars8,islot,idipver8))
      elseif(islot.le.mset0) then
        if(keycnt8(islot,6).eq.0) stop 'sparParTo5: empty slot'
        ipver = int(dparGetPar(pars8,islot,idipver8))
        if(ipver.ne.islot)  stop 'sparParTo5: problem with key'
      else
        stop 'sparParTo5: non-existing slot'
      endif

C--   Nothing to do (never for slot 1 = current params)
*mb      write(6,*) 'sparparto5 start'
      if(islot.ne.1 .and. ipver.eq.ipver5) return

*mb      write(6,*)  'sparParTo5 update par5 to version',ipver
      ia8     =  keyadr8(islot)
C--   Update non-threshold parameters
      iord5   = int(pars8(ia8+idiord8))
      nfix5   = int(pars8(ia8+idnfix8))
      aar5    =     pars8(ia8+idaar8)
      bbr5    =     pars8(ia8+idbbr8)
      iymac5  = int(pars8(ia8+idiymac8))
      itmic5  = int(pars8(ia8+iditmic8))
      itmac5  = int(pars8(ia8+iditmac8))
      izmic5  = int(pars8(ia8+idizmic8))
      izmac5  = int(pars8(ia8+idizmac8))
      itmin5  = int(pars8(ia8+iditmin8))
      nzz5    = int(pars8(ia8+idnzz8))
      ipver5  = ipver
C--   Conversions
      ixmic5  = nyy2(0)+1-iymac5
      ymac5   = ygrid2(iymac5)
      xmic5   = exp(-ymac5)
      tmic5   = tgrid2(itmic5)
      tmac5   = tgrid2(itmac5)
      qmic5   = exp(tmic5)
      qmac5   = exp(tmac5)
      ipver5  = ipver

C--   Threshold key
      itver   = int(pars8(ia8+iditver8))
C--   Nothing more to do (never for slot 1 = current params)
      if(islot.ne.1 .and. itver.eq.itver5) return

*mb      write(6,*) 'sparParTo5 update thr5 to version',itver

C--   Update threshold parameters
      id7     = 1000*islot+700
      iatfz   = iqcG7ij(pars8,0,id7+1)
      iazft   = iqcG7ij(pars8,0,id7+2)
      iz15(1) = int(pars8(ia8+idnfmin8))
      iz25(1) = int(pars8(ia8+idnfmax8))
      do nf = 3,6
        iz15(nf) = int(pars8(ia8+idz138-3+nf))
        iz25(nf) = int(pars8(ia8+idz238-3+nf))
      enddo
      do iz = 1,nzz5
        itfiz5( iz) = int(pars8(iatfz+iz))
        itfiz5(-iz) = int(pars8(iatfz-iz))
        izfit5( iz) = int(pars8(iazft+iz))
        izfit5(-iz) = int(pars8(iazft-iz))
      enddo
      itver5  = itver

      return
      end

C     ========================================
      subroutine sparListPar(islot,array,ierr)
C     ========================================

C--   Copy parameter list to a user array

C--   islot  (in) : slot in pars8 repository
C--   array (out) : dimensioned to at least 13 in the calling routine
C--   ierr  (out) : 0 = OK
C--                 1 = islot does not exist
C--                 2 = islot is empty

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'pstor8.inc'

      dimension array(*)

C--   Check
      if(islot.eq.1) then
        continue
      elseif(islot.le.mset0) then
        if(keycnt8(islot,6).eq.0) then
          ierr = 2
          return
        endif
      else
        ierr = 1
        return
      endif

      ierr      = 0
      ia0       = keyadr8(islot)

      array( 1) = pars8(ia0+idiord8)
      array( 2) = pars8(ia0+idalfq08)
      array( 3) = pars8(ia0+idq0alf8)
      array( 4) = pars8(ia0+idnfix8)
      if(array(4).lt.0) then
C--     MFNS
        array(5) = pars8(ia0+idrthrc8)
        array(6) = pars8(ia0+idrthrb8)
        array(7) = pars8(ia0+idrthrt8)
      else
C--     FFNS or VFNS
        array(5) = pars8(ia0+idqthrc8)
        array(6) = pars8(ia0+idqthrb8)
        array(7) = pars8(ia0+idqthrt8)
      endif
      array( 8) = pars8(ia0+idaar8)
      array( 9) = pars8(ia0+idbbr8)
      iymac     = int(pars8(ia0+idiymac8))
      itmic     = int(pars8(ia0+iditmic8))
      itmac     = int(pars8(ia0+iditmac8))
      ymax      = ygrid2(iymac)
      tmin      = tgrid2(itmic)
      tmax      = tgrid2(itmac)
      array(10) = exp(-ymax)
      array(11) = exp( tmin)
      array(12) = exp( tmax)
      array(13) = pars8(ia0+idievtyp8)

      return
      end

C==   ==================================================================
C==   Package routines =================================================
C==   ==================================================================

C     ===========================
      subroutine sparCountUp(key)
C     ===========================

C--   Increase all key counters

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

      if(key.le.1 .or. key.ge.mset0) stop 'sparCountUp: wrong key'

      ia0  = keyadr8(key)+idiover8-1

      do igrp = 1,mgrp0
        keygrp               = int(pars8(ia0+igrp))
        keycnt8(keygrp,igrp) = keycnt8(keygrp,igrp)+1
      enddo

      return
      end

C     ===========================
      subroutine sparCountDn(key)
C     ===========================

C--   Decrease all key counters

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

      if(key.eq.0) return
      if(key.le.1 .or. key.ge.mset0) stop 'sparCountDn: wrong key'

      ia0  = keyadr8(key)+idiover8-1

      do igrp = 1,mgrp0
        keygrp               = int(pars8(ia0+igrp))
        keycnt8(keygrp,igrp) = max(keycnt8(keygrp,igrp)-1,0)
      enddo

      return
      end

C==   ==================================================================

C     =========================================
      subroutine sparParAtoB(w1,kset1,w2,kset2)
C     =========================================

C--   Copy the parameter list

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

      dimension w1(*), w2(*)

      ia1 = iqcFirstWordOfParams(w1,kset1)
      ia2 = iqcFirstWordOfParams(w2,kset2)

      call smb_Vcopy( w1(ia1), w2(ia2), npcopy8 )

      return
      end

C     =========================================
      subroutine sparAlfAtoB(w1,kset1,w2,kset2)
C     =========================================

C--   Copy alfas tables

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w1(*), w2(*)

      do id = -mord0,mord0
        jd1 = 1000*kset1 + 601 + id + mord0
        jd2 = 1000*kset2 + 601 + id + mord0
        call sqcCopyType6(w1,jd1,w2,jd2)
      enddo

      return
      end

C     =========================================
      subroutine sparPntAtoB(w1,kset1,w2,kset2)
C     =========================================

C--   Copy subgrid pointer tables

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w1(*), w2(*)

      do id = 1,nsubt0
        jd1 = 1000*kset1 + 700 + id
        jd2 = 1000*kset2 + 700 + id
        call sqcCopyType7(w1,jd1,w2,jd2)
      enddo

      return
      end

C     ==============================
      subroutine sparPar6toA(w,kset)
C     ==============================

C--   Get params from /qpars6/ and /qgrid2/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'pstor8.inc'

      dimension w(*)

C--   Base address
      ia0 = iqcFirstWordOfParams(w,kset)-1

      w(ia0+idiord8 )  = dble(iord6)
      w(ia0+idalfq08)  = alfq06
      w(ia0+idq0alf8)  = q0alf6
      w(ia0+idnfix8 )  = dble(nfix6)
      w(ia0+idqthrc8)  = qthrs6(4)
      w(ia0+idqthrb8)  = qthrs6(5)
      w(ia0+idqthrt8)  = qthrs6(6)
      w(ia0+idtthrc8)  = tthrs6(4)
      w(ia0+idtthrb8)  = tthrs6(5)
      w(ia0+idtthrt8)  = tthrs6(6)
      w(ia0+idrthrc8)  = rthrs6(4)
      w(ia0+idrthrb8)  = rthrs6(5)
      w(ia0+idrthrt8)  = rthrs6(6)
      w(ia0+idaar8  )  = aar6
      w(ia0+idbbr8  )  = bbr6
      w(ia0+idiymac8)  = dble(iymac2)
      w(ia0+iditmic8)  = dble(itmic2)
      w(ia0+iditmac8)  = dble(itmac2)
      w(ia0+idaslim8)  = aslim6
      w(ia0+idnfmin8)  = dble(nfmin6)
      w(ia0+idnfmax8)  = dble(nfmax6)
      w(ia0+iditmin8)  = dble(itmin6)
      w(ia0+idizmic8)  = dble(izmic2)
      w(ia0+idizmac8)  = dble(izmac2)
      w(ia0+idnzz8)    = dble(nzz2)
      do i = 3,6
        w(ia0+idz138-3+i) = dble(izminf6(i))
        w(ia0+idz238-3+i) = dble(izmanf6(i))
      enddo

      return
      end

C     ==============================
      subroutine sparAtoPar6(w,kset)
C     ==============================

C--   Put params back to /qpars6/ and /qgrid2/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'pstor8.inc'

      dimension w(*)

C--   Base address
      ia0 = iqcFirstWordOfParams(w,kset)-1

      iord6     = int(w(ia0+idiord8))
      alfq06    =     w(ia0+idalfq08)
      q0alf6    =     w(ia0+idq0alf8)
      nfix6     = int(w(ia0+idnfix8))
      qthrs6(4) =     w(ia0+idqthrc8)
      qthrs6(5) =     w(ia0+idqthrb8)
      qthrs6(6) =     w(ia0+idqthrt8)
      tthrs6(4) =     w(ia0+idtthrc8)
      tthrs6(5) =     w(ia0+idtthrb8)
      tthrs6(6) =     w(ia0+idtthrt8)
      rthrs6(4) =     w(ia0+idrthrc8)
      rthrs6(5) =     w(ia0+idrthrb8)
      rthrs6(6) =     w(ia0+idrthrt8)
      aar6      =     w(ia0+idaar8)
      bbr6      =     w(ia0+idbbr8)
      iymac2    = int(w(ia0+idiymac8))
      itmic2    = int(w(ia0+iditmic8))
      itmac2    = int(w(ia0+iditmac8))
      aslim6    =     w(ia0+idaslim8)
      nfmin6    = int(w(ia0+idnfmin8))
      nfmax6    = int(w(ia0+idnfmax8))
      itmin6    = int(w(ia0+iditmin8))
      izmic2    = int(w(ia0+idizmic8))
      izmac2    = int(w(ia0+idizmac8))
      nzz2      = int(w(ia0+idnzz8))
      do i = 3,6
        izminf6(i) = int(w(ia0+idz138-3+i))
        izmanf6(i) = int(w(ia0+idz238-3+i))
      enddo

C--   Restore some more /qpars6/ variables
*      do i = 1,4
*        nfsubg6(i) = 0
*        it1sub6(i) = 0
*        it2sub6(i) = 0
*        iz1sub6(i) = 0
*        iz2sub6(i) = 0
*      enddo
*      ntsubg6 = 0
*      do nf = nfmin6,nfmax6
*        ntsubg6          = ntsubg6 + 1
*        nfsubg6(ntsubg6) = nf
*      enddo
*      it1sub6(1) = 1
*      do i = 2,ntsubg6
*        nf = nfsubg6(i)
*        it1sub6(i) = iqcItfrmT(tthrs6(nf))
*      enddo
*      do i = 1,ntsubg6-1
*        nf = nfsubg6(i)
*        it2sub6(i) = iqcItfrmT(tthrs6(nf+1))
*      enddo
*      it2sub6(ntsubg6) = ntt2
*      do i = 1,ntsubg6
*        nf = nfsubg6(i)
*        iz1sub6(i) = izminf6(nf)
*        iz2sub6(i) = izmanf6(nf)
*      enddo
*      do i = 1,10
*        itlim6(i) = 0
*      enddo
*      itlim6(1) = ntsubg6+1
*      itlim6(2) = 1
*      do i = 1,ntsubg6
*        itlim6(i+1) = it1sub6(i)
*        itlim6(i+6) = nfsubg6(i)
*      enddo
*      itlim6(ntsubg6+2) = ntt2
*      do i = 1,13
*        izlim6(i) = 0
*      enddo
*      izlim6(1) = ntsubg6
*      do i = 1,ntsubg6
*        izlim6(i+1) = iz1sub6(i)
*        izlim6(i+5) = iz2sub6(i)
*        izlim6(i+9) = nfsubg6(i)
*      enddo
*      ipver6 = int(w(ia0+idipver8))
*      itver6 = int(w(ia0+iditver8))

      return
      end

C     ===============================
      subroutine sparBaseToKey(islot)
C     ===============================

C--   Copy pars8 slot 1 to slot islot

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

C--   Check
      if(islot.le.1. .or. islot.gt.mpl0) then
        stop   'sparBaseToKey: wrong slot'
      endif

      call sparParAtoB(pars8,1,pars8,islot)
      call sparAlfAtoB(pars8,1,pars8,islot)
      call sparPntAtoB(pars8,1,pars8,islot)

      return
      end

C     ===============================
      subroutine sparKeyToBase(islot)
C     ===============================

C--   Copy pars8 slot key to slot 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

C--   Check
      if(islot.le.1. .or. islot.gt.mpl0) then
        stop   'sparKeyToBase: wrong slot'
      endif

      call sparParAtoB(pars8,islot,pars8,1)
      call sparAlfAtoB(pars8,islot,pars8,1)
      call sparPntAtoB(pars8,islot,pars8,1)

      return
      end

C==   ==================================================================

C     ================================================
      integer function iparMakeGroupKey(ipar,np,icntr)
C     ================================================

C--   Assign group key

C--   ipar   (in) : id of first parameter in group
C--   np     (in) : number of parameters in group
C--   icntr  (in) : group counter

C--   The routine looks for a slot with the same parameters as slot 1.
C--   If found, the key is set to that slot, if not it is set to the
C--   first empty slot (a slot with keycount 0).

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

      logical lmb_Vcomp

      dimension icntr(*)

C--   Slot not found
      iparMakeGroupKey = 0
C--   Ipar address slot 1
      ia1 = keyadr8(1)+ipar
C--   Loop backward over slots (but not slot 1)
      do i = mset0,2,-1
C--     Ipar address slot i
        iai = keyadr8(i)+ipar
        if(icntr(i).eq.0) then
C--       Empty slot
          iparMakeGroupKey = i
        elseif( lmb_Vcomp(pars8(ia1),pars8(iai),np,0.D0) ) then
C--       Found an identical slot
          iparMakeGroupKey = i
          return
        endif
      enddo

      if(iparMakeGroupKey.eq.0) stop
     +                         'iparMakeGroupKey: parameter store full'

      return
      end

C     ===========================
      subroutine sparMakeBaseKeys
C     ===========================

C--   Make and store all base keys

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

      ia0 = keyadr8(1)

C--   1 = Order
      key = iparMakeGroupKey(iparo8,nparo8,keycnt8(1,1))
      pars8(ia0+idiover8) = dble(key)
C--   2 = Alfas
      key = iparMakeGroupKey(ipara8,npara8,keycnt8(1,2))
      pars8(ia0+idiaver8) = dble(key)
C--   3 = Thresholds
      key = iparMakeGroupKey(ipart8,npart8,keycnt8(1,3))
      pars8(ia0+iditver8) = dble(key)
C--   4 = Scale
      key = iparMakeGroupKey(ipars8,npars8,keycnt8(1,4))
      pars8(ia0+idisver8) = dble(key)
C--   5 = Cuts
      key = iparMakeGroupKey(iparc8,nparc8,keycnt8(1,5))
      pars8(ia0+idicver8) = dble(key)
C--   6 = All parameters
      key = iparMakeGroupKey(iparp8,nparp8,keycnt8(1,6))
      pars8(ia0+idipver8) = dble(key)

      return
      end

C     =============================================
      integer function iparGetGroupKey(w,kset,igrp)
C     =============================================

C--   Get group key

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'pstor8.inc'

      dimension w(*)

      ia0 = iqcFirstWordOfParams(w,kset) + idiover8 - 2
      iparGetGroupKey =  int(w(ia0+igrp))

      return
      end

C==   ==================================================================

C     ======================================
      subroutine sparSetPar(w,kset,ipar,val)
C     ======================================

C--   Set parameter value

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*)

      ia0         = iqcFirstWordOfParams(w,kset)-1
      w(ia0+ipar) = val

      return
      end

C     =================================================
      double precision function dparGetPar(w,kset,ipar)
C     =================================================

C--   Get parameter value

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*)

      ia0        = iqcFirstWordOfParams(w,kset)-1
      dparGetPar = w(ia0+ipar)

      return
      end
