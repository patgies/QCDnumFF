
C-----------------------------------------------------------------------
      program testwspace
C-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      logical sjekhead

      parameter( nwords = 573 )
      parameter( ntags  = 7 )

      dimension w(nwords)
      dimension imi(2),ima(2)
      dimension kk(27),imin(25),imax(25)
      dimension index(2),row(10)
      character*80 fname
      dimension head(20),hstore(20,7)
      data iSet1/2/, iTab1/3/, iTab2/4/, iSet2/5/, iTab3/6/, iTab4/7/

      write(6,'(/A,I10)') ' WSTORE version = ', iws_Version()

      ndim   =  2
      imi(1) =  1
      ima(1) =  5
      imi(2) =  1
      ima(2) =  5

      write(6,*) ' '
      write(6,*) '=============== init workspace ======================'
      iaSet1 = iws_wsinit(w,nwords,ntags,' ')
      write(6,*) ' '
      write(6,*) '=============== add 1st table  ======================'
      iaTab1 = iws_wtable(w,imi,ima,ndim)
      write(6,*) ' '
      write(6,*) '=============== add 2nd table  ======================'
      iaTab2 = iws_wtable(w,imi,ima,ndim)
      write(6,*) ' '
      write(6,*) '=============== add 2nd tbset  ======================'
      iaSet2 = iws_newset(w)
      write(6,*) ' '
      write(6,*) '=============== add 3rd table  ======================'
      iaTab3 = iws_wtable(w,imi,ima,ndim)
      write(6,*) ' '
      write(6,*) '=============== add 4th table  ======================'
      iaTab4 = iws_wtable(w,imi,ima,ndim)
*      write(6,*) ' '
*      write(6,*) '=============== add 3rd tbset  ======================'
*      iaSet3 = iws_newset(w)
*      call sws_WShead(w,iaSet3)
*      stop

C--   Store headers
      call fillhead(w(1),hstore(1,1))
      call fillhead(w(iaSet1),hstore(1,iSet1))
      call fillhead(w(iaTab1),hstore(1,iTab1))
      call fillhead(w(iaTab2),hstore(1,iTab2))
      call fillhead(w(iaSet2),hstore(1,iSet2))
      call fillhead(w(iaTab3),hstore(1,iTab3))
      call fillhead(w(iaTab4),hstore(1,iTab4))

C--   Fill table 1
      index(1) = imi(1)
      index(2) = imi(2)
      ia = iws_TPoint(w,iaTab1,index,2)-1
      do j = imi(2),ima(2)
        do i = imi(1),ima(1)
          ia    = ia + 1
          w(ia) = dble(10*i+j)
        enddo
      enddo
C--   Print headers
      call sws_WShead(w,1)
      call sws_WShead(w,iaSet1)
      call sws_WShead(w,iaTab1)
      call sws_WShead(w,iaTab2)
      call sws_WShead(w,iaSet2)
      call sws_WShead(w,iaTab3)
      call sws_WShead(w,iaTab4)
C--   Print tree
      call sws_WStree(w)
C--   Print table 1
      write(6,*)," "
      write(6,*) 'Contents of table 1'
      do i = imi(1),ima(1)
        index(1) = i
        do j = imi(2),ima(2)
          index(2) = j
          ia = iws_TPoint(w,iaTab1,index,2)
          row(j) = w(ia)
        enddo
        write(6,*) (int(row(j)),j=imi(2),ima(2))
      enddo

C--   Last table and last set address
      Ltab  = iws_TFskip(w,iaTab3)
      Lset  = iws_SBskip(w,iaSet2)
      write(6,*) ' '
      write(6,*) 'ia last  table = ',Ltab + iaTab3
      write(6,*) 'ia first tbset = ',Lset + iaSet2

      write(6,*) ' '
      write(6,*) 'Check WsWipe, TsDump and TsRead'
      write(6,*) '-------------------------------'
C--   Remove last table
      call sws_WsWipe(w,iaTab4)
      write(6,*)  'Wipe  table 4   ',iaTab4
C--   Dump second table-set
      fname = '../weights/crap.wgt'
      idump = iaSet2
      call sws_TsDump(fname,0,w,idump,ierr)
      if(ierr.ne.0) then
        write(6,*) 'Error openening or writing file ',fname,' '
        stop
      else
        write(6,*) 'Dump  tableset 2',idump
      endif
C--   Wipe second table set
      call sws_WsWipe(w,iaSet2)
      write(6,*)  'Wipe  tableset 2',iaSet2
C--   Create new (empty) table set
      iadum = iws_NewSet(w)
      write(6,*)  'New   tableset  ',iadum
C--   Read second table set
      iaSet2 = iws_TsRead(fname,0,w,ierr)
      if(ierr.ne.0) then
        write(6,*) 'Error openening or reading file ',fname
        stop
      else
        write(6,*)  'Read  tableset 2',iaSet2
      endif
C--   Restore last table
      iaTab4 = iws_WTable(w,imi,ima,ndim)
      write(6,*)  'New   table 4   ',iaTab4
C--   Check headers
      write(6,*) ' '
      call fillhead(w(1),head)
      write(6,*) 'workspace header OK = ',sjekhead(head,hstore(1,1))
      call fillhead(w(iaSet2),head)
      write(6,*) 'set 2 header     OK = ',sjekhead(head,hstore(1,iSet2))
      call fillhead(w(iaTab3),head)
      write(6,*) 'tab 3 header     OK = ',sjekhead(head,hstore(1,iTab3))
      call fillhead(w(iaTab4),head)
      write(6,*) 'tab 4 header     OK = ',sjekhead(head,hstore(1,iTab4))

      write(6,*) ' '
      write(6,*) 'Check TbCopy'
      write(6,*) '------------'
      write(6,*) 'Contents of table 4 before copy'
      do i = imi(1),ima(1)
        index(1) = i
        do j = imi(2),ima(2)
          index(2) = j
          ia = iws_Tpoint(w,iaTab4,index,2)
          row(j) = w(ia)
        enddo
        write(6,*) (int(row(j)),j=imi(2),ima(2))
      enddo
      write(6,*) 'Copy table 1 to table 4'
      call sws_TbCopy(w,iaTab1,w,iaTab4,1)
      write(6,*) 'Contents of table 4 after copy'
      do i = imi(1),ima(1)
        index(1) = i
        do j = imi(2),ima(2)
          index(2) = j
          ia = iws_Tpoint(w,iaTab4,index,2)
          row(j) = w(ia)
        enddo
        write(6,*) (int(row(j)),j=imi(2),ima(2))
      enddo

      write(6,*) ' '
      write(6,*) 'Check WClone'
      write(6,*) '------------'
C--   Wipe second table set
      call sws_WsWipe(w,iaSet2)
      write(6,*)  'Wipe  tableset 2',iaSet2
C--   Create new (empty) table set
      iadum = iws_NewSet(w)
      write(6,*)  'New   tableset  ',iadum
C--   Clone table set 1
      iaSet2 = iws_WClone(w,iaSet1,w)
      write(6,*)  'Clone set 1 to 2',iaSet2
C--   Wipe table 3
      call sws_WsWipe(w,iaTab3)
      write(6,*)  'Wipe  table 3+4 ',iaTab3
C--   Clone table 2 to 3
      iaTab3 = iws_WClone(w,iaTab2,w)
      write(6,*)  'Clone tab 2 to 3',iaTab3
C--   Clone table 3 to 4
      iaTab4 = iws_WClone(w,iaTab3,w)
      write(6,*)  'Clone tab 3 to 4',iaTab4
C--   Check headers
      write(6,*) ' '
      call fillhead(w(1),head)
      write(6,*) 'workspace header OK = ',sjekhead(head,hstore(1,1))
      call fillhead(w(iaSet2),head)
      write(6,*) 'set 2 header     OK = ',sjekhead(head,hstore(1,iSet2))
      call fillhead(w(iaTab3),head)
      write(6,*) 'tab 3 header     OK = ',sjekhead(head,hstore(1,iTab3))
      call fillhead(w(iaTab4),head)
      write(6,*) 'tab 4 header     OK = ',sjekhead(head,hstore(1,iTab4))

C--   Print tree
      call sws_WsTree(w)

C--   Metadata field
      write(6,*) ' '
      write(6,*) 'Fingerprint and metadata table1 with address',iaTab1
      call swsGetMeta(w,iaTab1,nd,kk,imin,imax)
      write(6,'( '' ifp   = '', I15 )') kk(1)
      write(6,'( ''  nd   = '', I15 )') nd
      do i = 1,nd
        write(6,'( '' lim   = '', I15,2i10 )') i,imin(i),imax(i)
      enddo

C--   Build-in pointer function
      write(6,*) ' '
      write(6,*) 'Test iws_Tpoint function on table1'
      index(1) =  imi(1)
      index(2) =  imi(2)
      ia = iws_TPoint(w,iatab1,index,2)
      ib = iws_BeginTBody(w,iatab1)
      write(6,*) 'begin table = ',ia, ib
      index(1) =  ima(1)
      index(2) =  ima(2)
      ia = iws_TPoint(w,iatab1,index,2)
      ib = iws_EndTBody(w,iatab1)
      write(6,*) 'end   table = ',ia, ib

      call sws_WsMark(icWsp,icSet,icTab)

      write(6,*) ' '
      nsize = iws_TbSize(imi,ima,ndim)
      write(6,*) 'tsize  = ', nsize
      nsize = iws_HdSize()
      write(6,*) 'hsize  = ', nsize
      write(6,*) 'Wspace = ', icWsp
      write(6,*) 'TbSet  = ', icSet
      write(6,*) 'Table  = ', icTab
      iyesn = iws_IsaWorkspace(q)
      write(6,*) 'Yesno  = ', iyesn
      nword = iws_SizeOfW(w)
      write(6,*) 'Wsize  = ', nword
      nword = iws_WordsUsed(w)
      write(6,*) 'Wused  = ', nword
      nword = iws_Nheader(w)
      write(6,*) 'Nhead  = ', nword
      nword = iws_Ntags(w)
      write(6,*) 'Ntags  = ', nword
      nword = iws_HeadSkip(w,iaSet2)
      write(6,*) 'Hskip  = ', nword
      nword = iws_IaDrain(w,iatab2)
      write(6,*) 'Drain  = ', nword
      nword = iws_IaNull(w,iaSet1)
      write(6,*) 'Inull  = ', nword
      ityp0 = iws_ObjectType(w,66)
      ityp1 = iws_ObjectType(w,1)
      ityp2 = iws_ObjectType(w,iaSet1)
      ityp3 = iws_ObjectType(w,iaTab2)
      write(6,*) 'Itype  = ', ityp0,ityp1,ityp2,ityp3
      nword = iws_ObjectSize(w,1)
      write(6,*) 'WSsize = ', nword
      nword = iws_ObjectSize(w,iaSet1)
      write(6,*) 'TSsize = ', nword
      nword = iws_ObjectSize(w,iaTab1)
      write(6,*) 'TBsize = ', nword
      nobj  = iws_Nobjects(w,1)
      write(6,*) 'Nsets  = ', nobj
      nobj  = iws_Nobjects(w,iaset1)
      iobj  = iws_ObjectNumber(w,iaSet1)
      write(6,*) 'Ntabs  = ', iobj, nobj
      nobj  = iws_Nobjects(w,iaset2)
      iobj  = iws_ObjectNumber(w,iaSet2)
      write(6,*) 'Ntabs  = ', iobj, nobj
      ifp   = iws_FingerPrint(w,1)
      write(6,*) 'Wfprnt = ', ifp
      ifp   = iws_FingerPrint(w,iaSet1)
      write(6,*) 'S1fp   = ', ifp
      ifp   = iws_FingerPrint(w,iaSet2)
      write(6,*) 'S2fp   = ', ifp
      ifp   = iws_FingerPrint(w,iaTab1)
      write(6,*) 'T1fp   = ', ifp
      ifp   = iws_FingerPrint(w,iaTab2)
      write(6,*) 'T2fp   = ', ifp
      itag  = iws_IaFirstTag(w,1)
      write(6,*) 'Wtag   = ',itag
      mdim  = iws_TableDim(w,iaTab2)
      write(6,*) 'Ndim   = ',mdim
      iakk  = iws_IaKARRAY(w,iaTab1)
      write(6,*) 'KARR1  = ',iakk
      iadr  = iws_IaIMIN(w,iaTab1)
      write(6,*) 'IMIN1  = ',iadr
      iadr  = iws_IaIMAX(w,iaTab1)
      write(6,*) 'IMAX1  = ',iadr
      iadr = iws_BeginTbody(w,iaTab1)
      write(6,*) 'Body1  = ',iadr
      iadr = iws_EndTbody(w,iaTab1)
      write(6,*) 'Bend1  = ',iadr

      end

C     ===========================
      subroutine fillhead(w,head)
C     ===========================

      implicit double precision (a-h,o-z)

      dimension w(*), head(*)

      n = iws_HdSize()
      call smb_Vcopy(w,head,n)

      return
      end

C     ======================================
      logical function sjekhead(head1,head2)
C     ======================================

      implicit double precision (a-h,o-z)
      logical lmb_Vcomp

      dimension head1(*), head2(*)

      n        = iws_HdSize()
      sjekhead = lmb_Vcomp(head1,head2,n,-1.D-9)

      return
      end
