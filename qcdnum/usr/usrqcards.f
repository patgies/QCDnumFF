
C--   This is the file usrqcards.f with datacard routines

C--   subroutine qcards( usub, filename, iprint )
C--   subroutine qcbook(action,key)

C     ===========================================
      subroutine qcards( usub, filename, iprint )
C     ===========================================

C--   Read datacard file
C--
C--   usub     (in) : subroutine to be called when user key encountered
C--   filename (in) : name of the datacard file
C--   iprint   (in) : #0 = print, 0 = run quiet
C--
C--   NB: The datacard file will always be closed after reading

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      external usub
      character*(*) filename

      character*7   key7
      character*37  emsg(11)
C--                 1234567890123456789012345678901234567
      data emsg   /'       : error reading parameter list',     ! 1
     +             '       : problem processing keycard  ',     ! 2
     +             '       : dont know how to process    ',     ! 3
     +             'Error reading datacard file          ',     ! 4
     +             'Cannot open datacard file            ',     ! 5
     +             '       : empty parameter list        ',     ! 6
     +             '       : unbalanced quotes           ',     ! 7
     +             '       : too many words (>100)       ',     ! 8
     +             '       : wordlength exceeded (>120)  ',     ! 9
     +             '       : parameter string too large  ',     !10
     +             '       : format descriptor too large '   /  !11

      character*80 subnam
      data subnam /'QCARDS ( MYCARDS, FILENAME, IPRINT )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Free logical unit number
      lun = iqcLunFree(10)

C--   Read and process the datacard file (does not close lun)
      write(lunerr1,'(/'' QCARDS: read datacards from '',A)') filename
      call sqcQcards( usub, lun, filename, iprint, ierr, key7 )
      close(lun)

      if(ierr.eq.0) then
        return
      elseif(ierr.ne.4 .and. ierr.ne.5) then
        emsg(ierr)(1:7) = key7
      endif

C--   Error
      call sqcErrMsg(subnam,emsg(ierr))

      return
      end

C     =============================
      subroutine qcbook(action,key)
C     =============================

C--   Add, delete, list keys in common block qcard9
C--
C--   action   (in) :  'A...' = add, 'D...' = delete, 'L...' = list
C--   key      (in) :  key to add or delete
C--
C--   NB: Only admitted are keys with a length of 1-7 characters, with
C--       no embedded blanks.
C--   NB: Added keys are always of type USER.
C--   NB: One cannot add a key with a name that already exists.
C--   NB: Key names are always stored in upper case.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qcard9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*(*) action, key
      character*1   opt
      character*20  message
      data message /'    : Unknown action'/
      character*34  emsg(5)
C--                 1234567890123456789012345678901234
      data emsg   /'Input key is empty                ',
     +             'Input key longer than 7 characters',
     +             'Input key has embedded blanks     ',
     +             'Input key already exists          ',
     +             'List of known keys is full        '   /

      character*80 subnam
      data subnam /'QCBOOK ( ACTION, KEY )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Use mbutil for character string manipulations
      opt   = ' '
      i1    = imb_frstc(action)
      if(i1.ne.0) then
        opt = action(i1:i1)
        call smb_cltou(opt)
      endif

C--   Switch through options
      ierr = 0
      if(opt.eq.'L') then
        write(lunerr1,'(/''  List of predefined and user keys ''/
     +                   ''  -------------------------------- '')')
        j = 0
        do i = 1,mky0
          if(qkeys9(i)(9:12) .ne. 'FREE') then
            j = j+1
            write(lunerr1,'(I4,2X,A)') j, qkeys9(i)
          endif
        enddo
      elseif(opt.eq.'A' .or. opt.eq.'D') then
        call sqcQcBook( opt, key, ierr )
      elseif(i1.eq.0) then
        call sqcErrMsg(subnam,'Empty action string')
      else
        message(2:2) = opt
        call sqcErrMsg(subnam,message)
      endif

C--   Errors from sqcQcBook
      if(ierr.ne.0) then
        call sqcErrMsg(subnam,emsg(ierr))
      endif

      return
      end



