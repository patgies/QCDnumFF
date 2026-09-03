
C--   This is the file srcQcards.f with datacard routines

C--   subroutine sqcQcards( usub, lun, filename, iprint, ierr, key7 )
C--   subroutine sqcQcProc( key7, nk, pars, np, fmt, nf, ierr )
C--   subroutine sqcQcsplit( dcard, key7, nk, rest, nr )
C--   subroutine sqcQcBook( action, key, ierr )
C--   subroutine sqcChecKey( keyin, key7, nk, ierr )
C--   integer function iqcFindKey( key, kind )

C     ===============================================================
      subroutine sqcQcards( usub, lun, filename, iprint, ierr, key7 )
C     ===============================================================

C--   usub     (in) : subroutine called when USER card encountered
C--   lun      (in) : logical unit number  (will remain open at exit)
C--   filename (in) : input datacard file
C--   iprint   (in) : 0 not print; >0 print; <0 only print (dry run)
C--   ierr    (out) : 0 = OK
C--   key7    (out) : last key read

C--   ierr  (out)  0 = OK
C--                1 = error reading parameter list
C--                2 = problem processing keycard
C--                3 = unknown key
C--                4 = error reading input file
C--                5 = cannot open input file
C--                6 = empty parameter list
C--                7 = found unbalanced quotes
C--                8 = wordcount  exceeded
C--                9 = wordlength exceeded
C--               10 = not enough space in pars
C--               11 = not enough space in fmt

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      character*(*) filename
      character*120 line, rest
      character*120 pars
      character*200 fmt
      character*7   key7
      character*4   kind
      logical eof

      external usub

      jprint = iprint !avoid compiler warning
      eof    = .false.

      open(unit=lun,file=filename,status='old',err=500)
      rewind lun

      if(iprint.ne.0) write(lunerr1,'('' '')')

      do while(.not.eof)

        read(lun,*,end=300,err=400) line
        call sqcQcsplit( line, key7, nk, rest, nr )
        if(key7.eq.'QCSTOP') goto 300

        ifound = iqcFindKey( key7, kind )

        if(ifound .ne. 0) then

          call smb_sfmat( rest, pars, fmt, jerr )
          if(jerr.ge.2) then
            ierr = jerr+5
            return
          endif
          np = imb_lenoc(pars)
          nf = imb_lenoc(fmt)

          if(iprint.ne.0 .and. np.eq.0)     then
            write(lunerr1,'(A8)')           key7
          elseif(abs(iprint).eq.1) then
            write(lunerr1,'(A8,2X,A)')      key7, pars(1:np)
          elseif(abs(iprint).ge.2) then
            write(lunerr1,'(A8,2X,A,2X,A)') key7, pars(1:np), fmt(1:nf)
          endif

          if(iprint.ge.0) then
            if    (kind .eq. 'QKEY') then
              if(np.eq.0) then             !all cards must have params
                ierr = 6
                return
              endif
              call sqcQcProc( key7, nk, pars, np, fmt, nf, ierr )
              if(ierr .ne. 0) return
            elseif(kind .eq. 'USER') then
              call      usub( key7, nk, pars, np, fmt, nf, ierr )
              if(ierr .ne. 0) return
            else
              stop 'sqcQcards: unknown type of key'
            endif
          endif

        endif

      enddo

 300  continue
      ierr = 0
      return
 400  continue
      ierr = 4
      return
 500  continue
      ierr = 5
      return

      end


C     =========================================================
      subroutine sqcQcProc( key7, nk, pars, np, fmt, nf, ierr )
C     =========================================================

C--   Process predefined keycards
C--
C--   key7   (in)  keyword 7 chars max in upper case
C--   nk     (in)  number of characters in key
C--   pars   (in)  character string with arguments
C--   np     (in)  number of characters in pars 0=no argument list
C--   ierr  (out)  0 = OK
C--                1 = error reading parameter list
C--                2 = problem processing keycard
C--                3 = unknown key

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      character*7    key7
      character*4    key4
      character*(*)  pars, fmt
      character*120  fnam

      dimension xmi(mxg0), iwt(mxg0), qqi(20), qwt(20)

      data iwt /1, 2, 4, 8, 16/
      data qwt /20*1.D0/

C--   Nothing to process (change this if there are predefined keys
C--   without parameters to process)
      if(nk.eq.0 .or. np.eq.0 .or. nf.eq.0) then
        ierr = 0
        return
      endif

      if(key7(1:6).eq.'SETLUN') then
        read(unit=pars,fmt=fmt,err=10,end=10)   lun,fnam
        call SETLUN( lun, fnam )
        return
   10   read(unit=pars,fmt=fmt,err=100,end=100) lun
        call SETLUN( lun, ' ' )

      elseif(key7(1:6).eq.'SETVAL') then
        read(unit=pars,fmt=fmt,err=100,end=100) key4,val
        call SETVAL( key4, val )

      elseif(key7(1:6).eq.'SETINT') then
        read(unit=pars,fmt=fmt,err=100,end=100) key4,ival
        call SETINT( key4, ival )

      elseif(key7(1:6).eq.'GXMAKE') then
        read(unit=pars,fmt=fmt,err=100,end=100)
     +       ioy,nx,n,(xmi(i),i=1,min(n,mxg0))
        call GXMAKE( xmi, iwt, min(n,mxg0), nx, nxout, ioy )

      elseif(key7(1:6).eq.'GQMAKE') then
        read(unit=pars,fmt=fmt,err=100,end=100)
     +       nq,n,(qqi(i),i=1,min(n,20))
        call GQMAKE( qqi, qwt, min(n,20), nq, nqout )

      elseif(key7(1:6).eq.'FILLWT') then
        read(unit=pars,fmt=fmt,err=20,end=20) itype,fnam
        j1  = imb_frstc(fnam)
        j2  = imb_lenoc(fnam)
        lun = iqcLunFree(0)
        call READWT( lun, fnam(j1:j2), idmin, idmax, nwds, jerr )
        if(jerr.ne.0) then
          call FILLWT( itype, idmin, idmax, nwds )
          call DMPWGT( itype, lun, fnam(j1:j2) )
        endif
        return
   20   continue
        read(unit=pars,fmt=fmt,err=100,end=100) itype
        call FILLWT( itype, idmin, idmax, nwds )

      elseif(key7(1:6).eq.'SETORD') then
        read(unit=pars,fmt=fmt,err=100,end=100) iord
        call SETORD( iord )

      elseif(key7(1:6).eq.'SETALF') then
        read(unit=pars,fmt=fmt,err=100,end=100) as0,r20
        call SETALF( as0, r20 )

      elseif(key7(1:6).eq.'SETCBT') then
        read(unit=pars,fmt=fmt,err=100,end=100) nfix,q2c,q2b,q2t
        call SETCBT( nfix, iqfrmq(q2c), iqfrmq(q2b), iqfrmq(q2t) )

      elseif(key7(1:6).eq.'MIXFNS') then
        read(unit=pars,fmt=fmt,err=100,end=100) nfix,r2c,r2b,r2t
        call MIXFNS( nfix, r2c, r2b, r2t )

      elseif(key7(1:6).eq.'SETABR') then
        read(unit=pars,fmt=fmt,err=100,end=100) ar,br
        call SETABR( ar, br )

      elseif(key7(1:6).eq.'SETCUT') then
        read(unit=pars,fmt=fmt,err=30,end=30) xmi,q2mi,q2ma,ddum
        call SETCUT( xmi, q2mi, q2ma, ddum )
        return
   30   read(unit=pars,fmt=fmt,err=100,end=100) xmi,q2mi,q2ma
        call SETCUT( xmi, q2mi, q2ma, 0.D0 )


      else
        ierr = 3
        return
      endif

      ierr = 0
      return

 100  ierr = 1
      return
      end

C     ===================================================
      subroutine sqcQcsplit( dcard, key7, nk, rest, nr )
C     ===================================================

C--   Split a datacard into a keyword and a parameter list
C--   A keyword is the first word of dcard provided 1 <= #chars <= 7
C--   If there is no keyword then there is also no rest
C--
C--   dcard  (in) : datacard
C--   key7  (out) : key in upper case (always .le. 7 characters)
C--   nk    (out) : length of key  ; 0 = no  key found or dcard empty
C--   rest (out)  : parameter list
C--   nr    (out) : length of rest ; 0 = no rest found or dcard empty

      implicit double precision (a-h,o-z)

      character*(*) dcard
      character*(*) rest
      character*(*) key7

      call smb_cfill(' ',rest)
      call smb_cfill(' ',key7)
      nk = 0
      nr = 0

      i1c    = imb_frstc(dcard)
      i2c    = imb_lenoc(dcard)

      if(i2c.eq.0) then
C--     Empty dcard
        return
      else
C--     Find limits of first word
        i1w = i1c
        i2w = i1w
        do i = i1c,i2c
          i2w = i
          if(dcard(i:i).eq.' ') then
            i2w = i-1
            goto 10
          endif
        enddo
  10    continue
      endif
C--   Check first word is valid keyword
      lenw = i2w-i1w+1
      if(lenw.gt.7) then
C--     Invalid key
        return
      else
C--     Valid key
        key7 = dcard(i1w:i2w)
        call smb_cltou(key7)
        nk   = lenw
      endif
C--   Put remnant as rest
      i1p = i2w+1
      if(i1p.gt.i2c) then
C--     No parameter leist
        return
      else
C--     Parameter list found
        rest = dcard(i1p:i2c)
        nr   = imb_lenoc(rest)
      endif

      return
      end

C     =========================================
      subroutine sqcQcBook( action, key, ierr )
C     =========================================

C--   Add or delete a key in the qkeys9 list of keys

C--   action      (in)  'A' add or 'D' delete
C--   key         (in)  key to be added or deleted
C--   ierr       (out)  0 = OK
C--                     1 = input key empty
C--                     2 = input key longer than 7 characters
C--                     3 = input key has embedded blanks
C--                     4 = input key already exists
C--                     5 = list is full

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcard9.inc'

      character*(*) action, key
      character*4   kind
      character*7   key7

      if    (action(1:1) .eq. 'D') then
        call sqcChecKey(key,key7,nk,ierr)
        if(ierr.ne.0) return
        i = iqcFindKey(key7,kind)
        if(i.ne.0) then
          qkeys9(i)(1: 8) = '        '
          qkeys9(i)(9:12) = 'FREE'
        endif
        ierr = 0
      elseif(action(1:1) .eq. 'A') then
        call sqcChecKey(key,key7,nk,ierr)
        if(ierr.ne.0) then
          return
        elseif(iqcFindKey(key7,kind) .ne. 0) then
          ierr = 4
          return
        else
          ierr = 5
          do i = 1,mky0
            if(qkeys9(i)(9:12) .eq. 'FREE') then
              qkeys9(i)(1: 7) =  key7
              qkeys9(i)(9:12) = 'USER'
              ierr            = 0
              return
            endif
          enddo
        endif
      else
        stop 'sqcQcBook: unknown action'
      endif

      return
      end

C     ==============================================
      subroutine sqcChecKey( keyin, key7, nk, ierr )
C     ==============================================

C--   Strip leading and trailing blanks, transform to upper case
C--   and return the key in key7
C--
C--   keyin   (in)  input character string
C--   key7   (out)  output key
C--   nk     (out)  number of characters .le. 7 ; 0=empty
C--   ierr   (out)  0 = OK
C--                 1 = empty key
C--                 2 = input key longer than 7 characters
C--                 3 = input key has embedded blanks

      implicit double precision (a-h,o-z)

      character*(*) keyin
      character*7   key7

C--           1234567
      key7 = '       '
      nk   = 0
      i1   = imb_frstc(keyin)
      i2   = imb_lenoc(keyin)
      if(i2.eq.0) then
        ierr = 1
        return
      elseif(i2-i1+1 .gt. 7) then
        ierr = 2
        return
      else
        do i = i1,i2
          if(keyin(i:i) .eq. ' ') then
            ierr = 3
          endif
        enddo
      endif
      key7 = keyin(i1:i2)
      call smb_cltou(key7)

      return
      end

C     ========================================
      integer function iqcFindKey( key, kind )
C     ========================================

C--   Find a key in the qkeys9 list of keys
C--
C--   key         (in)   key to search for (in upper case)
C--   kind       (out)   kind of key QKEY or USER
C--   iqcFindKey (out)   position in the list;  0 = not found

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcard9.inc'

      character*(*) key
      character*4   kind

      iqcFindKey = 0
      kind       = '    '
      if(imb_lenoc(key).eq.0) return
      do i = 1,mky0
        if( key .eq. qkeys9(i)(1:7) )  then
          iqcFindKey = i
          kind       = qkeys9(i)(9:12)
          return
        endif
      enddo

      return
      end



