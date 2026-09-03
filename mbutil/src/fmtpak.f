
C--   This is the file fmtpak.f with string formatting routines

C--   subroutine smb_sfmat(stin, stout, fmt, ierr)
C--   subroutine sfmtPutwd(separator, word, cstring, ierr)
C--   subroutine sfmtRefmt(field, strin, strout, nw, nd)
C--   subroutine sfmtSform(field, nw, nd, cformat, nf)
C--   subroutine sfmtStype(cstring, ctype, field, nw, nd)
C--   integer function ifmtDEFIC(cstring, cformat)
C--   logical function sfmtDorE(cstring, i, idot)
C--   logical function sfmtReal(cstring, i)
C--   logical function sfmtInte(cstring, n)
C--   logical function sfmtUInt(cstring, n)
C--   subroutine sfmtParseit(idim, cstring, iw1, iw2, nw, ierr)
C--   subroutine sfmtGetWord(cstring, i1, k1, k2, ierr)
C--   integer function ifmtFstChar(cstring, i)
C--   integer function ifmtFstEofW(cstring, i1)

C     ============================================
      subroutine smb_sfmat(stin, stout, fmt, ierr)
C     ============================================

C--   Reformat an input string and return fortran format descriptor
C--
C--   stin    (in) : input string
C--   stout  (out) : reformatted string
C--   fmt    (out) : format descriptor
C--   ierr   (out) : 1 = empty input string
C--                  2 = unbalanced quotes in stin
C--                  3 = more than mxwd words in stin
C--                  4 = wordlength lenw exceeded or zero length word
C--                  5 = not enough room in stout
C--                  6 = not enough room in fmt
C--
C--   Limits hardwired and checked: mxwd=100 and lenw=120
C--   Limit  not checked          : lenf=30 (should be more than enough)
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) stin, stout, fmt
      parameter( mxwd = 100 )
      dimension iw1(mxwd),iw2(mxwd)
      parameter( lenw = 120 )
      character*120 wordout
      parameter( lenf = 30 )
      character*30  cformat
      character*4   ctype
      character*1   field

C--   Initialize
      ierr = 0
      call smb_cfill(' ',stout)
      call smb_cfill(' ',fmt)
C--   Check empty string
      nn = imb_lenoc(stin)
      if(nn.eq.0) then
        ierr = 1
        return
      endif
C--   Parse input string
      call sfmtParseit(mxwd,stin,iw1,iw2,nn,jerr)
      if(jerr.ne.0) then
        ierr = jerr+1
        return
      endif
      if(nn.eq.0)   then
        ierr = 1
        return
      endif
C--   Loop over words in stin
      do i = 1,nn
        i1 = iw1(i)
        i2 = iw2(i)
C--     Classify word
        call sfmtStype(stin(i1:i2),ctype,field,nw,nd)
C--     Reformat word
        call sfmtRefmt(field,stin(i1:i2),wordout,nw,nd)
C--     Too large or zero length word encountered
        if(nw.gt.lenw .or. nw.eq.0) then
          ierr = 4
          return
        endif
C--     Get format descriptor
        call sfmtSform(field,nw,nd,cformat,nf)
C--     Assemble stout and fmt
        if(nn.eq.1) then
          call sfmtPutwd(' '    , wordout, stout, ierr)
          if(ierr.ne.0) goto 100
          call sfmtPutwd('( 1X,', cformat, fmt  , ierr)
          if(ierr.ne.0) goto 200
          call sfmtPutwd(' '    , ')'    , fmt  , ierr)
          if(ierr.ne.0) goto 200
        elseif(i.eq.1) then
          call sfmtPutwd(' '    , wordout, stout, ierr)
          if(ierr.ne.0) goto 100
          call sfmtPutwd('( 1X,', cformat, fmt  , ierr)
          if(ierr.ne.0) goto 200
        elseif(i.eq.nn) then
          call sfmtPutwd(' '    , wordout, stout, ierr)
          if(ierr.ne.0) goto 100
          call sfmtPutwd(',1X,' , cformat, fmt  , ierr)
          if(ierr.ne.0) goto 200
          call sfmtPutwd(' '    , ')'    , fmt  , ierr)
          if(ierr.ne.0) goto 200
        else
          call sfmtPutwd(' '    , wordout, stout, ierr)
          if(ierr.ne.0) goto 100
          call sfmtPutwd(',1X,' , cformat, fmt  , ierr)
          if(ierr.ne.0) goto 200
        endif
C--   End of loop over words in stin
      enddo

      return
 100  ierr = 5
      return
 200  ierr = 6
      return

      end

C     ====================================================
      subroutine sfmtPutwd(separator, word, cstring, ierr)
C     ====================================================

C--   Append separator//word to cstring
C--
C--   separator  (in) : separator like ' ' or ',' or ',1X,'
C--   word       (in) : input word
C--   cstring (inout) : string with word to be added
C--   ierr      (out) : 1 = nothing added since not enough space
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) separator, word, cstring

      ierr = 0
      mx   = len(cstring)
      nc   = imb_lenoc(cstring)
      nw   = imb_lenoc(word)
      ns   = max(imb_lenoc(separator),1)

      if(nc+nw+ns.gt.mx) then
        ierr = 1
        return
      endif

      cstring(nc+1:)    = separator
      cstring(nc+1+ns:) = word

      return
      end

C     ==================================================
      subroutine sfmtRefmt(field, strin, strout, nw, nd)
C     ==================================================

C--   Bring numbers into standard format
C--   The input numbers should be of the correct type
C--
C--   field   (in) : field descriptor L, I, F, D, E, A
C--   strin   (in) : input number of correct type
C--   strout (out) : reformatted number
C--   nw     (out) : field width
C--   nd     (out) : decimal field width
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) strin, strout
      character*10  fmt           !temporary format to read exponent
      character*1   field, kwot
      data kwot /''''/

      nn = imb_lenoc(strin)

      if(field.eq.'A' .and.
     +   strin(1:1).eq.kwot .and. strin(nn:nn).eq.kwot) then
C--     Quoted string: strip off the quotes
        strout = strin(2:nn-1)
        nw     = nn-2
        nd     = 0

      elseif(field.eq.'L' .or. field.eq.'A') then
C--     Just copy
        strout = strin
        nw     = nn
        nd     = 0

      elseif(field.eq.'I') then
C--     Here leading zeros are stripped from an integer e.g. 002 --> 2
C--     Handle sign
        if(strin(1:1) .eq. '+') then
          strout(1:1) = '+'
          ipos = 2
          jpos = 2
        elseif(strin(1:1) .eq. '-') then
          strout(1:1) = '-'
          ipos = 2
          jpos = 2
        else
          ipos = 1
          jpos = 1
        endif
C--     Skip leading zeros
        do while(strin(ipos:ipos).eq.'0' .and. ipos.lt.nn )
          ipos = ipos+1
        enddo
C--     Now copy
        strout(jpos:) = strin(ipos:nn)
C--     Field descriptors
        nw = imb_lenoc(strout)
        nd = 0

      elseif(field.eq.'F') then
C--     Here leading zeros are stripped from a floating point number.
C--     A leading (trailing) decimal point will we peceeded (followed)
C--     by a zero, e.g. .12 --> 0.12 and 5. --> 5.0
        idot = index(strin,'.')
C--     Handle sign
        if(strin(1:1) .eq. '+') then
          strout(1:1) = '+'
          ipos = 2
          jpos = 2
        elseif(strin(1:1) .eq. '-') then
          strout(1:1) = '-'
          ipos = 2
          jpos = 2
        else
          ipos = 1
          jpos = 1
        endif
C--     Skip leading zeros
        do while(strin(ipos:ipos).eq.'0')
          ipos = ipos+1
        enddo
C--     Add leading zero
        if(idot.eq.ipos) then
          strout(jpos:jpos) = '0'
          jpos = jpos+1
        endif
C--     Copy rest of string
        strout(jpos:)  = strin(ipos:nn)
        jpos = imb_lenoc(strout)+1
C--     Check for trailing decimal point
        if(strin(nn:nn).eq.'.') then
          strout(jpos:jpos) = '0'
          jpos = jpos+1
        endif
C--     Now get the new field descriptors
        nw   = imb_lenoc(strout)
        idot = index(strout,'.')
        nd   = nw-idot

      elseif(field.eq.'D' .or. field.eq.'E') then
*         write(6,*) 'strin---',strin
C--      Now for the E and D formats
         iddee = index(strin,'D')
         if(iddee .eq. 0) iddee = index(strin,'E')
         if(iddee .eq. 0) iddee = index(strin,'d')
         if(iddee .eq. 0) iddee = index(strin,'e')
         idot = index(strin,'.')
         nn   = imb_lenoc(strin)
         ipos = 1
         jpos = 1
C--      Handle sign
         if(strin(1:1) .eq. '+') then
           strout(1:3) = '+0.'
           ipos        = 2
           jpos        = 4
         elseif(strin(1:1) .eq. '-') then
           strout(1:3) = '-0.'
           ipos        = 2
           jpos        = 4
         else
           strout(1:2) = '0.'
           ipos        = 1
           jpos        = 3
         endif
*        write(6,*) 'strout1---',strout
C--     Skip leading zeros
        do while(strin(ipos:ipos).eq.'0')
          ipos = ipos+1
        enddo
C--     Number of digits before the decimal point
        if(idot.ne.0) then
          nshft = idot - ipos
        else
          nshft = iddee - ipos
        endif
        if(nshft .ne. 0) then
          strout(jpos:jpos+nshft-1) = strin(ipos:ipos+nshft-1)
          ipos = ipos+nshft
          jpos = jpos+nshft
        endif
*        write(6,*) 'strout2---',strout
        if(idot.ne.0) then
C--       ipos is now at the position of the dot so skip
          if(ipos.ne.idot) stop 'sfmtRefmt: ipos not at idot'
          ipos = ipos+1
        endif
C--     copy the remaining digits upto E or D
        nleft = iddee-ipos
        if(nleft.ne.0) then
          strout(jpos:jpos+nleft-1) = strin(ipos:ipos+nleft-1)
          ipos = ipos+nleft
          jpos = jpos+nleft
        endif
*        write(6,*) 'strout3---',strout
C--     ipos is now at the position of E or D so skip
        if(ipos.ne.iddee) stop 'sfmtRefmt: ipos not at iddee'
        ipos = ipos+1
C--     write the field descriptor on the output string
        strout(jpos:jpos) = field
        jpos = jpos+1
*        write(6,*) 'strout4---',strout
C--     read the exponent including sign
        nexp = nn-ipos+1
        fmt(1:2) = '(I'
        call smb_itoch(nexp,fmt(3:),ll)
        ll = imb_lenoc(fmt)
        fmt(ll+1:ll+1) = ')'
        read(strin(ipos:),fmt) iexp
C--     Update exponent
        iexp = iexp+nshft
C--     Write new exponent
        call smb_itoch(iexp,strout(jpos:),ll)
*        write(6,*) 'strout5---',strout
C--     Now get the new field descriptors
        nw    = imb_lenoc(strout)
        ieedd = index(strout,field)
        idot  = index(strout,'.')
        nd    = ieedd-idot-1

      else
        stop 'sfmtRefmt: cannot handle this type'
      endif

      return
      end

C     ================================================
      subroutine sfmtSform(field, nw, nd, cformat, nf)
C     ================================================

C--   Make format string
C--
C--   field    (in) : field descriptor L, I, F, D, E, A
C--   nw       (in) : field width
C--   nd       (in) : decimal width
C--   cformat (out) : string with format descriptor
C--   nf      (out) : number of chars in cformat
C--
C--   It is not checked that the format descriptor fits into cformat
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cformat
      character*1   field

      if(field.eq.'L' .or. field.eq.'I' .or. field .eq. 'A' ) then
        cformat(1:1) = field
        call smb_itoch(nw,cformat(2:),nn)
        nf = imb_lenoc(cformat)
      else
        cformat(1:1) = field
        call smb_itoch(nw,cformat(2:),nn)
        cformat(2+nn:2+nn) = '.'
        call smb_itoch(nd,cformat(3+nn:),nn)
        nf = imb_lenoc(cformat)
      endif

      return
      end

C     ===================================================
      subroutine sfmtStype(cstring, ctype, field, nw, nd)
C     ===================================================

C--   Determine type of string
C--
C--   cstring  (in)  string without embedded, leading or trailing blanks
C--   ctype   (out)  LOGI, INTE, REAL, CHAR, VOID
C--   field   (out)  Field descriptor L, I, F, D, E, A, or ' ' = empty
C--   nw      (out)  Field width
C--   nd      (out)  Decimal filed width (irrelevant for I and A)
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15


      implicit double precision (a-h,o-z)

      character*(*) cstring
      character*4   ctype
      character*1   fmt, field
      logical sfmtInte, sfmtReal, sfmtDorE

      n = imb_lenoc(cstring)
      if(n.eq.0) then
        ctype = 'VOID'
        field = ' '
        nw    = 0
        nd    = 0
        return
      endif

C--   Preclassify
      i = ifmtDEFIC(cstring,fmt)

      if(fmt .eq. 'L') then
        ctype = 'LOGI'
        field = 'L'
        nw = n
        nd = 0
        return
      elseif(fmt .eq. 'C') then
        ctype = 'CHAR'
        field = 'A'
        nw = n
        nd = 0
        return
      elseif(fmt .eq. 'Q') then
        ctype = 'CHAR'
        field = 'A'
        nw = n
        nd = 0
        return
      elseif(fmt. eq. 'I' .and. sfmtInte(cstring,n)) then
        ctype = 'INTE'
        field = 'I'
        nw    = n
        nd    = 0
        return
      elseif(fmt. eq. 'F' .and. sfmtReal(cstring,i)) then
        ctype = 'REAL'
        field = 'F'
        nw    =  n
        nd    =  n-i
        return
      elseif(fmt. eq. 'D' .and. sfmtDorE(cstring,i,idot)) then
        ctype = 'REAL'
        field = 'D'
        nw    = n
        nd    = 0
        if(idot.ne.0) nd = i-idot-1
        return
      elseif(fmt. eq. 'E' .and. sfmtDorE(cstring,i,idot)) then
        ctype = 'REAL'
        field = 'E'
        nw    =  n
        nd    =  0
        if(idot.ne.0) nd = i-idot-1
        return
      else
        ctype = 'CHAR'
        field = 'A'
        nw    =  n
        nd    =  0
        return
      endif

      return
      end

C     ============================================
      integer function ifmtDEFIC(cstring, cformat)
C     ============================================

C--   Preclassifies string
C--
C--   cstring  (in)  string without embedded, leading or trailing blanks
C--   cformat (out) 'D' contains a 'D' and no 'E' return position of 'D'
C--                 'E' contains a 'E' and no 'D' return position of 'E'
C--                 'F' contains a '.' and no 'D' or 'E' return pos of '.'
C--                 'Q' quoted string return 1
C--                 'L' single T or F in upper case return 1
C--                 'C' defenitely a character string return 1
C--                 'I' not any of the above return 1
C--                 ' ' empty string return 0
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cstring
      character*1   cformat, kwot
      data kwot /''''/

      n    = imb_lenoc(cstring)

      if(n.eq.0) then
        cformat   = ' '
        ifmtDEFIC = 0
        return
      elseif(n.eq.1 .and. (cstring.eq.'T' .or. cstring.eq.'F')) then
        cformat   = 'L'
        ifmtDEFIC = 1
        return
      elseif(n.eq.1) then
        cformat   = 'I'
        ifmtDEFIC = 1
        return
      elseif(cstring(1:1).eq.kwot .and. cstring(n:n).eq.kwot) then
        cformat   = 'Q'
        ifmtDEFIC = 1
        return
      endif

      idot = index(cstring,'.')
      iddd = index(cstring,'d')
      if(iddd.eq.0) iddd = index(cstring,'D')
      ieee = index(cstring,'e')
      IF(ieee.eq.0)ieee = index(cstring,'E')

C--   Find format
      if(idot.ne.0 .and. iddd.eq.0 .and. ieee.eq.0) then
        cformat = 'F'
        ifmtDEFIC = idot
        return
      elseif(ieee.eq.0 .and. iddd.gt.1 .and. iddd.lt.n) then
        cformat = 'D'
        ifmtDEFIC = iddd
        return
      elseif(iddd.eq.0 .and. ieee.gt.1 .and. ieee.lt.n) then
        cformat = 'E'
        ifmtDEFIC = ieee
        return
      elseif(ieee.eq.1 .or. ieee.eq.n .or.
     +       iddd.eq.1 .or. iddd.eq.n ) then
        cformat = 'C'
        ifmtDEFIC = 1
        return
      else
        cformat = 'I'
        ifmtDEFIC = 1
        return
      endif

      return
      end

C     ===========================================
      logical function sfmtDorE(cstring, i, idot)
C     ===========================================

C--   True if string is a D or E format real number
C--   The string must have been preclassified with ifmtDEFIC
C--
C--   cstring  (in) string without embedded, leading or trailing blanks
C--   i        (in) position of D or E
C--   idot    (out) position of decimal dot (0 = no dot)
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cstring

      logical sfmtInte, sfmtReal

      n    = imb_lenoc(cstring)
C--   Empty string
      if(n.eq.0 .or. i.le.0 .or. i.gt.n) then
        sfmtDorE = .false.
        return
      endif

      idot = index(cstring(1:i-1),'.')

C--   Must be signed integer if no dot before E
      if(idot.eq.0 .and. .not.sfmtInte(cstring(1:i-1),i-1)) then
        sfmtDorE = .false.
        return
      endif
C--   Must be real if dot before E
      if(idot.ne.0 .and. .not.sfmtReal(cstring(1:i-1),idot)) then
        sfmtDorE = .false.
        return
      endif
C---  Must be signed integer after E
      if(i.lt.n .and. .not.sfmtInte(cstring(i+1:n),n-i)) then
        sfmtDorE = .false.
      return
      endif

      sfmtDorE = .true.

      return
      end

C     =====================================
      logical function sfmtReal(cstring, i)
C     =====================================

C--   True if string is a F format real number
C--   The string must have been preclassified with ifmtDEFIC
C--
C--   cstring  (in) string without embedded, leading or trailing blanks
C--   i        (in) position of decimal point
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cstring

      logical sfmtInte, sfmtUint

      n = imb_lenoc(cstring)
C--   Empty string
      if(n.eq.0 .or. i.le.0 .or. i.gt.n) then
        sfmtReal = .false.
C--   Not signed integer before the dot
      elseif(i.gt.1 .and. .not.sfmtInte(cstring(1:i-1),i-1) ) then
        sfmtReal = .false.
C---  Not unsigned integer after dot
      elseif(i.lt.n .and. .not.sfmtUint(cstring(i+1:n),n-i) ) then
        sfmtReal = .false.
      else
        sfmtReal = .true.
      endif

      return
      end

C     =====================================
      logical function sfmtInte(cstring, n)
C     =====================================

C--   True if string is an integer (no preclassification needed)
C--
C--   cstring  (in) string without embedded, leading or trailing blanks
C--   n        (in) number of characters in string
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cstring
      character*12  charset
C--                  123456789012
      data charset /'+-1234567890'/

      logical found

      if(n.le.0) stop 'sfmtInte: invalid string length'

      sfmtInte = .true.
      do i = 1,n
        j = 3
        if(i.eq.1) j = 1
        found = .false.
        do while(.not.found .and. j.le.12)
          if(cstring(i:i) .eq. charset(j:j)) found = .true.
          j = j+1
        enddo
        if(.not.found) then
          sfmtInte = .false.
          return
        endif
      enddo

      return
      end

C     =====================================
      logical function sfmtUInt(cstring, n)
C     =====================================

C--   True if string is an unsigned integer (no preclassification needed)
C--
C--   cstring  (in) string without embedded, leading or trailing blanks
C--   n        (in) number of characters in string
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cstring
      character*10  charset
C--                  1234567890
      data charset /'1234567890'/

      logical found

      if(n.le.0) stop 'sfmtUint: invalid string length'

      sfmtUint = .true.
      do i = 1,n
        j = 1
        found = .false.
        do while(.not.found .and. j.le.10)
          if(cstring(i:i) .eq. charset(j:j)) found = .true.
          j = j+1
        enddo
        if(.not.found) then
          sfmtUint = .false.
          return
        endif
      enddo

      return
      end

C     =========================================================
      subroutine sfmtParseit(idim, cstring, iw1, iw2, nw, ierr)
C     =========================================================

C--   Parse string into words (quotes are dealt with in sfmtGetWord)
C--
C--   idim     (in) dimension of iw1 and iw2
C--   cstring  (in) character string
C--   iw1(i)  (out) start of word i
C--   iw2(i)  (out) end of word i
C--   nw      (out) number of words  0=empty string
C--   ierr    (out) 1 = unbalanced quotes in string
C--                 2 = word count exceeded
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cstring
      dimension iw1(*), iw2(*)

      m1   = imb_frstc(cstring)
      m2   = imb_lenoc(cstring)
      nw   = 0
      ierr = 0
      if(m1.eq.0) then
        nw = 0
      else
        i1 = m1
        k1 = m1
        do while(k1 .ne. 0)
          call sfmtGetWord(cstring,i1,k1,k2,ierr)
          if(ierr.eq.1) return
          if(k1  .eq.0) return
          if(ierr.ne.2) then
            nw = nw+1
            if(nw.gt.idim) then
              ierr = 2
              return
            endif
            iw1(nw) = k1
            iw2(nw) = k2
          endif
          i1   = k2+1
          ierr = 0
        enddo
      endif

      return
      end

C     =================================================
      subroutine sfmtGetWord(cstring, i1, k1, k2, ierr)
C     =================================================

C--   Get limits of first word in string(i1:)
C--   These limits are the first and last nonblank characters of
C--   a continuous substring of nonblank characters; for quoted substrings
C--   the openening and closing quote are taken as the limits
C--
C--   cstring    (in) character string
C--   i1         (in) start of search range
C--   k1        (out) position of first character of word 0=no word
C--   k2        (out) position of last  character of word 0=no word
C--   ierr      (out) 1 = unbalanced quotes
C--                   2 = empty quoted string
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cstring
      character*1   kwot
      data kwot  /''''/
      character*2   kwotb
      data kwotb /''' '/

      k1   = 0
      k2   = 0
      ierr = 0
      m2   = imb_lenoc(cstring)

C--   Empty string
      if(m2.eq.0)  return

C--   Search out of range
      if(i1.gt.m2) return

C--   Start of first word in cstring(i1:)
      k1 = ifmtFstChar(cstring,i1)

C--   Empty cstring(i1:)
      if(k1.eq.0) return

      if(k1.eq.m2 .and. cstring(k1:k1).eq.kwot) then
C--     1-character string with quote --> unbalanced quote
        ierr = 1

      elseif(cstring(k1:k1).eq.kwot) then
C--     Found opening quote, search for closing quote
        iq = 0
        do i = k1+1,m2-1
          if(iq.eq.0 .and. cstring(i:i+1).eq.kwotb) iq = i
        enddo
        if(iq.eq.0 .and. cstring(m2:m2).eq.kwot)    iq = m2
C--     Unbalanced quote
        if(iq.eq.0) then
          k1   = 0
          k2   = 0
          ierr = 1
          return
        endif
        k2 = iq
C--     Flag empty quoted string
        if(k2.eq.k1+1) ierr = 2

      else
C--     Handle non quoted string
        k2 = ifmtFstEofW(cstring,i1)

      endif

      return
      end

C     ========================================
      integer function ifmtFstChar(cstring, i)
C     ========================================

C--   Position of first non-blank character in cstring(i:)
C--   Returns zero when the string is empty
C--
C--   cstring  (in) : input string
C--   i        (in) : start of search range
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      implicit double precision (a-h,o-z)

      character*(*) cstring

      ii = imb_frstc(cstring(i:))
      if(ii.eq.0) then
        ifmtFstChar = 0
      else
        ifmtFstChar = i+ii-1
      endif

      return
      end

C     =========================================
      integer function ifmtFstEofW(cstring, i1)
C     =========================================

C--   Find position of the first end-of-word in string(i1:)
C--   end-of-word = nonblank-blank or nonblank at end of string
C--   Returns zero when the string is empty or eow not found
C--
C--   cstring  (in) : input string
C--   i1       (in) : start of search range
C--
C--   Author: Michiel Botje h24@nikhef.nl   23-03-15

      character*(*) cstring

      ifmtFstEofW = 0

      m1 = imb_frstc(cstring)
      m2 = imb_lenoc(cstring)
      if(m2.eq.0) return
      j1 = max(i1,m1)
      j2 = m2

      do i = j1,j2
        if(cstring(i:i).ne.' ') then
          if(i.eq.m2) then
            ifmtFstEofW = i
            return
          elseif(cstring(i+1:i+1).eq.' ') then
            ifmtFstEofW = i
            return
          endif
        endif
      enddo

      return
      end

