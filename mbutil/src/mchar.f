
C--   This is the file mchar.f containing the MBUTIL character manipulation routines 
C--
C--   subroutine smb_cfill(char,cstring)
C--   subroutine smb_cleft(cstring)
C--   subroutine smb_crght(cstring)
C--   subroutine smb_cltou(cstring)
C--   subroutine smb_cutol(cstring)
C--   integer function imb_lastc(cstring)
C--   integer function imb_frstc(cstring)
C--   logical function lmb_compc(string1,string2,n1,n2)
C--   logical function lmb_match(string,substr,wdcard)
C--   subroutine smb_itoch(in,chout,leng)
C--   subroutine smb_dtoch(dd,n,chout,leng)
C--   subroutine smb_hcode(ihash,hcode)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  MBUTIL character routines from mchar.f                  */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  MBUTIL character routines from mchar.f                    */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  MBUTIL character routines from mchar.f                    */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

      
C     ==================================
      subroutine smb_cfill(char,cstring)
C     ==================================

C--   Input:  character string char.
C--           character string cstring.
C--           On exit all characters in cstring will be set to char(1:1).
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) char, cstring

      l = len(cstring)

      do i = 1,l
        cstring(i:i) = char(1:1)
      enddo

      return
      end

C     =============================
      subroutine smb_cleft(cstring)
C     =============================

C--   Left adjust character string cstring.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring

      max         = len(cstring)
      if(max.le.0)  return
      i1          = imb_frstc(cstring)
      i2          = imb_lastc(cstring)
      k           = 0
      do i = i1,i2
        k = k+1
        cstring(k:k) = cstring(i:i)
      enddo
      do i = k+1,max
        cstring(i:i) = ' '
      enddo

      return
      end

C     =============================
      subroutine smb_crght(cstring)
C     =============================

C--   Right adjust character string cstring.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring

      max         = len(cstring)
      if(max.le.0)  return
      i1          = imb_frstc(cstring)
      i2          = imb_lastc(cstring)
      k           = max+1
      do i = i2,i1,-1
        k = k-1
        cstring(k:k) = cstring(i:i)
      enddo
      do i = k-1,1,-1
        cstring(i:i) = ' '
      enddo

      return
      end

C     =============================
      subroutine smb_cltou(cstring)
C     =============================

C--   Input: character string cstring.
C--          On exit all characters in cstring will be set
C--          to upper case.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring
      character*26  charl, charu
      data charl   /'abcdefghijklmnopqrstuvwxyz'/
      data charu   /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      j          = len(cstring)

      do i = 1,j
        do k = 1,26
          if(cstring(i:i).eq.charl(k:k)) then
            cstring(i:i) = charu(k:k)
          endif
        enddo
      enddo

      return
      end

C     =============================
      subroutine smb_cutol(cstring)
C     =============================

C--   Input: character string cstring.
C--          On exit all characters in cstring will be set
C--          to lower case.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring
      character*26  charl, charu
      data charl   /'abcdefghijklmnopqrstuvwxyz'/
      data charu   /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      j          = len(cstring)

      do i = 1,j
        do k = 1,26
          if(cstring(i:i).eq.charu(k:k)) then
            cstring(i:i) = charl(k:k)
          endif
        enddo
      enddo

      return
      end

C     ===================================
      integer function imb_lastc(cstring)
C     ===================================

C--   Input:  cstring       Character string.
C--   Output: imb_lastc     Position of the last non-blank
C--                         character in cstring.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring

      j         = len(cstring)
      imb_lastc = 0

      do i = j,1,-1
        if(cstring(i:i).ne.' ') then
          imb_lastc = i
          return
        endif
      enddo

      return
      end

C     ===================================
      integer function imb_lenoc(cstring)
C     ===================================

C--   Function name kept for backward compatibility

      character*(*) cstring

      imb_lenoc = imb_lastc(cstring)

      return
      end

C     ===================================
      integer function imb_frstc(cstring)
C     ===================================

C--   Input:  cstring       Character string.
C--   Output: imb_frstc     Position of the first non-blank
C--                         character in cstring.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring

      j         = len(cstring)
      imb_frstc = 0

      do i = 1,j
        if(cstring(i:i).ne.' ') then
          imb_frstc = i
          return
        endif
      enddo

      return
      end

C     =================================================
      logical function lmb_compc(string1,string2,n1,n2)
C     =================================================

C--   Routine kept for backward compatibility
C--   Do a case independent comparision of the characters
C--   string1(n1:n2) and string2(n1:n2).
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) string1, string2
      character*1   ch1    , ch2

      lmb_compc = .false.
      if((n1.le.0).or.(n2.le.0).or.(n2.lt.n1)) return
      len1 = imb_lenoc(string1)
      if(len1.lt.n2)                           return
      len2 = imb_lenoc(string2)
      if(len2.lt.n2)                           return

      lmb_compc = .true.
      do i = n1,n2
        ch1 = string1(i:i)
        ch2 = string2(i:i)
        call smb_cltou(ch1)
        call smb_cltou(ch2)
        if(ch1.ne.ch2) goto 900
      enddo

      return

 900  continue

      lmb_compc = .false.

      return
      end

C     ==================================================
      logical function lmb_comps(string1,string2,istrip)
C     ==================================================

C--   Do a case independent comparision of string1 and string2
C--   with or without stripping
C--
C--   Author: Michiel Botje h24@nikhef.nl   15-12-18

      character*(*) string1, string2
      character*1   ch1    , ch2

       i1 = imb_frstc(string1)
       i2 = imb_frstc(string2)
       j1 = imb_lastc(string1)
       j2 = imb_lastc(string2)

      if(j1.eq.0 .and. j2.eq.0)      then
        lmb_comps = .true.
      elseif( (j1-i1).ne.(j2-i2) )   then
        lmb_comps = .false.
      else
        if(istrip.eq.0) then
          i1 = 1
          i2 = 1
        endif
        lmb_comps = .true.
        k1 = i1-1
        k2 = i2-1
        do while(lmb_comps.and.k1.le.j1)
          k1  = k1+1
          k2  = k2+1
          ch1 = string1(k1:k1)
          ch2 = string2(k2:k2)
          call smb_cltou(ch1)
          call smb_cltou(ch2)
          if(ch1.ne.ch2) lmb_comps = .false.
        enddo
      endif

      return
      end

C     ================================================
      logical function lmb_match(string,substr,wdcard)
C     ================================================

C--   Translates string, substr to upper case and returns
C--   .true. if substr is contained in string, .false. otherwise. 
C--   If string and/or substr are null strings, lmb_match = .FALSE.
C--   Leading blanks are stripped off.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09
C--   Modified:                             15-12-18

      character*(*) string,substr,wdcard
      character*80  str,sub,sbb
      character*1   wcd

      lmb_match = .false.
      ipos      = 0

      len1 = imb_lastc(string)
      if(len1.eq.0.or.len1.gt.80)   return
      len2 = imb_lastc(substr)
      if(len2.eq.0.or.len2.gt.80)   return

C--   Avoid modifying the input argument(s)
      call smb_cfill(' ',str)
      call smb_cfill(' ',sub)
      str(1:len1) = string(1:len1)
      sub(1:len2) = substr(1:len2)
      wcd         = wdcard(1:1)
      call smb_cltou(str)
      call smb_cltou(sub)
      call smb_cltou(wcd)
      call smb_cleft(str)
      call smb_cleft(sub)
      len1 = imb_lastc(str)
      len2 = imb_lastc(sub)
      if(len2.gt.len1) return

      i2 = len1-len2+1
      do i = 1,i2
        iend = i-1+len2
        do k = i,iend
          l = k-i+1
          sbb(l:l) = sub(l:l)
          if(wcd.ne.' '.and.sub(l:l).eq.wcd) sbb(l:l) = str(k:k)
        enddo
        if(str(i:iend).eq.sbb(1:len2)) lmb_match = .true.
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_itoch(int in, string &chout, int &leng);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_itochcpp FC_FUNC(smb_itochcpp,SMB_ITOCHCPP)
CXXHFW    void fsmb_itochcpp(int*, char*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_itoch(int in, string &chout, int &leng)
CXXWRP  {
CXXWRP    int ls = 20;
CXXWRP    char *mychar = new char[ls+1];
CXXWRP    fsmb_itochcpp(&in,mychar,&ls,&leng);
CXXWRP    chout = "";
CXXWRP    for(int i = 0; i < leng; i++) {
CXXWRP      chout = chout + mychar[i];
CXXWRP    }
CXXWRP    delete[] mychar;
CXXWRP  }
C-----------------------------------------------------------------------

C     =========================================
      subroutine smb_itochCPP(in,chout,ls,leng)
C     =========================================

C--   Wrapper for smb_itoch
C--
C--   Author: Michiel Botje h24@nikhef.nl   11-08-21

      implicit double precision (a-h,o-z)

      character*(100) chout

      if(ls.gt.100) stop 'smb_itochCPP: output chout > 100 characters'

      call smb_itoch(in,chout(1:ls),leng)

      return
      end
      
C     ===================================
      subroutine smb_itoch(in,chout,leng)
C     ===================================

C--   Format integer as character string of size leng
C--
C--   in       (in)    input integer
C--   chout   (out)    character variable
C--   leng    (out)    length of character string in chout
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision(a-h,o-z)
      
      character chout*(*), str*30

C--   Initialise
      call smb_cfill(' ',chout)
      lmax = len(chout)
C--   Write integer format into large enough string
      write(str,'(I30)') in
      i1   = imb_frstc(str)
      i2   = imb_lastc(str)
      leng = i2-i1+1
C--   Copy string into chout
      if(leng.gt.lmax) then
         call smb_cfill('*',chout)
         leng = lmax
      else
        chout(1:leng) = str(i1:i2)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_dtoch(double dd, int n, string &chout, int &leng);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_dtochcpp FC_FUNC(smb_dtochcpp,SMB_DTOCHCPP)
CXXHFW    void fsmb_dtochcpp(double*, int*, char*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_dtoch(double dd, int n, string &chout, int &leng)
CXXWRP  {
CXXWRP    int ls = 20;
CXXWRP    char *mychar = new char[ls+1];
CXXWRP    fsmb_dtochcpp(&dd,&n,mychar,&ls,&leng);
CXXWRP    chout = "";
CXXWRP    for(int i = 0; i < leng; i++) {
CXXWRP      chout = chout + mychar[i];
CXXWRP    }
CXXWRP    delete[] mychar;
CXXWRP  }
C-----------------------------------------------------------------------

C     ===========================================
      subroutine smb_dtochCPP(dd,n,chout,ls,leng)
C     ===========================================

C--   Wrapper for smb_dtoch
C--
C--   Author: Michiel Botje h24@nikhef.nl   11-08-21

      implicit double precision (a-h,o-z)

      character*(100) chout

      if(ls.gt.100) stop 'smb_dtochCPP: output chout > 100 characters'

      call smb_dtoch(dd,n,chout(1:ls),leng)

      return
      end

C     =====================================
      subroutine smb_dtoch(dd,n,chout,leng)
C     =====================================

C--   Put floating point number in character string
C--   Use the tightest format possible
C--
C--   dd       (in)    input double precision number
C--   n        (in)    number of digits to be kept on output [1-9]
C--   chout   (out)    character variable
C--   leng    (out)    length of character string in chout
C--
C--   Remark: smb_MantEx (see hash.f) brings n into range [1-9]
C--
C--   Author: Michiel Botje h24@nikhef.nl   09-08-21

      implicit double precision (a-h,o-z)

      character*(*) chout
      character*20  cm, ct, fmt, str

C--   Initialise
      call smb_cfill(' ',chout)
      lmax = len(chout)
C--   Catch zero
      if(dd.eq.0.D0) then
        chout(1:1) = '0'
        leng       = 1
        return
      endif
C--   Get mantissa and exponent
      call smbMantEx(dd,n,im,ie)
*      write(6,'('' dtoch dd = '',E20.9)') dd
*      write(6,'('' nround   = '',I20  )') n
*      write(6,'('' mantissa = '',I20  )') im
*      write(6,'('' exponent = '',I20  )') ie
      call smb_itoch(im,cm,nm)
      isign = 0
      if(cm(1:1).eq.'-') isign = 1
      nonzero = nm-isign
      i       = nm
C--   Strip trailing zero's from mantissa
      do while(cm(i:i).eq.'0')
        nonzero = nonzero-1
        i       = i-1
      enddo
      ndigits = max(nonzero,ie)
*      write(6,'('' sign     = '',I20  )') isign
*      write(6,'('' nonzero  = '',I20  )') nonzero
*      write(6,'('' digits   = '',I20  )') ndigits
C--   Catch rounding
      call smbMantEx(dd,9,jm,je)
      call smb_itoch(jm,ct,nt)
      iround = 0
      do i = nm+1,nt
        if(ct(i:i).ne.'0') iround = 1
      enddo
*      write(6,'('' iround   = '',I20  )') iround
C--   Format output: iw is number of digits behind the dot in F format
      iw = nonzero-ie
      if(ndigits.le.ie .and. ie.le.6 .and. iround.eq.0) then
C--     Integer format
        write(str,'(I20)') int(dd)
        i1   = imb_frstc(str)
        i2   = imb_lastc(str)
*      elseif(iw.ge.0 .and. iw.le.ndigits) then
      elseif(iw.ge.0 .and. iw.le.5) then
C--     F format
        write(fmt,'(''(F20.'',I1,'')'')') iw
        write(str,fmt) dd
        i1   = imb_frstc(str)
        i2   = imb_lastc(str)
      else
C--     E format
        write(fmt,'(''(1PE20.'',I1,'')'')') nonzero-1
        write(str,fmt) dd
        i1   = imb_frstc(str)
        i2   = imb_lastc(str)
      endif

      leng = i2-i1+1
      if(leng.gt.lmax) then
        call smb_cfill('*',chout)
        leng = lmax
      else
        chout(1:leng) = str(i1:i2)
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_hcode(int ihash, string &hcode);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_hcodecpp FC_FUNC(smb_hcodecpp,SMB_HCODECPP)
CXXHFW    void fsmb_hcodecpp(int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_hcode(int ihash, string &hcode)
CXXWRP  {
CXXWRP    int ls = 15;
CXXWRP    char *mychar = new char[ls+1];
CXXWRP    fsmb_hcodecpp(&ihash,mychar,&ls);
CXXWRP    hcode = "";
CXXWRP    for(int i = 0; i < ls; i++) {
CXXWRP      hcode = hcode + mychar[i];
CXXWRP    }
CXXWRP    delete[] mychar;
CXXWRP  }
C-----------------------------------------------------------------------


C     =======================================
      subroutine smb_hcodeCPP(ihash,hcode,ls)
C     =======================================

C--   Extract the 4 byte-values of ihash and put them into hcode
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      character*(100) hcode

      if(ls.gt.100) stop 'smb_hcodeCPP: output hcode > 100 characters'

      call smb_hcode(ihash,hcode(1:ls))

      return
      end


C     =================================
      subroutine smb_hcode(ihash,hcode)
C     =================================

C--   Extract the 4 byte-values of ihash and put them into hcode
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      character*(*) hcode
      character*15  hash
      character*3   byte

      lhash = len(hcode)
      if(lhash.lt.15) then
        call smb_cfill('*',hcode)
        return
      else
        call smb_cfill(' ',hcode)
      endif

      jbyte = 0
      hash  = '000-000-000-000'

C--   Byte 4
      call smb_cbyte(ihash,4,jbyte,1)
      call smb_itoch(jbyte,byte,lbyte)
      hash(04-lbyte:03) = byte(1:lbyte)
C--   Byte 3
      call smb_cbyte(ihash,3,jbyte,1)
      call smb_itoch(jbyte,byte,lbyte)
      hash(08-lbyte:07) = byte(1:lbyte)
C--   Byte 2
      call smb_cbyte(ihash,2,jbyte,1)
      call smb_itoch(jbyte,byte,lbyte)
      hash(12-lbyte:11) = byte(1:lbyte)
C--   Byte 1
      call smb_cbyte(ihash,1,jbyte,1)
      call smb_itoch(jbyte,byte,lbyte)
      hash(16-lbyte:15) = byte(1:lbyte)

      hcode(1:15) = hash

      return
      end

