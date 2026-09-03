
C--   This is the file mbits.f with the MBUTIL bitwise routines
C--
C--   subroutine smb_sbit1(i,n)
C--   subroutine smb_sbit0(i,n)
C--   integer function imb_gbitn(i,n)
C--   subroutine smb_cbyte(i1,ibyte1,i2,ibyte2)
C--   integer function imb_sbits(cpatt)
C--   subroutine smb_gbits(i,cpatt)
C--   integer function imb_test0(mask,i)
C--   integer function imb_test1(mask,i)
C--   integer function imbAllZero()
C--   integer function imbAllOne()
C--   integer function imbGimmeOne()

C     =========================
      subroutine smb_sbit1(i,n)
C     =========================

C--   Set bit n of i to one
C--
C--   Author: Michiel Botje h24@nikhef.nl   26-03-19

      implicit double precision (a-h,o-z)

      if(n.ge.1 .and. n.le.32) then
        j = ishft(imbGimmeOne(), n-1)
        i = ior( i, j )
      else
        i = imbAllOne()
      endif

      return
      end

C     =========================
      subroutine smb_sbit0(i,n)
C     =========================

C--   Set bit n of i to zero
C--
C--   Author: Michiel Botje h24@nikhef.nl   26-03-19

      implicit double precision (a-h,o-z)

      if(n.ge.1 .and. n.le.32) then
        j = ishft(imbGimmeOne(),n-1)
        i = iand( i, not( j ) )
      else
        i = imbAllZero()
      endif

      return
      end


C     ===============================
      integer function imb_gbitn(i,n)
C     ===============================

C--   Get bit n of i
C--
C--   Author: Michiel Botje h24@nikhef.nl   26-03-19

      implicit double precision (a-h,o-z)

      if(n.ge.1 .and. n.le.32) then
        j         = ishft(i,1-n)
        imb_gbitn = iand( j,imbGimmeOne() )
      else
        imb_gbitn = -1
      endif

      return
      end

C     =========================================
      subroutine smb_cbyte(i1,ibyte1,i2,ibyte2)
C     =========================================

C--   Copy ibyte1 of i1 to ibyte2 of i2
C--
C--   Author: Michiel Botje h24@nikhef.nl   26-03-19

      implicit double precision (a-h,o-z)

      dimension ib1(4), ib2(4)
      data ib1/1,9,17,25/, ib2/8,16,24,32/

      if(ibyte1.lt.1 .or. ibyte1.gt.4) return
      if(ibyte2.lt.1 .or. ibyte2.gt.4) return

      mask  = ishft( imbAllOne(), -24 )
      j     = ishft( i1, 1-ib1(ibyte1) )
      ibyte = iand( j, mask )
      j     = ishft( mask, ib1(ibyte2)-1 )
      mask  = not( j )
      i2    = iand( i2, mask )
      j     = ishft( ibyte, ib1(ibyte2)-1 )
      i2    = ior( i2, j )

      return
      end

C     =================================
      integer function imb_sbits(cpatt)
C     =================================

C--   Set bitpattern of 32-bit integer i
C--
C--   cpatt   (in) : character*32 string containing the bitpattern
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cpatt

      leng = len(cpatt)
      if(leng.lt.32) stop 'IMB_SBITS: input string < 32 characters'

      imb_sbits = 0      !avoid compiler warning

      do j = 1,32
        k = 33-j
        if(cpatt(k:k).eq.'0') then
          call smb_sbit0(imb_sbits,j)
        else
          call smb_sbit1(imb_sbits,j)
        endif
      enddo

      return
      end

C     =============================
      subroutine smb_gbits(i,cpatt)
C     =============================

C--   Get bitpattern of 32-bit integer i
C--
C--   i      (in): input 32-bit integer i
C--   cpatt (out): character*32 string with the bitpattern of i
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cpatt

      leng = len(cpatt)
      if(leng.lt.32) stop 'SMB_GBITS: output string < 32 characters'
      call smb_cfill(' ',cpatt)

      do j = 1,32
        k = 33-j
        l = imb_gbitn(i,j)
        if(l.eq.0) then
          cpatt(k:k) = '0'
        else
          cpatt(k:k) = '1'
        endif
      enddo

      return
      end

C     ===================================
      subroutine smb_bytes(cbitsi,cbitso)
C     ===================================

C--   Insert blanks at the bite boundaries of bitpattern cbitsi
C--
C--   Author: Michiel Botje h24@nikhef.nl   26-03-19

      implicit double precision (a-h,o-z)

      character*(*) cbitsi, cbitso
      character*36  cbits

      leng = len(cbitsi)
      if(leng.lt.32)  stop 'SMB_BYTES: input string < 32 characters'
      leng = len(cbitso)
      if(leng.lt.35)  stop 'SMB_BYTES: output string < 35 characters'

      ibiti = 0
      ibito = 0
      do ibyte = 1,4
        do ibit = 1,8
          ibiti = ibiti+1
          ibito = ibito+1
          cbits(ibito:ibito) = cbitsi(ibiti:ibiti)
        enddo
        ibito = ibito+1
        cbits(ibito:ibito) = ' '
      enddo

      cbitso(1:35) = cbits(1:35)
      do i = 36,leng
        cbitso(i:i) = ' '
      enddo

      return
      end

C     ==================================
      integer function imb_test0(mask,i)
C     ==================================

C--   Test pattern of 'zero' bits in 32-bit integer i
C--
C--   Input:  integer mask with bit n = 0(1) -> ignore(test) bit n of i
C--           integer i containing the bitpattern to be tested
C--   Output: integer imb_test0 = 0 if all tested bits of i are 0
C--                             # 0 if not all tested bits of i are 0
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      imb_test0 = iand(mask,i)

      return
      end

C     ==================================
      integer function imb_test1(mask,i)
C     ==================================

C--   Test pattern of 'one' bits in 32-bit integer i
C--
C--   Input:  integer mask with bit n = 0(1) -> ignore(test) bit n of i
C--           integer i containing the bitpattern to be tested
C--   Output: integer imb_test1 = 0 if all tested bits of i are 1
C--                             # 0 if not all tested bits of i are 1
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      imb_test1 = iand(mask,not(i))

      return
      end

C     =============================
      integer function imbAllZero()
C     =============================

C--   Returns integer with all bits set to zero
C--
C--   Author: Michiel Botje h24@nikhef.nl   26-03-19

      implicit double precision (a-h,o-z)

      imbAllZero = iand(kdum,not(kdum))

      return
      end


C     ============================
      integer function imbAllOne()
C     ============================

C--   Returns integer with all bits set to one
C--
C--   Author: Michiel Botje h24@nikhef.nl   26-03-19

      implicit double precision (a-h,o-z)

      imbAllOne = ior(kdum,not(kdum))

      return
      end

C     ==============================
      integer function imbGimmeOne()
C     ==============================

C--   Returns integer with all bits zero except the first bit (LSB)
C--
C--   Author: Michiel Botje h24@nikhef.nl   26-03-19

      implicit double precision (a-h,o-z)

      j           = ishft(imbAllOne(),1)
      imbGimmeOne = not( j )

      return

      end

