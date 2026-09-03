
C-- This is the file usrchecks.f with user interface checks
C--
C--   integer function iqcChkLmij(subnam,w,id,ix,iq,ichk)
C--   integer function iqcIxInside(subnam,ix,ichk)
C--   integer function iqcIqInside(subnam,iq,ichk)
C--   double precision function dqcXInside(subnam,xx,ichk)
C--   double precision function dqcQInside(subnam,qq,ichk)
C--
C--   subroutine sqcIlele(subnam,parnam,imin,ival,imax,comment)
C--   subroutine sqcIlelt(subnam,parnam,imin,ival,imax,comment)
C--   subroutine sqcIltle(subnam,parnam,imin,ival,imax,comment)
C--   subroutine sqcIltlt(subnam,parnam,imin,ival,imax,comment)
C--   subroutine sqcDlele(subnam,parnam,dmi,val,dma,comment)
C--   subroutine sqcDlelt(subnam,parnam,dmi,val,dma,comment)
C--   subroutine sqcDltle(subnam,parnam,dmi,val,dma,comment)
C--   subroutine sqcDltlt(subnam,parnam,dmi,val,dma,comment)

C==   ==================================================================
C==   ix,iq and x,q range checks
C==   ==================================================================

C==   ==================================================================
C==   Check cuts =======================================================
C==   ==================================================================


C     ===================================================
      integer function iqcChkLmij(subnam,w,id,ix,iq,ichk)
C     ===================================================

C--   Check if pdf(id) in workspace w exists and is non-empmty,
C--   and if (ix,iq) is within point5 limits
C--
C--   subnam (in)  : subroutine name for error message
C--   w      (in)  : store with pdf tables
C--   id     (in)  : pdf identifier in global format
C--   ix     (in)  : value of ix
C--   iq     (in)  : value of iq
C--   ichk   (in)  : 1=fatal error when ix,iq, out of range
C--
C--   iqcChkLmij = -1 : id does not exist or is empty (fatal error)
C--                 0 : OK
C--                 1 : ix > ixmax (fatal error when ichk = 1)
C--                 2 : ix < ixmin (fatal error when ichk = 1)
C--                 3 : iq < iqmin (fatal error when ichk = 1)
C--                 4 : iq > iqmax (fatal error when ichk = 1)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

      logical lqcIsFilled

      dimension w(*)

      character*(*) subnam
      character*20  etxt
      character*80  emsg

C--   Default
      iqcChkLmij = 0
C--   Check type-5 table
      if(int(iqcGetLocalId(id)/100).ne.5)  iqcChkLmij = -1
C--   Check table filled
      if(.not.lqcIsFilled(w,id))           iqcChkLmij = -1
C--   Pdf does not exist or is empty
      if(iqcChkLmij.eq.-1) then
        call smb_itoch(id,etxt,ltxt)
        write(emsg,
     +      '(''Pdf id = '',A,'' does not exist or is empty'')')
     +          etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif
C--   Check cuts
      if(ix.lt.ixmic5)      iqcChkLmij = 1
      if(ix.gt.nyy2(0))     iqcChkLmij = 2
      if(abs(iq).lt.itmic5) iqcChkLmij = 3
      if(abs(iq).gt.itmac5) iqcChkLmij = 4
C--   Yes/no fatal error, thats the question
      if(ichk.eq.0 .or. iqcChkLmij.eq.0) then
        return
      else
        call sqcIlele(subnam,'IX',
     +                ixmic5,  ix  ,nyy2(0),'IX outside grid or cuts')
        call sqcIlele(subnam,'IQ',
     +                itmic5,abs(iq),itmac5,'IQ outside grid or cuts')
      endif

      return
      end

C     ============================================
      integer function iqcIxInside(subnam,ix,ichk)
C     ============================================

C--   If inside grid or cuts then iy = iqcIxInside(subnam,ix,ichk)
C--                          else iy = -1
C--
C--   ix   (in) : x-grid point
C--   ichk (in) : 0 = return zero if outside; 1 = fatal error if outside

      implicit double precision(a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'

      character*(*) subnam

      iy          = nyy2(0)+1-ix
      iqcIxInside = iy

C--   First for the check
      if(iy.ge.0 .and. iy.le.iymac5) return
C--   Outside grid or cuts
      iqcIxInside = -1
      if(ichk.eq.0) return
C--   Now for the fatal error
      call sqcIlele(subnam,
     +             'IX',ixmic5,ix,nyy2(0),'IX outside grid or cuts')

      return
      end

C     ============================================
      integer function iqcIqInside(subnam,iq,ichk)
C     ============================================

C--   If inside grid or cuts then it = iqcIqInside(subnam,iq,ichk)
C--                          else it = 0
C--
C--   iq   (in) : q-grid point
C--   ichk (in) : 0 = return zero if outside; 1 = fatal error if outside

      implicit double precision(a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'

      character*(*) subnam

      iqcIqInside = iq
      it          = abs(iq)

C--   First for the check
      if(it.ge.itmic5 .and. it.le.itmac5) return
C--   Outside grid or cuts
      iqcIqInside = 0
      if(ichk.eq.0) return
C--   Now for the fatal error
      call sqcIlele(subnam,
     +             'IQ',itmic5,it,itmac5,'IQ outside grid or cuts')

      return
      end

C     ====================================================
      double precision function dqcXInside(subnam,xx,ichk)
C     ====================================================

C--   If inside grid or cuts then yy = dqcXInside(subnam,xx,ichk)
C--                          else yy = -1.D0
C--
C--   xx   (in) : x-value
C--   ichk (in) : 0 = return zero if outside; 1 = fatal error if outside

      implicit double precision(a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

      character*(*) subnam
      logical lmb_le, lmb_lt, lmb_eq

C--   First for the check
      if(lmb_le(xmic5,xx,-aepsi6) .and. lmb_lt(xx,1.D0,-aepsi6)) then
        dqcXInside = -log(xx)
        return
      elseif(lmb_eq(xx,1.D0,-aepsi6)) then
        dqcXInside = 0.D0
        return
      endif
C--   Outside grid or cuts
      dqcXInside = -1.D0
      if(ichk.eq.0) return
C--   Now for the fatal error
      call sqcDlele(subnam,
     +             'X',xmic5,xx,1.D0,'X outside grid or cuts')

      return
      end

C     ====================================================
      double precision function dqcQInside(subnam,qq,ichk)
C     ====================================================

C--   If inside grid or cuts then tt = dqcQInside(subnam,qq,ichk)
C--                          else tt = 0.D0
C--
C--   qq   (in) : q-value
C--   ichk (in) : 0 = return zero if outside; 1 = fatal error if outside

      implicit double precision(a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

      character*(*) subnam
      logical lmb_le

C--   First for the check
      if(lmb_le(qmic5,qq,-aepsi6) .and. lmb_le(qq,qmac5,-aepsi6)) then
        dqcQInside = log(qq)
        return
      endif
C--   Outside grid or cuts
      dqcQInside = 0.D0
      if(ichk.eq.0) return
C--   Now for the fatal error
      call sqcDlele(subnam,
     +             'Q2',qmic5,qq,qmac5,'Q2 outside grid or cuts')

      return
      end

C==   ==================================================================
C==   Range checks and associated messages  ============================
C==   ==================================================================

C     =========================================================
      subroutine sqcIlele(subnam,parnam,imin,ival,imax,comment)
C     =========================================================

C--   OK if imin .le. ival .le. imax

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*10  cmin, cval, cmax
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if( (imin.le.ival) .and. (ival.le.imax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      call smb_ItoCh(imin,cmin,lmin)
      call smb_ItoCh(ival,cval,lval)
      call smb_ItoCh(imax,cmax,lmax)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) parnam,' = ',cval(1:lval),' not in range [ ',
     +                 cmin(1:lmin),' , ',cmax(1:lmax),' ]'
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     =========================================================
      subroutine sqcIlelt(subnam,parnam,imin,ival,imax,comment)
C     =========================================================

C--   OK if imin .le. ival .lt. imax

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*10  cmin, cval, cmax
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if( (imin.le.ival) .and. (ival.lt.imax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      call smb_ItoCh(imin,cmin,lmin)
      call smb_ItoCh(ival,cval,lval)
      call smb_ItoCh(imax,cmax,lmax)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) parnam,' = ',cval(1:lval),' not in range [ ',
     +                 cmin(1:lmin),' , ',cmax(1:lmax),' )'
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     =========================================================
      subroutine sqcIltle(subnam,parnam,imin,ival,imax,comment)
C     =========================================================

C--   OK if imin .lt. ival .le. imax

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*10  cmin, cval, cmax
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if( (imin.lt.ival) .and. (ival.le.imax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      call smb_ItoCh(imin,cmin,lmin)
      call smb_ItoCh(ival,cval,lval)
      call smb_ItoCh(imax,cmax,lmax)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) parnam,' = ',cval(1:lval),' not in range ( ',
     +                 cmin(1:lmin),' , ',cmax(1:lmax),' ]'
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     =========================================================
      subroutine sqcIltlt(subnam,parnam,imin,ival,imax,comment)
C     =========================================================

C--   OK if imin .lt. ival .lt. imax

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*10  cmin, cval, cmax
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if( (imin.lt.ival) .and. (ival.lt.imax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      call smb_ItoCh(imin,cmin,lmin)
      call smb_ItoCh(ival,cval,lval)
      call smb_ItoCh(imax,cmax,lmax)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) parnam,' = ',cval(1:lval),' not in range ( ',
     +                 cmin(1:lmin),' , ',cmax(1:lmax),' )'
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     ======================================================
      subroutine sqcDlele(subnam,parnam,dmi,val,dma,comment)
C     ======================================================

C--   OK if dmi .le. val .le. dma
C--   We use lmb_eq to check equality within relative tolerance aepsi6

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'

      logical       lmb_le
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if(lmb_le(dmi,val,-aepsi6) .and. lmb_le(val,dma,-aepsi6)) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,'( 1X,A,'' = '',G11.4,'' not in range [ '',G11.4,
     +                  '' , '',G11.4,'' ]'')') parnam,val,dmi,dma
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     ======================================================
      subroutine sqcDlelt(subnam,parnam,dmi,val,dma,comment)
C     ======================================================

C--   OK if dmi .le. val .lt. dma
C--   We use lmb_eq to check equality within relative tolerance aepsi6

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'

      logical       lmb_le, lmb_lt
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if(lmb_le(dmi,val,-aepsi6) .and. lmb_lt(val,dma,-aepsi6)) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,'( 1X,A,'' = '',G11.4,'' not in range [ '',G11.4,
     +                  '' , '',G11.4,'' )'')') parnam,val,dmi,dma
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     ======================================================
      subroutine sqcDltle(subnam,parnam,dmi,val,dma,comment)
C     ======================================================

C--   OK if dmi .lt. val .le. dma
C--   We use lmb_eq to check equality within relative tolerance aepsi6

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'

      logical       lmb_lt, lmb_le
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if(lmb_lt(dmi,val,-aepsi6) .and. lmb_le(val,dma,-aepsi6)) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,'( 1X,A,'' = '',G11.4,'' not in range ( '',G11.4,
     +                  '' , '',G11.4,'' ]'')') parnam,val,dmi,dma
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     ======================================================
      subroutine sqcDltlt(subnam,parnam,dmi,val,dma,comment)
C     ======================================================

C--   OK if dmi .lt. val .lt. dma
C--   We use lmb_eq to check equality within relative tolerance aepsi6

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'

      logical       lmb_lt
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if(lmb_lt(dmi,val,-aepsi6) .and. lmb_lt(val,dma,-aepsi6)) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,'( 1X,A,'' = '',G12.4,'' not in range ( '',G12.4,
     +                  '' , '',G12.4,'' )'')') parnam,val,dmi,dma
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

