
C--   file srcTboxWeights.f containing toolbox weight routines
C--
C--   subroutine sqcSetKey(keyin,keyout)
C--   logical function lqcSjekey(key1,key2)
C--   subroutine sqcDumpTab(w,kset,lun,key,ierr)
C--   subroutine sqcReadTab(w,nw,lun,key,new,kset,lastwd,ierr)

C--   subroutine sqcUweitA(w,id,afun,achi,ierr)
C--   subroutine sqcUweitB(w,id,bfun,achi,idel,ierr)
C--   subroutine sqcUwgtRS(w,id,rfun,sfun,achi,idel,ierr)
C--   subroutine sqcUweitD(w,id,dfun,achi,ierr)
C--   subroutine sqcUweitX(w,id,ierr)

C--   subroutine sqcScaleWt(w,c,id)
C--   subroutine sqcCopyWt(w1,id1,w2,id2,iadd)
C--   subroutine sqcWcrossW(wa,ida,wb,idb,wc,idc,igt1,igt2,iadd)
C--   subroutine sqcWtimesF(fun,w1,id1,w2,id2,iadd)

C--   double precision function dqcUIgauss(pfun,ti,nf,a,b)
C--   double precision function 
C--  +   dqcUAgauss(idk,afun,yi,ti,nf,a,b,del)
C--   double precision function 
C--  +   dqcUBgauss(idk,bfun,yi,ti,nf,a,b,del)
C--   double precision function 
C--  +   dqcURSgaus(idk,rfun,sfun,yi,ti,nf,a,b,del)

C     ==================================
      subroutine sqcSetKey(keyin,keyout)
C     ==================================

C--   Left adjust, truncate to 50 chars and convert to upper case
C--
C--   ierr   (out)  0 = OK
C--                 1 = key truncated

      implicit double precision (a-h,o-z)
      
      character*(*) keyin
      character*50  keyout
      
      call smb_cfill(' ',keyout)
      i1   = imb_frstc(keyin)
      i2   = imb_lenoc(keyin)
      if(i1.eq.i2) return           !empty string
      leng = min(i2-i1+1,50)
      keyout = keyin(i1:i1+leng-1)
      call smb_cltou(keyout)
      
      return
      end
      
C     =====================================      
      logical function lqcSjekey(key1,key2)
C     =====================================

C--   True if two keys match

      implicit double precision(a-h,o-z)
      
      character*(*) key1,key2
      character*50  kkk1,kkk2              

      call sqcSetKey(key1,kkk1)
      call sqcsetKey(key2,kkk2)
      
      lqcSjekey = .false.
      if(kkk1.eq.kkk2) lqcSjekey = .true.
      
      return
      end
      
C     ==========================================
      subroutine sqcDumpTab(w,kset,lun,key,ierr)
C     ==========================================

C--   Dump set of tables on logical unit number lun
C--
C--   ierr = 0  : all OK
C--          1  : write error

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qvers1.inc'
      include 'qgrid2.inc'

      dimension w(*)
      character*(*) key
      character*50  keyout

      dimension itypes(mtyp0)

C--   Initialize
      ierr = 0
      call sqcSetKey(key,keyout)      
C--   Dump QCDNUM version
      write(lun,err=500) cvers1, cdate1
C--   Dump key
      write(lun,err=500) keyout
C--   Dump some array sizes
      write(lun,err=500) mxg0, mxx0, mqq0, mst0
      write(lun,err=500) mord0, mnf0, mbp0, mpp0, maa0, mtyp0, mchk0
C--   Dump relevant grid parameters
      write(lun,err=500) nyy2, nyg2, ioy2, dely2
      write(lun,err=500) ntt2
      write(lun,err=500) (tgrid2(i),i=1,ntt2)
C--   Fill itypes array
      do ityp = 1,mtyp0
        itypes(ityp) = iqcSgnNumberOfTables(w,kset,ityp)
      enddo
C--   Number of pars, words and position of kset
      npar   = iqcGetNumberOfParams(w,kset)     !number of evol params
      nusr   = iqcGetNumberOfUParam(w,kset)     !number of user params
      ifirst = iqcFirstWordOfSet(w,kset)        !first word of kset
      nwdset = iqcGetNumberOfWords(w(ifirst))   !number of words in kset
      ilast  = ifirst + nwdset - 1              !last word of set
C--   Control word
      icword = 123456
C--   Go...
      write(lun,err=500) icword
      write(lun,err=500) nwdset,itypes,npar,nusr
      write(lun,err=500) (w(i),i=ifirst,ilast)

      return

 500  continue
C--   Write error
      ierr = 1
      return

      end

C     ========================================================
      subroutine sqcReadTab(w,nw,lun,key,new,kset,lastwd,ierr)
C     ========================================================

C--   Read table set from logical unit number lun
C--
C--   w            (in)   store
C--   nw           (in)   number of words in the store
C--   lun          (in)   logical unit number
C--   key          (in)   key character string
C--   new          (in)   0=add new set 1=overwrite existing sets
C--   kset         (out)  table set number in w
C--   lastwd       (out)  last word used in the store < 0 no space
C--   ierr         (out)  0 = all OK
C--                       1 = read error
C--                       2 = problem with QCDNUM version
C--                       3 = key mismatch
C--                       4 = x-mu2 grid not the same
C--                       5 = not enough space
C--                       6 = iset count exceeded [1,mst0]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qvers1.inc'
      include 'qgrid2.inc'

      character*10  cversr
      character*8   cdater
      character*50  keyred
      character*(*) key
      logical lqcSjekey
      dimension nyyr(0:mxg0),delyr(0:mxg0)
      dimension tgridr(mqq0)
      dimension itypes(mtyp0)

      dimension w(*)

C--   Initialize
      nwords = 0
      ierr   = 0
C--   QCDNUM version
      read(lun,err=500,end=500) cversr, cdater
      if(cversr.ne.cvers1 .or. cdater.ne.cdate1) then
        ierr = 2
        return
      endif
C--   Key
      read(lun,err=500,end=500) keyred
      if(.not.lqcSjekey(key,keyred)) then
        ierr = 3
        return
      endif
C--   Some array sizes
      read(lun,err=500,end=500) mxgr, mxxr, mqqr, mstr
      if(mxgr.ne.mxg0 .or. mxxr.ne.mxx0 .or.
     +   mqqr.ne.mqq0 .or. mstr.ne.mst0) then
        ierr = 2
        return
      endif
C--   More array sizes
      read(lun,err=500,end=500)
     +       mordr, mnfr, mbpr, mppr, maar, mtypr, mchkr
      if(mordr.ne.mord0 .or. mnfr .ne.mnf0 .or. mbpr.ne.mbp0 .or.
     +   mppr .ne.mpp0  .or. maar .ne.maa0 .or.
     +   mtypr.ne.mtyp0 .or. mchkr.ne.mchk0      ) then
        ierr = 2
        return
      endif
C--   Relevant x-grid parameters
      read(lun,err=500,end=500) nyyr, nygr, ioyr, delyr
      if(nygr.ne.nyg2 .or. ioyr.ne.ioy2) then
        ierr = 4
        return
      endif
      do i = 0,mxg0
        if(nyyr(i).ne.nyy2(i) .or. delyr(i).ne.dely2(i)) then
          ierr = 4
          return
        endif
      enddo
C--   Mu2 grid
      read(lun,err=500,end=500) nttr
      if(nttr .ne. ntt2) then
        ierr = 4
        return
      endif
      read(lun,err=500,end=500) (tgridr(i),i=1,ntt2)
      do i = 1,ntt2
        if(tgrid2(i).ne.tgridr(i)) then
          ierr = 4
          return
        endif
      enddo
C--   Read header information
      read(lun,err=500,end=500) icword
      if(icword.ne.123456) goto 500
      read(lun,err=500,end=500) nwords,itypes,npar,nusr
C--   Book weight tables in the store
      call sqcMakeTab(w,nw,itypes,npar,nusr,new,kset,lastwd)
      if(kset.gt.0) then
        ierr = 0
      elseif(kset.eq.-1) then
        stop 'sqcReadTab empty set encountered'
      elseif(kset.eq.-2) then
        ierr = 5
        return
      elseif(kset.eq.-3) then
        ierr = 6
        return
      else
        stop 'sqcReadTab unknown error from sqcMakeTab'
      endif
C--   Number of words and position of kset
      ifirst = iqcFirstWordOfSet(w,kset)        !first word of kset
      nwdset = iqcGetNumberOfWords(w(ifirst))   !number of words in kset
      ilast  = ifirst + nwdset - 1              !last word of set
C--   Check
      if(nwdset.ne.nwords) goto 500
C--   Now read the store
      read(lun,err=500,end=500) (w(i),i=ifirst,ilast)

      return

  500 continue
C--   Read error
      ierr = 1
      return

      end

C==   ==============================================================
C==   Weight calculation ===========================================
C==   ==============================================================

C     =========================================
      subroutine sqcUweitA(w,id,afun,achi,ierr)
C     =========================================

C--   Calculate weights for regular piece (id in global format)

C--   ierr = 0  : OK
C--          1  : achi(qmu2) < 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(6),ma(6)

      external afun, achi

      logical lmb_eq

      kset = id/1000
      ityp = (id-1000*kset)/100
      ierr = 0
     
C--   Only table types 1,2,3,4, thank you
      if(ityp.lt.1 .or. ityp.gt.4) stop 'sqcUweitA: invalid table type'

C--   Get index limits      
      call sqcGetLimits(w,id,mi,ma,jmax)

C--   Loop over ioy2 range
      iorem = ioy2

      do ioy2 = mi(6),ma(6)

C--     Validate table
        call sqcValidate(w,id)

C--     Prepare fast indexing
        inc1 = iqcGaddr(w,1,0,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc2 = iqcGaddr(w,0,1,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc3 = iqcGaddr(w,0,0,1,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc4 = iqcGaddr(w,0,0,0,1,id)-iqcGaddr(w,0,0,0,0,id)
        ia4  = iqcGaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--     Nested loop: depending on the table type, several of these are
C--     one-trip loops. The index ranges are set by sqcGetLimits above
C--     Loop over subgrids
        do ig = mi(4),ma(4)
          ia4 = ia4+inc4
          ia3 = ia4-inc3
          del = dely2(ig)
C--       Loop over number of flavors
          do nf = mi(3),ma(3)
            ia3 = ia3+inc3
            ia2 = ia3-inc2
C--         Loop over t-grid
            do it = mi(2),ma(2)
              ia2 = ia2+inc2
              ia1 = ia2-inc1
              ti  = tgrid2(it)
              aa  = achi(exp(ti))
C--           Check
              if(lmb_eq(aa,1.D0,aepsi6)) then
                aa = 1.D0
              elseif(aa.lt.1.D0)           then
                ierr = 1
                return
              endif
              bb  = log(aa)
C--           Loop over y-subgrid
              do iy = 1,nyy2(ig)
                yi   = iy*del
                yb   = yi-bb
                weit = 0.D0
                if(yb.gt.0.D0) then
                  a    = 0.D0              !lower integration limit
                  b    = min(yb,ioy2*del)  !upper integration limit
                  weit = dqcUAgauss(ioy2-1,afun,yb,ti,nf,a,b,del)
                  weit = weit/aa
                endif
C--             Add weight
                ia1    = ia1+inc1
                w(ia1) = w(ia1)+weit
C--           End of loop over y-subgrid
              enddo
C--         End of loop over t-grid
            enddo
C--       End of loop over flavors
          enddo
C--     End of loop over subgrids
        enddo

C--   End of loop over ioy2 range
      enddo

C--   Validate table
      call sqcValidate(w,id)

      ioy2 = iorem

      return
      end

C     ==============================================
      subroutine sqcUweitB(w,id,bfun,achi,idel,ierr)
C     ==============================================

C--   Calculate weights for singular piece (id in global format)

C--   idel = 0  : dont include delta(1-x) piece
C--             : otherwise include this piece
C--   ierr = 0  : OK
C--          1  : achi(qmu2) < 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(6),ma(6)

      external bfun, achi

      logical lmb_eq

      kset = id/1000
      ityp = (id-1000*kset)/100
      ierr = 0

C--   Only table types 1,2,3,4, thank you
      if(ityp.lt.1 .or. ityp.gt.4) stop 'sqcUweitB: ivalid table type'

C--   Get index limits      
      call sqcGetLimits(w,id,mi,ma,jmax)

C--   Loop over ioy2 range
      iorem = ioy2

      do ioy2 = mi(6),ma(6)

C--     Prepare fast indexing
        inc1 = iqcGaddr(w,1,0,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc2 = iqcGaddr(w,0,1,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc3 = iqcGaddr(w,0,0,1,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc4 = iqcGaddr(w,0,0,0,1,id)-iqcGaddr(w,0,0,0,0,id)
        ia4  = iqcGaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--     Nested loop: depending on the table type, several of these are
C--     one-trip loops. The index ranges are set by sqcGetLimits above
C--     Loop over subgrids
        do ig = mi(4),ma(4)
          ia4 = ia4+inc4
          ia3 = ia4-inc3
          del = dely2(ig)
C--       Loop over number of flavors
          do nf = mi(3),ma(3)
            ia3 = ia3+inc3
            ia2 = ia3-inc2
C--         Loop over t-grid
            do it = mi(2),ma(2)
              ia2 = ia2+inc2
              ia1 = ia2-inc1
              ti = tgrid2(it)
              aa  = achi(exp(ti))
C--           Check
              if(lmb_eq(aa,1.D0,aepsi6)) then
                aa = 1.D0
              elseif(aa.lt.1.D0)           then
                ierr = 1
                return
              endif
              bb  = log(aa)
C--           Loop over y-subgrid
              do iy = 1,nyy2(ig)
                yi     = iy*del
                yb     = yi-bb
                weit   = 0.D0
                if(yb.gt.0.D0) then
                  xi     = exp(-yi)
                  ax     = aa*xi
                  a      = 0.D0              !lower integration limit
                  b      = min(yb,ioy2*del)  !upper integration limit
                  wgt1   = dqcUBgauss(ioy2-1,bfun,yb,ti,nf,a,b,del)
                  if(idel.ne.0) then
                    wgt2   = dqcBsplyy(ioy2-1,1,yb/del) *
     +                       dqcUIgauss(bfun,ti,nf,0.D0,ax)
                  else
                    wgt2 = 0.D0
                  endif
                  weit   = (wgt1-wgt2)/aa
                endif
C--             Add weight
                ia1    = ia1+inc1
                w(ia1) = w(ia1)+weit
C--           End of loop over y-subgrid
              enddo
C--         End of loop over t-grid
            enddo
C--       End of loop over flavors
          enddo
C--     End of loop over subgrids
        enddo

C--   End of loop over ioy2 range
      enddo

C--   Validate table
      call sqcValidate(w,id)

      ioy2 = iorem

      return
      end

C     ===================================================
      subroutine sqcUwgtRS(w,id,rfun,sfun,achi,idel,ierr)
C     ===================================================

C--   Calculate weights for RS piece (id in global format)

C--   idel = 0  : dont include delta(1-x) piece
C--             : otherwise include this piece
C--   ierr = 0  : OK
C--          1  : achi(qmu2) < 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(6),ma(6)

      external rfun, sfun, achi

      logical lmb_eq

      kset = id/1000
      ityp = (id-1000*kset)/100
      ierr = 0

C--   Only table types 1,2,3,4, thank you
      if(ityp.lt.1 .or. ityp.gt.4) stop 'sqcUwgtRS: invalid table type'

C--   Get index limits      
      call sqcGetLimits(w,id,mi,ma,jmax)

C--   Loop over ioy2 range
      iorem = ioy2

      do ioy2 = mi(6),ma(6)

C--     Prepare fast indexing
        inc1 = iqcGaddr(w,1,0,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc2 = iqcGaddr(w,0,1,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc3 = iqcGaddr(w,0,0,1,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc4 = iqcGaddr(w,0,0,0,1,id)-iqcGaddr(w,0,0,0,0,id)
        ia4  = iqcGaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--     Nested loop: depending on the table type, several of these are
C--     one-trip loops. The index ranges are set by sqcGetLimits above
C--     Loop over subgrids
        do ig = mi(4),ma(4)
          ia4  = ia4+inc4
          ia3  = ia4-inc3
          del  = dely2(ig)
C--       Loop over number of flavors
          do nf = mi(3),ma(3)
            ia3 = ia3+inc3
            ia2 = ia3-inc2
C--         Loop over t-grid
            do it = mi(2),ma(2)
              ia2 = ia2+inc2
              ia1 = ia2-inc1
              ti = tgrid2(it)
              aa  = achi(exp(ti))
C--           Check
              if(lmb_eq(aa,1.D0,aepsi6)) then
                aa = 1.D0
              elseif(aa.lt.1.D0)           then
                ierr = 1
                return
              endif
              bb  = log(aa)
C--           Loop over y-subgrid
              do iy = 1,nyy2(ig)
                yi   = iy*del
                yb   = yi-bb
                weit = 0.D0
                if(yb.gt.0.D0) then
                  xi   = exp(-yi)
                  ax   = aa*xi
                  a    = 0.D0              !lower integration limit
                  b    = min(yb,ioy2*del)  !upper integration limit
                  wgt1 = dqcURSgaus(
     +                   ioy2-1,rfun,sfun,yb,ti,nf,a,b,del)
                  if(idel.ne.0) then
                    wgt2 = rfun(1.D0,ti,nf)          *
     +                     dqcBsplyy(ioy2-1,1,yb/del) *
     +                     dqcUIgauss(sfun,ti,nf,0.D0,ax)
                  else
                    wgt2 = 0.D0
                  endif
                  weit = (wgt1-wgt2)/aa
                endif
C--             Add weight
                ia1    = ia1+inc1
                w(ia1) = w(ia1)+weit
C--           End of loop over y-subgrid
              enddo
C--         End of loop over t-grid
            enddo
C--       End of loop over flavors
          enddo
C--     End of loop over subgrids
        enddo

C--   End of loop over ioy2 range
      enddo

C--   Validate table
      call sqcValidate(w,id)

      ioy2 = iorem

      return
      end

C     =========================================
      subroutine sqcUweitD(w,id,dfun,achi,ierr)
C     =========================================

C--   Calculate weights for D(x)*delta(1-x) (id in global format)
C--
C--   ierr = 0  : OK
C--          1  : achi(qmu2) < 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(6),ma(6)

      external dfun, achi

      logical lmb_eq

      kset = id/1000
      ityp = (id-1000*kset)/100
      ierr = 0

C--   Only table types 1,2,3,4, thank you
      if(ityp.lt.1 .or. ityp.gt.4) stop 'sqcUweitD: invalid table type'

C--   Get index limits
      call sqcGetLimits(w,id,mi,ma,jmax)

C--   Loop over ioy2 range
      iorem = ioy2

      do ioy2 = mi(6),ma(6)

C--     Prepare fast indexing
        inc1 = iqcGaddr(w,1,0,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc2 = iqcGaddr(w,0,1,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc3 = iqcGaddr(w,0,0,1,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc4 = iqcGaddr(w,0,0,0,1,id)-iqcGaddr(w,0,0,0,0,id)
        ia4  = iqcGaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--     Nested loop: depending on the table type, several of these are
C--     one-trip loops. The index ranges are set by sqcGetLimits above
C--     Loop over subgrids
        do ig = mi(4),ma(4)
          ia4  = ia4+inc4
          ia3  = ia4-inc3
          del  = dely2(ig)
C--       Loop over number of flavors
          do nf = mi(3),ma(3)
            ia3 = ia3+inc3
            ia2 = ia3-inc2
C--         Loop over t-grid
            do it = mi(2),ma(2)
              ia2 = ia2+inc2
              ia1 = ia2-inc1
              ti  = tgrid2(it)
              aa  = achi(exp(ti))
C--           Check
              if(lmb_eq(aa,1.D0,aepsi6)) then
                aa = 1.D0
              elseif(aa.lt.1.D0)           then
                ierr = 1
                return
              endif
              bb  = log(aa)
C--           Loop over y-subgrid
              do iy = 1,nyy2(ig)
                yi   = iy*del
                yb   = yi-bb
                wgt  = 0.D0
                if(yb.gt.0.D0) then
                  wgt  = dfun(exp(-yb),exp(ti),nf) *
     +                   dqcBsplyy(ioy2-1,1,yb/del)
                  wgt  = wgt/aa
                endif
C--             Add weight
                ia1    = ia1+inc1
                w(ia1) = w(ia1)+wgt
C--           End of loop over y-subgrid
              enddo
C--         End of loop over t-grid
            enddo
C--       End of loop over flavors
          enddo
C--     End of loop over subgrids
        enddo

C--   End of loop over ioy2 range
      enddo

C--   Validate table
      call sqcValidate(w,id)

      ioy2 = iorem

      return
      end

C     ===============================
      subroutine sqcUweitX(w,id,ierr)
C     ===============================

C--   Calculate weights for convolution F cross F (id in global format)

C--   ierr = 0  : OK
C--          1  : error (should never occur)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(6),ma(6)

      kset = id/1000
      ityp = (id-1000*kset)/100

      ierr = 0

C--   Only table types 1,2,3,4, thank you
      if(ityp.lt.1 .or. ityp.gt.4) stop 'sqcUweitX: invalid table type'

C--   Get index limits      
      call sqcGetLimits(w,id,mi,ma,jmax)

C--   Loop over ioy2 range
      iorem = ioy2

      do ioy2 = mi(6),ma(6)

C--     Prepare fast indexing
        inc1 = iqcGaddr(w,1,0,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc2 = iqcGaddr(w,0,1,0,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc3 = iqcGaddr(w,0,0,1,0,id)-iqcGaddr(w,0,0,0,0,id)
        inc4 = iqcGaddr(w,0,0,0,1,id)-iqcGaddr(w,0,0,0,0,id)
        ia4  = iqcGaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--     Nested loop: depending on the table type, several of these are
C--     one-trip loops. The index ranges are set by sqcGetLimits above
C--     Loop over subgrids
        do ig = mi(4),ma(4)
          ia4 = ia4+inc4
          ia3 = ia4-inc3
          del = dely2(ig)
C--       Loop over number of flavors
          do nf = mi(3),ma(3)
            ia3 = ia3+inc3
            ia2 = ia3-inc2
C--         Loop over t-grid
            do it = mi(2),ma(2)
              ia2 = ia2+inc2
              ia1 = ia2-inc1
              ti  = tgrid2(it)
C--           Loop over y-subgrid
              do iy = 1,nyy2(ig)
                yi   = iy*del
                weit = 0.D0
                a    = 0.D0              !lower integration limit
                b    = min(yi,ioy2*del)  !upper integration limit
                weit = dqcUXgauss(ioy2-1,yi,a,b,del)
C--             Store weight
                ia1    = ia1+inc1
                w(ia1) = weit
C--           End of loop over y-subgrid
              enddo
C--         End of loop over t-grid
            enddo
C--       End of loop over flavors
          enddo
C--     End of loop over subgrids
        enddo

C--   End of loop over ioy2 range
      enddo

C--   Validate table
      call sqcValidate(w,id)

      ioy2 = iorem

      return
      end

C==   ==============================================================
C==   Operations on weight tables ==================================
C==   ==============================================================      

C     =============================
      subroutine sqcScaleWt(w,c,id)
C     =============================

C--   Multiply weight table (id in global format) by constant factor c

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(6),ma(6)

C--   Get index limits
      call sqcGetLimits(w,id,mi,ma,jmax)

C--   Loop over iosp range
      iorem = ioy2

      do ioy2 = mi(6),ma(6)

C--     First and last word of table
        ia1 = iqcGaddr(w,mi(1),mi(2),mi(3),mi(4),id)
        ia2 = iqcGaddr(w,ma(1),ma(2),ma(3),ma(4),id)
C--     Multiply with scalefactor
        do ia = ia1,ia2
          w(ia) = c*w(ia)
        enddo

C--   End of loop over iosp range
      enddo

C--   Validate the table
      call sqcValidate(w,id)

      ioy2 = iorem

      return
      end
      
C     ========================================
      subroutine sqcCopyWt(w1,id1,w2,id2,iadd)
C     ========================================

C--   Copy table id1 in w1 to table id2 in w2 (id in global format)
C--   iadd = -1,0,1    subtract/copy/add  id1 to id2
C--   Type of id2 should be .ge. type id1 but this is not checked here

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w1(*),w2(*)
      dimension mi1(6),ma1(6),mi2(6),ma2(6)

C--   Get index limits
      call sqcGetLimits(w1,id1,mi1,ma1,jmax1)
      call sqcGetLimits(w2,id2,mi2,ma2,jmax2)

C--   Loop over iosp range
      iorem = ioy2
      io1   = max(mi1(6),mi2(6))
      io2   = min(ma1(6),ma2(6))

      do ioy2 = io1,io2

C--     Loop over output table
        do ig2 = mi2(4),ma2(4)              !loop over grids
C--       Bracket ig1 limits
          ig1 = max(mi1(4),ig2)
          ig1 = min(ma1(4),ig1)
          do nf2 = mi2(3),ma2(3)            !loop over flavors
C--         Bracket nf1  limits
            nf1 = max(mi1(3),nf2)
            nf1 = min(ma1(3),nf1)
            do iq2 = mi2(2),ma2(2)          !loop over mu2
C--           Bracket iq1 limits
              iq1 = max(mi1(2),iq2)
              iq1 = min(ma1(2),iq1)
C--           Base address of first x-bin
              ia1 = iqcGaddr(w1,mi1(1),iq1,nf1,ig1,id1)-1
              ia2 = iqcGaddr(w2,mi2(1),iq2,nf2,ig2,id2)-1
              if(iadd.eq.-1) then
C--             id2 = id2 - id1
                do i = mi2(1),ma2(1)
                  ia1     = ia1+1
                  ia2     = ia2+1
                  w2(ia2) = w2(ia2) - w1(ia1)
                enddo
              elseif(iadd.eq.0) then
C--             id2 = id1
                do i = mi2(1),ma2(1)
                  ia1    = ia1+1
                  ia2    = ia2+1
                  w2(ia2) = w1(ia1)
                enddo
              elseif(iadd.eq.+1) then
C--             id2 = id2 + id1
                do i = mi2(1),ma2(1)
                  ia1 = ia1+1
                  ia2 = ia2+1
                  w2(ia2) = w2(ia2) + w1(ia1)
                enddo
              else
                stop 'sqcCopyWt: invalid iadd'
              endif
            enddo
          enddo
        enddo

C--   End of loop over iosp range
      enddo

C--   Copy entries of satellite table
      ia1 = iqcGSij(w1,1,id1)-1
      ia2 = iqcGSij(w2,1,id2)-1
      do j = 1,jmax1
        ia1 = ia1+1
        ia2 = ia2+1
        w2(ia2) = w1(ia1)
      enddo

      ioy2 = iorem

      return
      end
      
C     ==========================================================
      subroutine sqcWcrossW(wa,ida,wb,idb,wc,idc,igt1,igt2,iadd)
C     ==========================================================

C--   Store in idc the weight of ida * idb :   Wc = Wa S^{-1} Wb
C--
C--   wa     (in) :  workspace containing Wa
C--   ida    (in) :  global id of Wa
C--   wb     (in) :  workspace containing Wb
C--   idb    (in) :  global id of Wb
C--   wc     (in) :  workspace containing Wc
C--   idc    (in) :  global id of Wc
C--   igt1,2 (in) :  global ids of scratch pdf tables in stor7
C--   iadd   (in) : -1/0/1 subtract/put/add result to Wc

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension wa(*), wb(*), wc(*)
      dimension mia(6),maa(6),mib(6),mab(6),mic(6),mac(6)
      
C--   Get index limits
      call sqcGetLimits(wa,ida,mia,maa,jmaxa)
      call sqcGetLimits(wb,idb,mib,mab,jmaxb)
      call sqcGetLimits(wc,idc,mic,mac,jmaxc)

C--   Loop over iosp range
      iorem = ioy2
      io1   = max(mia(6),mib(6),mic(6))
      io2   = min(maa(6),mab(6),mac(6))

      do ioy2 = io1,io2

C--     Loop over output table
        do igc = mic(4),mac(4)              !loop over grids
C--       Bracket iga and igb limits
          iga = max(mia(4),igc)
          iga = min(maa(4),iga)
          igb = max(mib(4),igc)
          igb = min(mab(4),igb)
          do nfc = mic(3),mac(3)            !loop over flavors
C--         Bracket nfa and nfb limits
            nfa = max(mia(3),nfc)
            nfa = min(maa(3),nfa)
            nfb = max(mib(3),nfc)
            nfb = min(mab(3),nfb)
            do iqc = mic(2),mac(2)          !loop over mu2
C--           Bracket iqa and iqb limits
              iqa = max(mia(2),iqc)
              iqa = min(maa(2),iqa)
              iqb = max(mib(2),iqc)
              iqb = min(mab(2),iqb)
C--           Address of first x-bin
              ia1 = iqcGaddr(wa,mia(1),iqa,nfa,iga,ida)
              ib1 = iqcGaddr(wb,mib(1),iqb,nfb,igb,idb)
              ic1 = iqcGaddr(wc,mic(1),iqc,nfc,igc,idc)
              it1 = iqcG5ijk(stor7,1,1,igt1)
              it2 = iqcG5ijk(stor7,1,1,igt2)
C--           First do Z = S^{-1}B   (store in it1)
              call sqcABmult(sinvy2(1,ioy2),wb(ib1),stor7(it1),nyy2(0))
C--           Now do C = AZ (store in it2)
              call sqcABmult(wa(ia1),stor7(it1),stor7(it2),nyy2(0))
              ic1 = ic1-1
              it2 = it2-1
              if(iadd.eq.-1) then
C--             C = C - convol
                do i = mic(1),mac(1)
                  ic1     = ic1+1
                  it2     = it2+1
                  wc(ic1) = wc(ic1) - stor7(it2)
                enddo
              elseif(iadd.eq.0) then
C--             C = convol
                do i = mic(1),mac(1)
                  ic1     = ic1+1
                  it2     = it2+1
                  wc(ic1) = stor7(it2)
                enddo
              elseif(iadd.eq.+1) then
C--             C = C + convol
                do i = mic(1),mac(1)
                  ic1     = ic1+1
                  it2     = it2+1
                  wc(ic1) = wc(ic1) + stor7(it2)
                enddo
              else
                stop 'sqcWcrossW: invalid iadd'
              endif
            enddo
          enddo
        enddo

C--   End op loop over iosp range
      enddo

C--   Validate output table
      call sqcValidate(wc,idc)

      ioy2 = iorem

      return
      end

C     =============================================
      subroutine sqcWtimesF(fun,w1,id1,w2,id2,iadd)
C     =============================================

C--   Multiply id1 by fun(iq,nf) and store result in id2
C--   iadd = -1,0,1   subtract,store,add  result to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension w1(*),w2(*)
      dimension mi1(6),ma1(6),mi2(6),ma2(6)
      
      external fun

C--   Get index limits
      call sqcGetLimits(w1,id1,mi1,ma1,jmax1)
      call sqcGetLimits(w2,id2,mi2,ma2,jmax2)

C--   Loop over iosp range
      iorem = ioy2
      io1   = max(mi1(6),mi2(6))
      io2   = min(ma1(6),ma2(6))

      do ioy2 = io1,io2

C--     Loop over output table
        do ig2 = mi2(4),ma2(4)              !loop over grids
C--       Bracket ig1 limits
          ig1 = max(mi1(4),ig2)
          ig1 = min(ma1(4),ig1)
          do nf2 = mi2(3),ma2(3)            !loop over flavors
C--         Bracket nf1  limits
            nf1 = max(mi1(3),nf2)
            nf1 = min(ma1(3),nf1)
            do iq2 = mi2(2),ma2(2)          !loop over mu2
C--           Bracket iq1 limits
              iq1 = max(mi1(2),iq2)
              iq1 = min(ma1(2),iq1)
C--           Base address of 1st x-bin
              ia1 = iqcGaddr(w1,mi1(1),iq1,nf1,ig1,id1)-1
              ia2 = iqcGaddr(w2,mi2(1),iq2,nf2,ig2,id2)-1
              fac = fun(iq2,nf2)
              if(iadd.eq.-1) then
C--             id2 = id2 - f*id1
                do i = mi2(1),ma2(1)
                  ia1     = ia1+1
                  ia2     = ia2+1
                  w2(ia2) = w2(ia2) - fac*w1(ia1)
                enddo
              elseif(iadd.eq.0) then
C--             id2 = f*id1
                do i = mi2(1),ma2(1)
                  ia1     = ia1+1
                  ia2     = ia2+1
                  w2(ia2) = fac*w1(ia1)
                enddo
              elseif(iadd.eq.+1) then
C--             id2 = id2 + f*id1
                do i = mi2(1),ma2(1)
                  ia1     = ia1+1
                  ia2     = ia2+1
                  w2(ia2) = w2(ia2) + fac*w1(ia1)
                enddo
              else
                stop 'sqcWtimesF: invalid iadd'
              endif
            enddo
          enddo
        enddo

C--   End op loop over iosp range
      enddo

C--   Validate output table
      call sqcValidate(w2,id2)

      ioy2 = iorem

      return
      end

C==   ==============================================================
C==   Gauss integration ============================================
C==   ==============================================================

C     ====================================================
      double precision function dqcUIgauss(pfun,ti,nf,a,b)
C     ====================================================
 
C--   Integrate coeffient function pfun(x,q,nf) from a to b.
C--   NB: integral over x and not over y = -ln x!
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      external pfun
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
      f(u) = pfun(u,exp(ti),nf)
      
      eps  = gepsi6
C--   End of MB extension
 
      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
      S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3 CONTINUE
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
      S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4 CONTINUE
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcUIgauss: too high accuracy required'','//
     +         '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcUIgauss=H

      RETURN
      END

C     ========================================
      double precision function 
     +   dqcUAgauss(idk,afun,yi,ti,nf,a,b,del)
C     ========================================

C--   Integrate regular piece (A) of coeffient function.
C--
C--   I(yi) = Int_a^b dz Abar(yi-z) B1(z) with Abar(z) = exp(-z)A(exp(-z)).
C--
C--   idk        = spline order iord-1: 1=linear and 2=quadratic
C--   afun       = coefficient function versus x = exp(-y)
C--   yi         = upper limit of convolution (exp(-yi) is argument of afun)
C--   ti         = t-value (exp(t) is argument of afun)      
C--   nf         = number of flavors (argument of afun) 
C--   a          = lower limit of integral = 0.D0
C--   b          = upper limit of integral = min[yi,iord*del]
C--   del        = grid spacing
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      external afun
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
C--   Define function to integrate
      f(u) = dqcBsplyy(idk,1,u/del) * exp(-(yi-u)) * 
     +       afun(exp(-(yi-u)),exp(ti),nf)
C--   Check limits
      if(b.le.a) then
        dqcUAGauss = 0.D0
        return
      endif
      eps  = gepsi6
C--   End of MB extension

      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
      S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3 CONTINUE
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
      S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4 CONTINUE
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcUAgauss: too high accuracy required'','//
     +               '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcUAgauss=H

      RETURN
      END

C     ========================================
      double precision function 
     +   dqcUBgauss(idk,bfun,yi,ti,nf,a,b,del)
C     ========================================
 
C--   Integrate singular piece (B) of splitting or coeffient function.
C--
C--   I(yi) = Int_a^b dz Bbar(yi-z)[B1(z)-B1(yi)]
C--   with Bbar(z) = exp(-z)B(exp(-z)).
C--
C--   idk        = spline order iord-1: 1=linear and 2=quadratic
C--   bfun       = coefficient function versus x = exp(-y)
C--   yi         = upper limit of convolution (exp(-yi) is argument of bfun)
C--   ti         = t-value (exp(t) is argument of bfun)      
C--   nf         = number of flavors (argument of bfun) 
C--   a          = lower limit of integral = 0.D0
C--   b          = upper limit of integral = min[yi,iord*del]
C--   del        = grid spacing
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      external bfun
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
C--   Define function to integrate
      f(u) = (dqcBsplyy(idk,1,u/del) - dqcBsplyy(idk,1,yi/del)) * 
     +  exp(-(yi-u)) * bfun(exp(-(yi-u)),exp(ti),nf)
C--   Check limits
      if(b.le.a) then
        dqcUBGauss = 0.D0
        return
      endif
      eps  = gepsi6
C--   End of MB extension
 
      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
      S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3 CONTINUE
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
      S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4 CONTINUE
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcUBgauss: too high accuracy required'','//
     +         '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcUBgauss=H

      RETURN
      END

C     =============================================
      double precision function 
     +   dqcURSgaus(idk,rfun,sfun,yi,ti,nf,a,b,del)
C     =============================================
 
C--   Integrate singular*regular piece (R*S) of splitting function.
C--
C--   I(yi) = Int_a^b dz Sbar(yi-z)[Rbar(yi-z)B1(z)-Rbar(0)B1(y)]
C--   with Sbar(z) = exp(-z)S(exp(-z)) and Rbar(z) = R(exp(-z)).
C--
C--   idk  = spline order iord-1: 1=linear and 2=quadratic
C--   rfun = coefficient function versus x = exp(-y)
C--   sfun = coefficient function versus x = exp(-y)
C--   yi   = upper limit of convolution (exp(-yi) is argument of rfun)
C--   ti   = t-value (exp(t) is argument of rfun)      
C--   nf   = number of flavors (argument of rfun) 
C--   a    = lower limit of integral = 0.D0
C--   b    = upper limit of integral = min[yi,iord*del]
C--   del  = grid spacing
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      external rfun,sfun
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
C--   Define function to integrate
      f(u) = 
     +  (rfun(exp(-(yi-u)),exp(ti),nf)*dqcBsplyy(idk,1,u/del) -
     +   rfun(1.D0,exp(ti),nf)*dqcBsplyy(idk,1,yi/del)) * 
     +   exp(-(yi-u)) * sfun(exp(-(yi-u)),exp(ti),nf)
C--   Check limits
      if(b.le.a) then
        dqcURSgaus = 0.D0
        return
      endif
      eps  = gepsi6
C--   End of MB extension
 
      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
      S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3 CONTINUE
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
      S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4 CONTINUE
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcURSgaus: too high accuracy required'','//
     +         '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcURSgaus=H

      RETURN
      END
      
C     ====================================================
      double precision function dqcUXgauss(idk,yi,a,b,del)
C     ====================================================

C--   Integrate product of B-spline functions for convolution FcrossF.
C--
C--   I(yi) = Int_a^b dz B1(z) B1(yi-z).
C--
C--   idk        = spline order iord-1: 1=linear and 2=quadratic
C--   yi         = upper limit of convolution       
C--   a          = lower limit of integral = 0.D0
C--   b          = upper limit of integral = min[yi,iord*del]
C--   del        = grid spacing
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
C--   Define function to integrate
      f(u) = dqcBsplyy(idk,1,u/del) * dqcBsplyy(idk,1,(yi-u)/del)
C--   Check limits
      if(b.le.a) then
        dqcUXGauss = 0.D0
        return
      endif
      eps  = gepsi6
C--   End of MB extension

      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
      S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3 CONTINUE
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
      S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4 CONTINUE
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcUXgauss: too high accuracy required'','//
     +               '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcUXgauss=H

      RETURN
      END
       

