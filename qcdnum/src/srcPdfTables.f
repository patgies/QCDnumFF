
C--   This is the file srcPdfTables.f with pdf table routines

C--   subroutine sqcPdfBook(iset,n,idfst,noalf,nwlast,ierr)
C--
C--   subroutine sqcPreset(idg,val)
C--   subroutine sqcPSetjj(idg,it,val)
C--   subroutine sqcPCopjj(ig1,j1,ig2,j2)
C--   subroutine sqcPdfCop(idg1,idg2)
C--   subroutine sqcT1toT2(idgin,idgout,iy1,iy2,iz1,iz2)
C--
C--   subroutine sqcPdfCopy(w1,ig1,w2,ig2,coef,iadd)
C--
C--   subroutine sqcG0toGi(igg0,iggi,ig,nyg,iz)
C--   subroutine sqcGitoG0(iggi,ig,igg0)
C--   subroutine sqcGiFtoA(ig1,ig2,nyg,iz1,iz2)
C--   subroutine sqcGiAtoF(ig1,ig2,nyg,iz1,iz2)
C--   subroutine sqcAitoF0(iggi,ig,nyg,iz1,iz2,igg0)
C--
C--   subroutine sqcGiLtoQ(idg1,idg2,nyg,iz1,iz2)
C--   subroutine sqcGiQtoL(idg1,idg2,nyg,iz1,iz2)
C--
C--   subroutine sqcGetSplA(ww,id,iy,iz,ig,iyg,aout)

C=======================================================================
C==   Book pdf and alfas tables in the internal store ==================
C=======================================================================

C     =====================================================
      subroutine sqcPdfBook(iset,n,idfst,noalf,nwlast,ierr)
C     =====================================================

C--   Book n pdf tables + alfas tables (optional) in internal store
C--
C--   iset    (in)  [1,mset0] pdf set identifier
C--   n       (in)  [1,mpdf0] number of pdfs to book (usually 13+)
C--   idfst   (in)  user index of first table (e.g. 0=gluon)
C--   noalf   (in)  .ne.0 then do not book alfas tables
C--   nwlast (out)  number of words used (<0 no space)
C--   ierr   (out)  0 = OK
C--                -1 = empty set of tables (never occurs)
C--                -2 = not enough space
C--                -3 = iset count MST0 exceeded
C--                -4 = iset exist but with smaller n
C--                -5 = iset exist but with insufficient pointer tables
C--
C--  NB: if iset already exits, sqcPdfBook checks if all tables will fit

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension itypes(mtyp0)

      call smb_IFill(itypes,mtyp0,0)

      if(isetf7(iset).eq.0) then
C--     Book new set of tables
        new       = 0
        npar      = mpar0
        nusr      = 0
        itypes(5) = n                                        !pdf tables
        if(noalf.eq.0) itypes(6) = 2*mord0+1               !alfas tables
        itypes(7) = nsubt0                               !pointer tables
        call sqcMakeTab(stor7,nwf0,itypes,npar,nusr,new,kset,nwlast)
        if(kset.lt.0) then
C--       Maketab error
          ierr = kset
          return
        else
C--       Maketab went OK
          ierr         =  0
          isetf7(iset) =  kset
          ifrst7(iset) =  idfst
          ilast7(iset) =  idfst+n-1
          Lfill7(iset) = .false.
        endif
      elseif(iqcGetNumberOfTables(stor7,isetf7(iset),5).lt.n) then
C--     Iset exists but with less than n tables
        ierr = -4
      elseif(iqcGetNumberOfTables(stor7,isetf7(iset),7).lt.nsubt0) then
C--     Iset exists but with less than nsubt0 pointer tables
        ierr = -5
      else
C--     Iset exists with n or more tables
        nwlast       = iqcGetNumberOfWords(stor7)
        ierr         = 0
        ifrst7(iset) =  idfst
        ilast7(iset) =  idfst+n-1
        Lfill7(iset) = .false.
      endif

      return
      end

C=======================================================================
C==   Basic operations on the store ====================================
C=======================================================================

C     =============================
      subroutine sqcPreset(idg,val)
C     =============================

C--   Set all entries of idg to val

C--   idg  : stor7 pdf identifier in global format
C--   val  : input value

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

C--   Calculate base address
      call sqcPdfLims(idg,iy1,iy2,it1,it2,jmax)
      ia1 = iqcG5ijk(stor7,iy1,it1,idg)-1
      nwd = (iy2-iy1+1)*(it2-it1+1)
C--   Set nwd words to val
      do i = 1,nwd
        ia1        = ia1+1
        stor7(ia1) = val
      enddo

      return
      end
      
C     ================================
      subroutine sqcPSetjj(idg,it,val)
C     ================================

C--   Set all entries of id at bin it to val

C--   idg  : stor7 pdf identifier in global format
C--   it   : t index
C--   val  : input value

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

C--   Calculate base address
      call sqcPdfLims(idg,iy1,iy2,it1,it2,jmax)
      ia1 = iqcG5ijk(stor7,iy1,it,idg)-1
C--   Set value for all iy
      do i = iy1,iy2
        ia1        = ia1+1
        stor7(ia1) = val
      enddo

      return
      end      

C     ===================================
      subroutine sqcPCopjj(ig1,j1,ig2,j2)
C     ===================================

C--   Copy column j1 of table id1 to column j2 of table id2

C--   ig1,2 : stor7 pdf identifiers in global format
C--   j1,2  : t-index

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

C--   Calculate base addresses
      call sqcPdfLims(ig1,iy1,iy2,it1,it2,jmax)
      ia1 = iqcG5ijk(stor7,iy1,j1,ig1)-1
      ia2 = iqcG5ijk(stor7,iy1,j2,ig2)-1
C--   Copy column
      do i = iy1,iy2
        ia1        = ia1+1
        ia2        = ia2+1
        stor7(ia2) = stor7(ia1)
      enddo

      return
      end

C     ===============================
      subroutine sqcPdfCop(idg1,idg2)
C     ===============================

C--   Copy table ig1 to table ig2
C--   Copy also the entries in the satellite table
C--
C--   idg1,2 : stor7 pdf identifiers in global format

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      if(idg1.eq.idg2) return
C--   Calculate base addresses
      call sqcPdfLims(idg1,iy1,iy2,it1,it2,jmax)
      ia1 = iqcG5ijk(stor7,iy1,it1,idg1)-1
      ia2 = iqcG5ijk(stor7,iy1,it1,idg2)-1
      nwd = (iy2-iy1+1)*(it2-it1+1)
C--   Copy nwd words
      do i = 1,nwd
        ia1        = ia1+1
        ia2        = ia2+1
        stor7(ia2) = stor7(ia1)
      enddo
C--   Copy entries in satellite table
      ia1 = iqcGSij(stor7,1,idg1)-1
      ia2 = iqcGSij(stor7,1,idg2)-1
C--   Copy jmax words
      do i = 1,jmax
        ia1        = ia1+1
        ia2        = ia2+1
        stor7(ia2) = stor7(ia1)
      enddo

      return
      end
      
C     ==================================================
      subroutine sqcT1toT2(idgin,idgout,iy1,iy2,iz1,iz2)
C     ==================================================

C--   Copy given range of table idin to table idout
C--
C--   idgin    (in)   stor7 pdf identifier in global format
C--   idgout   (in)   stor7 pdf identifier in global format
C--   iy1,iy2  (in)   iy range to copy (iy1,2 in ascending order)
C--   iz1,iz2  (in)   iz range to copy (iz1,2 in ascending order)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      if(iy2.lt.iy1) stop 'sqcT1toT2: iy2 .lt. iy1'
      if(iz2.lt.iz1) stop 'sqcT1toT2: iz2 .lt. iz1'

      if(idgin.eq.idgout) return
C--   Calculate increments
      ia0  = iqcG5ijk(stor7,1,1,idgin)
      incy = iqcG5ijk(stor7,2,1,idgin) - ia0
      incz = iqcG5ijk(stor7,1,2,idgin) - ia0
      iain = iqcG5ijk(stor7,iy1,iz1,idgin)
      idif = iqcG5ijk(stor7,iy1,iz1,idgout)-iain
      iaz  = iain-incz
      do iz = iz1,iz2
        iaz = iaz + incz
        iay = iaz - incy
        do iy = iy1,iy2 
          iay        = iay + incy
          ia2        = iay + idif
          stor7(ia2) = stor7(iay)
        enddo  
      enddo

      return
      end

C=======================================================================
C==   More advanced copy routine
C=======================================================================

C     ==============================================
      subroutine sqcPdfCopy(w1,ig1,w2,ig2,coef,iadd)
C     ==============================================

C--   Copy coef*ig1 (type5) in w1 to ig2 (type5) in w2
C--
C--   w1        (in) input workspace
C--   ig1       (in) identifier of pdf table in w1 in global format
C--   w2        (in) target workspace
C--   ig2       (in) identifier of pdf table in w2 in global format
C--   coef(3:6) (in) input pdf is multiplied by coef(nf)
C--   iadd      (in) -1, 0, +1  =  substract, store , add

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension w1(*), w2(*), coef(3:6)
      dimension imin(6),imax(6)

C--   Find nf range
      nfmin = itfiz5(-1   )
      nfmax = itfiz5(-nzz5)
C--   Loop over nf
      do nf = nfmin,nfmax
        iz1 = izminf6(nf)
        iz2 = izmanf6(nf)
        cc  = coef(nf)
C--     Loop over z-grid
        do iz = iz1,iz2
          ia1 = iqcG5ijk(w1,1,iz,ig1)-1
          ia2 = iqcG5ijk(w2,1,iz,ig2)-1
C--       Inner loop over y-grid
          if(iadd.eq.-1)    then
            do iy = 1,nyy2(0)
              w2(ia2+iy) = w2(ia2+iy) - cc*w1(ia1+iy)
            enddo
          elseif(iadd.eq.0) then
            do iy = 1,nyy2(0)
              w2(ia2+iy) = cc*w1(ia1+iy)
            enddo
          elseif(iadd.eq.1) then
            do iy = 1,nyy2(0)
              w2(ia2+iy) = w2(ia2+iy) + cc*w1(ia1+iy)
            enddo
          else
            stop 'sqcPdfCopy: invalid iadd'
          endif
        enddo
      enddo

C--   Copy entries in satellite table
      call sqcGetLimits(w1,ig1,imin,imax,jmax)
      ia1 = iqcGSij(w1,1,ig1)-1
      ia2 = iqcGSij(w2,1,ig2)-1
C--   Copy jmax words
      do i = 1,jmax
        ia1     = ia1+1
        ia2     = ia2+1
        w2(ia2) = w1(ia1)
      enddo

*C--   Validate ig2
*      call sqcValidate(w2,ig2)

      return
      end

C=======================================================================
C==   Transformations between interpolation and evolution grids ========
C=======================================================================

C     =========================================
      subroutine sqcG0toGi(igg0,iggi,ig,nyg,iz)
C     =========================================

C--   Copy t-bin iz from interpolation grid G0 to subgrid Gi.
C--
C--   (in)  igg0  identifier grid G0 in global format
C--   (in)  iggi  identifier evolution subgrid inglobal format
C--   (in)  ig    index (i) = 1,...,nyg2  of the subgrid
C--   (in)  nyg   Upper yloop index in subgrid
C--   (in)  iz    z-bin in G0   to be copied 
C--
C--   For the y-grid index parameters see sqcGyMake in qcdgrd.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

C--   Base addresses
      iai = iqcG5ijk(stor7,1,iz,iggi)-1
      ia0 = iqcG5ijk(stor7,1,iz,igg0)-1
C--   Loop over ybins in Grid i
      do iyi = 1,nyg
        iy0            = iy0fiyg2(iyi,ig)
        stor7(iai+iyi) = stor7(ia0+iy0)
      enddo
     
      return
      end

C     ==================================
      subroutine sqcGitoG0(iggi,ig,igg0)
C     ==================================

C--   Copy F coefficients from Gi to interpolation grid G0
C--   The bin it=0 is not copied
C--
C--   (in)  iggi  identifier of evolution subgrid in global format
C--   (in)  ig    index (i) = 1,...,nyg2  of the subgrid
C--   (in)  igg0  identifier of grid G0 in global format
C--
C--   For the y-grid index parameters see sqcGyMake in qcdgrd.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

C--   Loop over z-grid
      do iz = 1,nzz2
C--     Base address in grid G0 and Gi
        ia0 = iqcG5ijk(stor7,1,iz,igg0)-1
        iai = iqcG5ijk(stor7,1,iz,iggi)-1
C--     Loop over subgrid points to copy
        do iyi = jymi2(ig),nyy2(ig)
          iy0            = iy0fiyg2(iyi,ig)
          stor7(ia0+iy0) = stor7(iai+iyi)
        enddo                          
      enddo

      return
      end

C     =========================================
      subroutine sqcGiFtoA(ig1,ig2,nyg,iz1,iz2)
C     =========================================

C--   Convert t-bins iz1 to iz2 from F to A and store result in id2
C--
C--   ig1   Input  subgrid table id in global format
C--   ig2   Output subgrid table id in global format
C--   nyg   Upper limit yloop
C--   iz1   First t-bin to convert (iz1,2 in ascending order)
C--   iz2   Last  t-bin to convert
C--   
C--   It is allowed to have id1 = id2 (in-place conversion)
C--   To convert all t-bins, set iz1 (iz2) very large negative (positive)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      if(iz2.lt.iz1) stop 'sqcGiFtoA: iz2 .lt. iz1'

C--   Index ranges
      call sqcPdfLims(ig1,iya,iyb,ita,itb,jmax)
      iy1 = 1
      iy2 = nyg
      it1 = max(ita,iz1)
      it2 = min(itb,iz2)
C--   Calculate base addresses
      inc = iqcG5ijk(stor7,iy1,it1+1,ig1)-
     +      iqcG5ijk(stor7,iy1,it1  ,ig1)
      ia1 = iqcG5ijk(stor7,iy1,it1,ig1)-inc
      ia2 = iqcG5ijk(stor7,iy1,it1,ig2)-inc
C--   Convert columns f --> a by solving f = Sa
      do j = it1,it2
        ia1 = ia1+inc
        ia2 = ia2+inc
        call sqcNSeqs(smaty2(1,ioy2),nmaty2(ioy2),stor7(ia2),
     +                stor7(ia1),iy2)
      enddo

      return
      end

C     =========================================
      subroutine sqcGiAtoF(ig1,ig2,nyg,iz1,iz2)
C     =========================================

C--   Convert t-bins iz1 to iz2 from A to F and store result in id2
C--
C--   ig1   Input  subgrid table id in global formet
C--   ig2   Output subgrid table id in global format
C--   nyg   Upper limit yloop
C--   iz1   First t-bin to convert (iz1,2 in ascending order)
C--   iz2   Last  t-bin to convert
C--   
C--   It is allowed to have id1 = id2 (in-place conversion via buffer)
C--   To convert all t-bins, set iz1 (iz2) very large negative (positive)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      if(iz2.lt.iz1) stop 'sqcGiAtoF: iz2 .lt. iz1'

C--   Index ranges
      call sqcPdfLims(ig1,iya,iyb,ita,itb,jmax)
      iy1 = 1
      iy2 = nyg
      it1 = max(ita,iz1)
      it2 = min(itb,iz2)
C--   Calculate base addresses
      inc = iqcG5ijk(stor7,iy1,it1+1,ig1)-
     +      iqcG5ijk(stor7,iy1,it1  ,ig1)
      ia1 = iqcG5ijk(stor7,iy1,it1,ig1)-inc
      ia2 = iqcG5ijk(stor7,iy1,it1,ig2)-inc
C--   Convert columns a --> f by multiplying f = Sa
      do j = it1,it2
        ia1 = ia1+inc
        ia2 = ia2+inc
        call sqcNSmult(smaty2(1,ioy2),nmaty2(ioy2),stor7(ia1),bufy7,iy2)
        do i = 1,iy2
          stor7(ia2-1+i) = bufy7(i)
        enddo
      enddo

      return
      end

C     ==============================================
      subroutine sqcAitoF0(iggi,ig,nyg,iz1,iz2,igg0)
C     ==============================================

C--   Transform A coefficients in Gi to F values in G0
C--   The bin it=0 is not copied
C--
C--   (in)  iggi  id of the equidistant subgrid in global format
C--   (in)  ig    index of the subgrid
C--   (in)  nyg   upper loop index in subgrid
C--   (in)  iz1   lower z index (iz1,2 in ascending order)
C--   (in)  iz2   upper z index
C--   (in)  igg0  identifier of the grid G0 in global format
C--
C--   For the y-grid index parameters see sqcGyMake in qcdgrd.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      if(iz2.lt.iz1) stop 'sqcAitoF0: iz2 .lt. iz1'

C--   Loop over z-grid
      do iz = iz1,iz2
C--     Base address in grid G0 and Gi
        ia0 = iqcG5ijk(stor7,1,iz,igg0)-1
        iai = iqcG5ijk(stor7,1,iz,iggi)
        call sqcNSmult(smaty2(1,ioy2),nmaty2(ioy2),stor7(iai),bufy7,nyg)
        do iyi = jymi2(ig),nyg
          iy0            = iy0fiyg2(iyi,ig)
          stor7(ia0+iy0) = bufy7(iyi)
        enddo  
      enddo      
      
      return
      end

C     ===========================================
      subroutine sqcGiLtoQ(idg1,idg2,nyg,iz1,iz2)
C     ===========================================

C--   Switch to quadratic interpolation and convert A coefficients
C--   of bins iz1 to iz2 from lin to quad 
C--
C--   idg1  Input  subgrid table id in global format (A *must* be lin)
C--   idg2  Output subgrid table id in global format (A will be quad)
C--   nyg   Upper limit yloop
C--   iz1   First t-bin to convert (iz1,2 in ascending order)
C--   iz2   Last  t-bin to convert
C--   
C--   It is allowed to have id1 = id2 (in-place conversion)
C--   To convert all t-bins, set iz1 (iz2) very large negative (positive)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      if(iz2.lt.iz1) stop 'sqcGiLtoQ: iz2 .lt. iz1'

C--   Convert lin A to F
      ioy2 = 2  !Tell AtoF that A is lin
      call sqcGiAtoF(idg1,idg2,nyg,iz1,iz2)
C--   Convert F to quad A
      ioy2 = 3
      call sqcGiFtoA(idg2,idg2,nyg,iz1,iz2)
      
      return
      end

C     ===========================================
      subroutine sqcGiQtoL(idg1,idg2,nyg,iz1,iz2)
C     ===========================================

C--   Switch to linear interpolation and convert A coefficients
C--   of bins iz1 to iz2 from quad to lin 
C--
C--   idg1  Input  subgrid table id in global format (A *must* be quad)
C--   idg2  Output subgrid table id in global format (A will be lin)
C--   nyg   Upper limit yloop
C--   iz1   First t-bin to convert (iz1,2 in ascending order)
C--   iz2   Last  t-bin to convert
C--   
C--   It is allowed to have id1 = id2 (in-place conversion)
C--   To convert all t-bins, set iz1 (iz2) very large negative (positive)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      if(iz2.lt.iz1) stop 'sqcGiQtoL: iz2 .lt. iz1'

C--   Convert quad A to F
      ioy2 = 3   !Tell AtoF that A is quad
      call sqcGiAtoF(idg1,idg2,nyg,iz1,iz2)
C--   Convert F to lin A
      ioy2 = 2
      call sqcGiFtoA(idg1,idg2,nyg,iz1,iz2)
      
      return
      end

C=======================================================================
C==   Spline coefficient transformations for convolution routines ======
C=======================================================================

C     ==============================================
      subroutine sqcGetSplA(ww,id,iy,iz,ig,iyg,aout)
C     ==============================================

C--   ww     (in)   store with pdf tables
C--   id     (in)   pdf table index (type 5) in global format
C--   iy     (in)   y-grid index in main grid G0
C--   iz     (in)   z-grid index in main grid G0
C--   ig     (out)  subgrid index
C--   iyg    (out)  y-grid index in subgrid Gi
C--   aout   (out)  array aout(1,...,iyg) of spline coefficients
C--                 defined on the subgrid Gi 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      dimension ww(*)
      dimension aout(*),buf(mxx0)

C--   Find subgrid
      ig = 1
      do i = 2,nyg2
        if(iy.gt.iyma2(i-1)) ig = i
      enddo
C--   Find y-grid index in subgrid (basically iy = y/del)
      iyg = iqcIyfrmY(ygrid2(iy),dely2(ig),nyy2(ig))
C--   Cross-check
      if(iy.ne.iy0fiyg2(iyg,ig)) then
       stop 'sqcGetSplA: problem y index in subgrid'
      endif
C--   Base address
      ia0 = iqcG5ijk(ww,1,iz,id)-1
C--   Copy pdf values to buffer
      do jy = 1,iyg
        iy0     = iy0fiyg2(jy,ig)
        buf(jy) = ww(ia0+iy0)
      enddo
C--   Now convert pdf values to spline coefficients
      call sqcNSeqs(smaty2(1,ioy2),nmaty2(ioy2),aout,buf,iyg)
      
      return
      end


