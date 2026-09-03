
C--   file srcFastEngine.f containing fast convolution engine routines

C--   integer function iqcTbufIjk(iy,it,ibufg)
C--   integer function iqcSbufIj(i,ibufg)
C--   subroutine sqcValidateBuf(ibufg)

C--   subroutine sqcFastBook(nwlast,ierr)

C--   subroutine sqcSetMark(xlist,qlist,n,margin,ierr)
C--   subroutine sqcFastInp(w,idf,ibg,iadd,factor,idense)
C--   subroutine sqcFastPdf(idg,coef,ibg,idense)
C--   subroutine sqcFastAdd(idg1,wt,n,idg2,
C--                         nzlist,izlist,nylist,iylist)
C--   subroutine sqcFastFxK(w,idwt,ibgi,ibgo,idense,ierr)
C--   subroutine sqcFccAtIt(ww,idwt,wi,ibgi,wo,ibgo,list,nl,iz)
C--   subroutine sqcFastWgt(w,idwt,iz,nf,ig,wmat)
C--   subroutine sqcFastFxF(w,idx,ibga,ibgb,ibgo,idense)
C--   subroutine sqcFcFAtIt(w,idx,ida,idb,ido,list,nl,iz)
C--   subroutine sqcFastKin(ibg,fun,idense)
C--   subroutine sqcFastCpy(ibg1,ibg2,iadd,idense)
C--   subroutine sqcFastFxq(w,idg,stf,n)

C--   subroutine sqcFastInt(w,idlst,coefs,m,ibg,x,q,f,n,ierr)
C--   subroutine sqcFastIntMpt(w,idlst,coefs,m,ibg,x,q,f,n,ierr)

C     ==================================================================
C     Fast convolution  ================================================
C     ==================================================================


C     ========================================
      integer function iqcTbufIjk(iy,it,ibufg)
C     ========================================

C--   Return address in stor7 of fast buffer table element
C--
C--   iy    (in) : y index [1-nyy2(0)]
C--   it    (in) : t index [1-ntt2]
C--   ibufg (in) : buffer id in global format

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'

      iqcTbufIjk = iqcG5ijk(stor7,iy,it,ibufg)

      return
      end

C     ===================================
      integer function iqcSbufIj(i,ibufg)
C     ===================================

C--   Return address in stor7 of fast buffer satellite element
C--
C--   i     (in) : attribute index 1=validation, 2=ymax, 3=tmin, 4=tmax
C--   ibufg (in) : buffer id in global format

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'

      iqcSbufIj = iqcGSij(stor7,i,ibufg)

      return
      end


C     ================================
      subroutine sqcValidateBuf(ibufg)
C     ================================

C--   Validate scratch buffer ibuf
C--
C--   ibufg (in) : buffer id in global format

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'

      call sqcValidate(stor7,ibufg)

      return
      end

C     ===================================
      subroutine sqcFastBook(nwlast,ierr)
C     ===================================

C--   Book set of fast engine scratch tables if not already booked
C--
C--   nwlast (out) : last word used in the store (<0 no space)
C--   ierr   (out) : 0 = OK
C--                 -1 = empty set of tables (never occurs)
C--                 -2 = not enough space
C--                 -3 = iset count MST0 exceeded
C--                 -4 = iset exist but with smaller n or different idfst
C--                 -5 = iset exist but with insufficient pointer tables

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

C--   Check tables already booked
      if(isetf7(-1).ne.0) then
        nwlast = iqcGetNumberOfWords(stor7)
        ierr   = 0
        return
      endif

C--   Now book pdf tables (w/o alfas tables); do nothing if tables exist
      iset  = -1
      ntab  =  mbf0
      ifrst =  1
      noalf =  1
      call sqcPdfBook(iset,ntab,ifrst,noalf,nwlast,ierr)
C--   Lfill7 has no meaning for this set but better set it to true
      Lfill7(iset) = .true.

      return
      end

C     ================================================
      subroutine sqcSetMark(xlist,qlist,n,margin,ierr)
C     ================================================

C--   Mark grid points
C--
C--   Input:  xlist      =  list of x   points
C--           qlist      =  list of mu2 points
C--           n          =  number of items in xlist and qlist
C--           margin     =  points to stay away from threshold [0,1]
C--
C--   Output: ierr       =  0    all OK
C--                         1    at least one x,q outside cuts
C--   Output in qfast9.inc:
C--
C--       mark9(0:mxx0,0:mqq0+7) = logical array with marked grid points
C--       xlst9,qlst9,nxq9       = copy of xlist, qlist and n
C--       ylst9,tlst9,nyt9       = contains only points inside grid
C--       ixqfyt9                = index from ylst,tlst --> xlst,qlst
C--       iy19,iy29,iz19,iz29    = interpolation mesh for each entry
C--                                in ylst,tlst
C--       nyy9, nzz9             = mesh size = interpolation order for
C--                                each entry in ylst,tlst
C--       wy9,wz9                = interpolation weights for each entry
C--                                in ylst,tlst

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qfast9.inc'

      dimension xlist(*), qlist(*)
      logical   lqcInside

C--   Check
      if(n.gt.mpt0) stop 'sqcSetMark: too many points n ---> STOP'

C--   Copy to common block
      nxq9 = n
      nyt9 = 0
      ierr = 0
      do i = 1,n
        xlst9(i) = xlist(i)
        qlst9(i) = qlist(i)
C--     Process only points inside grid        
        if(lqcInside(xlist(i),qlist(i))) then
          nyt9         =  nyt9+1   
          ylst9(nyt9)  = -log(xlist(i))
          tlst9(nyt9)  =  log(qlist(i))
C--       Remember position in original list
          ixqfyt9(nyt9) = i
        else
          ierr = 1            
        endif
      enddo

C--   Clear markers
      do j = 0,mqq0+7
        do i = 0,mxx0
          mark9(i,j) = .false.
        enddo
      enddo

C--   Put markers in the grid
      call sqcMarkyt(
     +     mark9,ylst9,tlst9,margin,iy19,iy29,iz19,iz29,it19,nyt9)

C--   Calculate interpolation weights
      do i = 1,nyt9
        nyy9(i) = iy29(i)-iy19(i)+1
        nzz9(i) = iz29(i)-iz19(i)+1
        call sqcIntWgt(
     +       iy19(i),nyy9(i),it19(i),nzz9(i),ylst9(i),tlst9(i),
     +       wy9(1,i),wz9(1,i))
      enddo
 
C--   Set up sparse and dense lists
C--   Loop over z points
      nz    = 0
      iymax = 0
      do iz = 1,nzz5
        ny = 0
        do iy = 1,nyy2(0)
          if(mark9(iy,iz)) then
            ny    = ny+1
            iymax = iy
          endif  
        enddo
        if(ny.ne.0) then
          nz           = nz+1
          izlist9(nz)  = iz
          nyslist9(nz) = ny
          nydlist9(nz) = iymax
        endif
        ny = 0
        do iy = 1,iymax
          iydlist9(iy,nz) = iy
          if(mark9(iy,iz)) then
            ny              = ny+1
            iyslist9(ny,nz) = iy
          endif                    
        enddo
      enddo
      nzlist9 = nz

C--   Number of markers set
*      nm = 0
*      do iz = 1,nzz5
*        do iy = 1,nyy2(0)
*          if(mark9(iy,iz)) nm = nm+1
*        enddo
*      enddo
*      write(6,*) 'sqcSETMARK nmark = ',nm
        
      return
      end

C     ===================================================
      subroutine sqcFastInp(w,idf,ibg,iadd,factor,idense)
C     ===================================================

C--   Multiply a pdf in w by factor and copy or add to ibuf in stor7
C--
C--   w                (in) : workspace
C--   idf              (in) : pdf identifier in w (global format)
C--   ibg              (in) : buffer id in stor7 (global format)
C--   iadd             (in) : store (0), add (1) or substract (-1)
C--   factor(3:6)      (in) : nf dependent weight factor
C--   idense           (in) : 0/1 no/yes dense table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      dimension w(*), factor(3:6)
      dimension sgn(-1:1)
      data sgn/-1.D0, 1.D0, 1.D0/

C--   Initialize
      if(iadd.eq.0) call sqcPreset(ibg,0.D0)

C--   Dense
      if(idense.eq.1) then
        do j = 1,nzlist9
          iz    = izlist9(j)
          nf    = itfiz5(-iz)
          fac   = sgn(iadd)*factor(nf)
          ia1   = iqcG5ijk(w,1,iz,idf)-1
          ia2   = iqcG5ijk(stor7,1,iz,ibg)-1
          if(iadd.eq.0) then
            do i = 1,nydlist9(j)
              iy = iydlist9(i,j)
              stor7(ia2+iy) = fac*w(ia1+iy)
            enddo
          else
            do i = 1,nydlist9(j)
              iy = iydlist9(i,j)
              stor7(ia2+iy) = stor7(ia2+iy) + fac*w(ia1+iy)
            enddo
           endif
        enddo
C--   Sparse
      else
        do j = 1,nzlist9
          iz    = izlist9(j)
          nf    = itfiz5(-iz)
          fac   = sgn(iadd)*factor(nf)
          ia1   = iqcG5ijk(w,1,iz,idf)-1
          ia2   = iqcG5ijk(stor7,1,iz,ibg)-1
          if(iadd.eq.0) then
            do i = 1,nyslist9(j)
              iy = iyslist9(i,j)
              stor7(ia2+iy) = fac*w(ia1+iy)
            enddo
          else
            do i = 1,nyslist9(j)
              iy = iyslist9(i,j)
              stor7(ia2+iy) = stor7(ia2+iy) + fac*w(ia1+iy)
            enddo
          endif
        enddo
      endif

      return
      end

C     ==========================================
      subroutine sqcFastPdf(idg,coef,ibg,idense)
C     ==========================================

C--   Copy linear combination of pdfs to ibg
C--
C--   idg              (in) : global id of gluon in stor7
C--   coef(3:6,0:12)   (in) : coefficients of the linear combination
C--   ibg              (in) : output buffer id in stor7 (global format)
C--   idense           (in) : 0/1 no/yes dense table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical lmb_ne
      dimension coef(3:6,0:12),cvec(3:6,12),idvec(12)

C--   Initialize
      call sqcPreset(ibg,0.D0)
C--   Weedout zero coefficients
      nvec = 0
      do i = 0,12
        istore = 0
        do  nf = 3,6
          if(lmb_ne(coef(nf,i),0.D0,aepsi6)) istore = 1
        enddo
        if(istore.eq.1) then
          nvec        = nvec+1
          if(nvec.gt.12) stop 'sqcFastPdf: nvec larger than 12'
          idvec(nvec) = idg+i
          do nf = 3,6
            cvec(nf,nvec) = coef(nf,i)
          enddo
        endif
      enddo 
C--   Summitup
      if(nvec.ne.0) then
        if(idense.eq.1) then
          call sqcFastAdd(idvec,cvec,nvec,ibg,
     +                    nzlist9,izlist9,nydlist9,iydlist9)
        else
          call sqcFastAdd(idvec,cvec,nvec,ibg,
     +                    nzlist9,izlist9,nyslist9,iyslist9)     
        endif
      endif

      return
      end

C     ==================================================
      subroutine sqcFastAdd(idg1,wt,n,idg2,
     +                      nzlist,izlist,nylist,iylist)
C     ==================================================

C--   Weighted sum of stor7 pdfs with weights dependent on nf.

C--   idg1(n)       list of input pdf ids in stor7 (global format)
C--   wt(3:6,n)     list of weights wt(nf,i) of pdf i
C--   n             number of pdfs to add
C--   idg2          output pdf id in stor7 (global format)
C--   nzlist        number of z grid points to loop over
C--   izlist(j)     grid point jz
C--   nylist(j)     number of y grid points to loop over at jz
C--   iylist((i,j)  grid point iy at jz

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension idg1(*),wt(3:6,*)
      dimension izlist(*),nylist(*),iylist(mxx0,*)

C--   Protect
      do i = 1,n
        if(idg1(i).eq.idg2) stop
     +      'sqcFastAdd: attempt to overwrite input id ---> STOP'
      enddo
C--   Initialize output table
      call sqcPreset(idg2,0.D0)
C--   Loop over table id's and fill target id
      do k = 1,n
        do j = 1,nzlist
          iz    = izlist(j)
          idgk  = idg1(k)
          ia1   = iqcG5ijk(stor7,1,iz,idgk)-1
          ia2   = iqcG5ijk(stor7,1,iz,idg2)-1
          nf    = itfiz5(-iz)
          wgt   = wt(nf,k)
          do i = 1,nylist(j)
            iy = iylist(i,j)
            stor7(ia2+iy) = stor7(ia2+iy) + wgt*stor7(ia1+iy)
          enddo
        enddo
      enddo
      
      return
      end

C     ===================================================
      subroutine sqcFastFxK(w,idwt,ibgi,ibgo,idense,ierr)
C     ===================================================

C--   Calculate FcrossC for a list of y,t values
C--
C--   w        (in)  : store with weight tables
C--   idwt(5)  (in)  : id(1)-(3) table ids LO, NLO, NNLO (0=no table)
C--                    id(4)     leading power of alfas  (0 or 1)
C--                    id(5)     key in pars8 where to get as table from
C--   ibgi     (in)  : buffer id with input pdf (global format)
C--   ibgo     (in)  : buffer id with output convolution (global format)
C--   idense   (in)  : 0/1 = no/yes dense output
C--   ierr     (out) : 0 = OK
C--                    1 = at least one grid point below alphaslim
C--
C--   idwt must be in global format

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      include 'qfast9.inc'

      dimension w(*),idwt(*)

C--   Loop over iz-bins
      ierr = 0
      nfxk = 0
      do j = 1,nzlist9
        iz = izlist9(j)
        it = itfiz5(iz)
C--     No alphas at such low t-bin
        itmin = int(dparGetPar(pars8,idwt(5),iditmin8))
        if(it.lt.itmin) ierr = 1
C--     Calculate stf for marked points and put result in ido
        if(idense.eq.0) then
          call sqcFccAtIt(w,idwt,
     +              stor7,ibgi,stor7,ibgo,iyslist9(1,j),nyslist9(j),iz)
          nfxk = nfxk + nyslist9(j)
        else
          call sqcFccAtIt(w,idwt,
     +              stor7,ibgi,stor7,ibgo,iydlist9(1,j),nydlist9(j),iz)
          nfxk = nfxk + nydlist9(j)
        endif
      enddo

C--   Number of convolutions
*      write(6,*) 'sqcFASTFXK nconvol = ',nfxk

      return
      end

C     =======================================================
      subroutine sqcFccAtIt(ww,idwt,wi,idi,wo,ido,list,nl,iz)
C     =======================================================

C--   Calculate FcrossC for a list of grid points iy at fixed iz

C--   ww      (in) : workspace containing weight tables
C--   idwt(5) (in) : weight table ids in global format
C--   wi      (in) : workspace with input pdf tables
C--   idi     (in) : id of input pdf table in global format
C--   wo      (in) : workspace with output pdf tables
C--   ido     (in) : identifier of output fxk table in global format
C--   list    (in) : list of grid points in y
C--   nl      (in) : number of grid points in the list
C--   iz      (in) : t-grid index

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension ww(*),wi(*),wo(*),coef(mxx0)
      dimension list(mxx0)
      dimension idwt(*)
      dimension wmat(mxx0)

C--   Value of t
      it = itfiz5( iz)
      nf = itfiz5(-iz)
      tt = tgrid2(it)

C--   Base address
      iao = iqcG5ijk(wo,1,iz,ido)-1

C--   Check list in ascending order
      if(list(nl).lt.list(1)) stop 'sqcFccAtIt: wrong y-loop'
C--   Loop over all requested grid points
      iglast = 0
C--   Must loop from high iy to low iy!      
      do ii = nl,1,-1
        jy  = list(ii)
        if(jy.eq.0) then
C--       Lower edge of y-grid: set convolution to zero
          fyj = 0.D0
        else
C--       Gridpoint above lower edge of y-grid
          yj  = ygrid2(jy)
          ig  = iqcFindIg(yj)
C--       New subgrid encountered
          if(ig.ne.iglast) then
C--         Fill weight table
            call sqcFastWgt(ww,idwt,iz,nf,ig,wmat)
C--         Convert F values to A values
            call sqcGetSplA(wi,idi,jy,iz,ig,iyg,coef)
            iglast = ig
          endif
C--       Done with setting-up subgrid; get y-index in subgrid 
          ky  = iqcIyfrmY(yj,dely2(ig),nyy2(ig))
C--       Convolution loop in subgrid
          fyj = 0.D0
          do i = 1,ky
            fyj = fyj + wmat(ky+1-i)*coef(i)
          enddo
        endif
C--     Store FxK in output table
        if(it.ge.itmin6) then
          wo(iao+jy) = fyj      !alphas available for this t-bin
        else
          wo(iao+jy) = qnull6   !no alphas for this t-bin
        endif
C--     End of loop over requested grid points
      enddo

      return
      end

C     ===========================================
      subroutine sqcFastWgt(w,idwt,iz,nf,ig,wmat)
C     ===========================================

C--   Return weight table or perturbative expansion of tables
C--
C--   w        (in) : store containing weight tables
C--   idwt(5)  (in) : id(1)-(3) table ids LO, NLO, NNLO (0=no table)
C--                   id(4)     leading power of alfas  (0 or 1)
C--                   id(5)     pars8 key to get the as table from
C--   iz       (in) : z-grid index
C--   nf       (in) : number of flavors
C--   ig       (in) : subgrid index
C--   wmat    (out) : weight table
C--
C--   NB: w must be in global format

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      include 'qpars6.inc'

      logical lqcIdExists

      dimension w(*),idwt(*),wmat(*),mi(6),ma(6)
      dimension ialf(-mord0:mord0)

      it = itfiz5(iz)

C--   Alfas table base address
      key  = idwt(5)
      do i = -mord0,mord0
        jd      = 1000*key + 601 + i + mord0
        ialf(i) = iqcG6ij(pars8,1,jd)-1
      enddo
      
      do iy = 1,nyy2(ig)
        wmat(iy) = 0.D0
      enddo
      
      do io = 1,iord6
        id = idwt(io)
        if(id.ne.0) then
          if(.not.lqcIdExists(w,id))
     +         stop 'sqcFastWgt: no weight table of this type'
          call sqcGetLimits(w,id,mi,ma,jmax)
          jt = max(it,mi(2))
          jt = min(jt,ma(2))
          kf = max(nf,mi(3))
          kf = min(kf,ma(3))
          ia = iqcGaddr(w,1,jt,kf,ig,id)-1
C--       Fill weight table
          if(idwt(4).eq.0) then
C--         Multiply LO,NLO,NNLO by 1,as,as2
            if(io.eq.1) then
              do iy = 1,nyy2(ig)
                ia = ia + 1
                wmat(iy) = wmat(iy) + w(ia)
              enddo  
            else
              do iy = 1,nyy2(ig)
                ia = ia + 1
                wmat(iy) = wmat(iy) + pars8(ialf(-(io-1))+iz)*w(ia)
              enddo
            endif
          else
C--         Multiply LO,NLO,NNLO by as,as2,as3
            do iy = 1,nyy2(ig)
              ia = ia + 1
              wmat(iy) = wmat(iy) + pars8(ialf(io)+iz)*w(ia)
            enddo                 
          endif
        endif
      enddo

      return
      end
      
C     ==================================================
      subroutine sqcFastFxF(w,idx,ibga,ibgb,ibgo,idense)
C     ==================================================

C--   Calculate FcrossF for a list of y,t values
C--
C--   w       (in)  : store with weight tables
C--   idx     (in)  : id (global format) of MakeWtX weight table
C--   ibga,b  (in)  : id of buffer with input pdf (global format)
C--   ibgo    (in)  : id of buffer with output convolution (global format)
C--   idense  (in)  : 0/1 = no/yes dense output
C--   ierr    (out) : 0 = OK
C--                    1 = at least one grid point below alphaslim

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qfast9.inc'
      include 'qstor7.inc'
      
      dimension w(*)

      if(idense.lt.0 .or. idense.gt.1) stop 'sqcFastFxF wrong idense'

C--   Loop over iz-bins
      do j = 1,nzlist9
        iz = izlist9(j)
        it = itfiz5(iz)
C--     Calculate stf for marked points and put result in idout
        if(idense.eq.0) then
          call sqcFcFAtIt(w,idx,stor7,ibga,stor7,ibgb,stor7,ibgo,
     +                    iyslist9(1,j),nyslist9(j),iz)
        else
          call sqcFcFAtIt(w,idx,stor7,ibga,stor7,ibgb,stor7,ibgo,
     +                    iydlist9(1,j),nydlist9(j),iz) 
        endif
      enddo
      
      return
      end

C     =============================================================
      subroutine sqcFcFAtIt(wx,idx,wa,ida,wb,idb,wo,ido,list,nl,iz)
C     =============================================================

C--   Calculate FcrossF for a list of grid points iy at fixed iz

C--   Input: wx       store containing weight tables
C--          idx      weight table id in global format
C--          wa,ida   pdf store and pdf identifier in global format
C--          wb,idb   pdf store and pdf identifier in global format
C--          wo,ido   pdf store and identifier of output fxf table
C--          list     list of grid points in y
C--          nl       number of grid points in the list
C--          iz       t-grid index

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension wx(*),wa(*),wb(*),wo(*),coefa(mxx0),coefb(mxx0)
      dimension list(mxx0)

C--   Value of t
      it = itfiz5( iz)
      nf = itfiz5(-iz)
      tt = tgrid2( it)

C--   Base address
      iao = iqcG5ijk(wo,1,iz,ido)-1
      
C--   Check list in ascending order
      if(list(nl).lt.list(1)) stop 'sqcFcFAtIt: wrong y-loop'
C--   Loop over all requested grid points
      iglast = 0
C--   Must loop from high iy to low iy!      
      do ii = nl,1,-1
        jy  = list(ii)
        if(jy.eq.0) then
C--       Lower edge of y-grid: set convolution to zero
          fxf = 0.D0
        else
C--       Gridpoint above lower edge of y-grid
          yj  = ygrid2(jy)
          ig  = iqcFindIg(yj)
          iwx = 0                   !avoid compiler warning
C--       New subgrid encountered
          if(ig.ne.iglast) then
C--         Convert F values to A values
            call sqcGetSplA(wa,ida,jy,iz,ig,iyg,coefa)
            call sqcGetSplA(wb,idb,jy,iz,jg,jyg,coefb)
            iglast = ig
C--         Weight table base address
            iwx = iqcGaddr(wx,1,it,nf,ig,idx)-1
          endif
C--       Done with setting-up subgrid; get y-index in subgrid 
          ky  = iqcIyfrmY(yj,dely2(ig),nyy2(ig))
C--       Convolution loop in subgrid
          fxf = 0.D0
          do j = 1,ky
            Aj = coefa(j)
            do k = 1,ky-j+1
              Bk  = coefb(k) 
              fxf = fxf + Aj*Bk*wx(iwx+ky-j-k+2)
            enddo
          enddo
        endif
C--     Store FxF in output table
        wo(iao+jy) = fxf
C--     End of loop over requested grid points
      enddo

      return
      end
      
C     =====================================
      subroutine sqcFastKin(ibg,fun,idense)
C     =====================================

C--   Multiply stf with x-q dependent factor
C--
C--   Input:  ibg    =  scratch buffer id in global format
C--           fun    =  function of ix, iq, nf and ithres
C--           idense =  0/1 no/yes dense output

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      external fun

C--   Loop over grid points
      if(idense.eq.0) then                           !sparse output
        do j     = 1,nzlist9
          iz     = izlist9(j)
          it     = itfiz5( iz)
          nf     = itfiz5(-iz)
          ithres = 0
          if(iz.ne.nzz5) then
            nfp1 = itfiz5(-(iz+1))
            if(nfp1.eq.nf+1) ithres = -1
          endif
          if(iz.ne.   1) then
            nfm1 = itfiz5(-(iz-1))
            if(nfm1.eq.nf-1) ithres =  1
          endif
          ia = iqcG5ijk(stor7,1,iz,ibg)-1
          do i = 1,nyslist9(j)
            iy = iyslist9(i,j)
            ix = nyy2(0)-iy+1
            stor7(ia+iy) = stor7(ia+iy)*fun(ix,it,nf,ithres)
          enddo
        enddo
      else                                           !dense output
        do j = 1,nzlist9
          iz = izlist9(j)
          it = itfiz5( iz)
          nf = itfiz5(-iz)
          ithres = 0
          if(iz.ne.nzz5) then
            nfp1 = itfiz5(-(iz+1))
            if(nfp1.eq.nf+1) ithres = -1
          endif
          if(iz.ne.   1) then
            nfm1 = itfiz5(-(iz-1))
            if(nfm1.eq.nf-1) ithres =  1
          endif
          ia = iqcG5ijk(stor7,1,iz,ibg)-1
          do i = 1,nydlist9(j)
            iy = iydlist9(i,j)
            ix = nyy2(0)-iy+1
            stor7(ia+iy) = stor7(ia+iy)*fun(ix,it,nf,ithres)
          enddo
        enddo  
      endif  

      return
      end

C     ============================================
      subroutine sqcFastCpy(ibg1,ibg2,iadd,idense)
C     ============================================

C--   Add content of buffer ibg1 to that of ibg2 (ibg in global format)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'
      
      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      if(idense.lt.0 .or. idense.gt.1) stop 'sqcFastCpy wrong idense'

      if(idense.eq.0) then                        !sparse output
      
        if(iadd.eq.-1) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcG5ijk(stor7,1,iz,ibg1)-1
            ia2   = iqcG5ijk(stor7,1,iz,ibg2)-1
            do i = 1,nyslist9(j)
              iy = iyslist9(i,j) 
              stor7(ia2+iy) = stor7(ia2+iy)-stor7(ia1+iy)
            enddo
          enddo
        elseif(iadd.eq.0) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcG5ijk(stor7,1,iz,ibg1)-1
            ia2   = iqcG5ijk(stor7,1,iz,ibg2)-1
            do i = 1,nyslist9(j)
              iy = iyslist9(i,j)
              stor7(ia2+iy) = stor7(ia1+iy)
            enddo
          enddo  
        elseif(iadd.eq.1) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcG5ijk(stor7,1,iz,ibg1)-1
            ia2   = iqcG5ijk(stor7,1,iz,ibg2)-1
            do i = 1,nyslist9(j)
              iy = iyslist9(i,j)
              stor7(ia2+iy) = stor7(ia2+iy)+stor7(ia1+iy)
            enddo
          enddo  
        else
          stop 'sqcFastCpy: invalid iadd'   
        endif
      
      else                                        !dense output
      
        if(iadd.eq.-1) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcG5ijk(stor7,1,iz,ibg1)-1
            ia2   = iqcG5ijk(stor7,1,iz,ibg2)-1
            do i = 1,nydlist9(j)
              iy = iydlist9(i,j) 
              stor7(ia2+iy) = stor7(ia2+iy)-stor7(ia1+iy)
            enddo
          enddo
        elseif(iadd.eq.0) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcG5ijk(stor7,1,iz,ibg1)-1
            ia2   = iqcG5ijk(stor7,1,iz,ibg2)-1
            do i = 1,nydlist9(j)
              iy = iydlist9(i,j)
              stor7(ia2+iy) = stor7(ia1+iy)
            enddo
          enddo  
        elseif(iadd.eq.1) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcG5ijk(stor7,1,iz,ibg1)-1
            ia2   = iqcG5ijk(stor7,1,iz,ibg2)-1
            do i = 1,nydlist9(j)
              iy = iydlist9(i,j)
              stor7(ia2+iy) = stor7(ia2+iy)+stor7(ia1+iy)
            enddo
          enddo  
        else
          stop 'sqcFastCpy: invalid iadd'   
        endif              
        
      endif

      return
      end
      
C     ==================================
      subroutine sqcFastFxq(w,idg,stf,n)
C     ==================================

C--   Do n interpolations on a pdf table

C--   w    (in)  : workspace with pdf tables
C--   idg  (in)  : global identifier of pdf table in w
C--   stf  (out) : array with interpolated results
C--   n    (in)  : requested number of interpolations
C--
C--   The routine fills min(n,nyt9) entries in the array stf which
C--   must be dimensioned to at least stf(n) in the calling routine
C--
C--   The interpolation mesh should have been set-up by a call to
C--   sqcSetMark before the call to sqcFastFxq. The routine sqcSetMark
C--   passes the following information via the common block qfast9.inc:
C--
C--   ylst9(i)  (in) :  list of y points
C--   tlst9(i)  (in) :  list of t points
C--   iy19(i)   (in) :  list of lower mesh limits in y
C--   iy29(i)   (in) :  list of upper mesh limits in y
C--   iz19(i)   (in) :  list of lower mesh limits in z
C--   iz29(i)   (in) :  list of upper mesh limits in z
C--   nyy9(i)   (in) :  mesh size = interpolation order in y
C--   nzz9(i)   (in) :  mesh size = interpolation order in z
C--   wy9(j,i)  (in) :  interpolation weights in y
C--   wz9(j,i)  (in) :  interpolation weights in z
C--   ixqfyt(i) (in) :  pointer to index in stf(n)
C--   nyt9      (in) :  number of items in the lists above

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qfast9.inc'

      dimension w(*),stf(*)

      if(n.le.0) stop 'sqcFastFxq wrong n'

C--   Preset for points outside grid
      do i = 1,n
        stf(i) = qnull6
      enddo   
C--   Now interpolate for points inside grid
      do i = 1,min(n,nyt9)
        ia = iqcG5ijk(w,iy19(i),iz19(i),idg)
        ny = nyy9(i)
        nz = nzz9(i)
        stf(ixqfyt9(i)) = dqcPdfPol(w,ia,ny,nz,wy9(1,i),wz9(1,i))
      enddo
      
      return
      end

C=======================================================================
C==   Higher level routines ============================================
C=======================================================================

C     =======================================================
      subroutine sqcFastInt(w,idlst,coefs,m,ibg,x,q,f,n,ierr)
C     =======================================================

C--   Interpolate weighted sum of pdfs (wrapper routine)
C--
C--   w               (in) : store
C--   idlst(m)        (in) : list of global pdf identifiers
C--   coefs(3:6,m)    (in) : list of coefficients
C--   m               (in) : number of pdfs to add
C--   ibg             (in) : global id of buffer for temporary storage
C--   x,q             (in) : interpolation points in x and mu2
C--   f              (out) : list of interpolated results
C--   n               (in) : number of points in x,q,f (no mpt0 limit)
C--   ierr            (in) : 0=OK, 1=at least one x,q outside cuts
C--
C--   NB: is is assumed that fast buffers are already booked

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension w(*), idlst(*), coefs(3:6,*), x(*), q(*), f(*)

      if(isetf7(-1).eq.0) stop 'sqcFastInt: no fast buffers booked'

      ierr  = 0
      nlast = 0
      ntodo = min(n,mpt0)
      do while( ntodo.gt.0 )
        i1 = nlast+1
        call sqcFastIntMpt(
     +       w,idlst,coefs,m,ibg,x(i1),q(i1),f(i1),ntodo,jerr)
        nlast = nlast+ntodo
        ntodo = min(n-nlast,mpt0)
        ierr  = max(ierr,jerr)
      enddo

      return
      end

C     ==========================================================
      subroutine sqcFastIntMpt(w,idlst,coefs,m,ibg,x,q,f,n,ierr)
C     ==========================================================

C--   Interpolate weighted sum of pdfs (with mpt0 limit)
C--
C--   w               (in) : store
C--   idlst(m)        (in) : list of global pdf identifiers
C--   coefs(3:6,m)    (in) : list of coefficients
C--   m               (in) : number of pdfs to add
C--   ibg             (in) : global id of buffer for temporary storage
C--   x,q             (in) : interpolation points in x and mu2
C--   f              (out) : list of interpolated results
C--   n               (in) : number of points in x,q,f (<= mpt0)
C--   ierr            (in) : 0=OK, 1=at least one x,q outside cuts
C--
C--   NB: is is assumed that fast buffers are already booked

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension w(*), idlst(*), coefs(3:6,*), x(*), q(*), f(*)

      call sqcSetMark(x,q,n,0,ierr)

      iadd   =  0
      idense =  0
      call sqcFastInp(w,idlst(1),ibg,iadd,coefs(3,1),idense)
      iadd   =  1
      do i = 2,m
        call sqcFastInp(w,idlst(i),ibg,iadd,coefs(3,i),idense)
      enddo
      call sqcFastFxq(stor7,ibg,f,n)

      return
      end

