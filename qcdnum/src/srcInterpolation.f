
C--   This is the file srcInterpolation.f with interpolation routines 
C--
C--   subroutine sqcLstIni(yy,tt,npt,ww,ndim,nused,ierr)
C--   subroutine sqcFillBuffer(fun,wtb,par,n,ww,ierr,jerr)
C--   subroutine sqcFillBuffij(fun,ww,ierr)
C--   subroutine sqcLstFun(ww,ff,n,nout,ierr)
C--   subroutine sqcTabIni(yy,ny,tt,nt,ww,ndim,nused,ierr)
C--   subroutine sqcTabFun(ww,ff,ierr)
C--
C--   function   dqcPdfSum(w,par,np,ix,iq,nf,ithr,ia,ff,first,jerr)
C--
C--   subroutine sqcIntWgt(iy1,ny,it1,nt,y,t,wy,wz)
C--   function   dqcPdfPol(w,ia0,ny,nz,wy,wz)
C--
C--   subroutine sqcZmesh(y,t,margin,iymi,iyma,izmi,izma,itmi)
C--   subroutine sqcZmeshy(y,iymi,iyma)
C--   subroutine sqcZmesht(t,margin,izmi,izma,itmi)
C--   subroutine sqcMarkyt(mark,ylst,tlst,margin,iy1,iy2,iz1,iz2,it1,n)
C--   subroutine sqcMarkyy(marky,ylst,iy1,iy2,n)
C--   subroutine sqcMarktt(markz,tlst,margin,iz1,iz2,it1,n)
C--
C--   subroutine sqcPolint(xa,ya,n,x,y,dy)                      not used
C--   subroutine sqcPolin2(xa,nx,ya,ny,za,x,y,z)                not used
C--   subroutine sqcGetABC(u,v,w,delta,a,b,c)

C=======================================================================
C==== Fast interpolation routines ======================================
C=======================================================================

C     ==================================================
      subroutine sqcLstIni(yy,tt,npt,ww,ndim,nused,ierr)
C     ==================================================

C--   Initialise list-interpolation by setting up two arrays in ww:
C--
C--   ArrayL(17,n) list of meshes and weights for all points
C--   ArrayF( 3,m) list of gridpoints to be filled with data
C--
C--   yy,tt  (in) : interpolation points (all inside grid or cuts)
C--   npt    (in) : number of interpolation points
C--   ww    (out) : working array
C--   ndim   (in) : dimension of ww declared in the calling routine
C--   nused (out) : number of words used in ww
C--   ierr  (out) : 0 = OK
C--                 1 = Not enough space in ww (fatal error)
C--                 2 = No temporary buffer available
C--
C--   Layout of ww
C--   ------------
C--   ww(01)    =  Initialisation flag (123456)
C--   ww(02)    =  Evolution parameter key
C--   ww(03)    =  Scratch buffer id (0 = no buffer assigned)
C--   ww(04)    =  Number of points ArrayF
C--   ww(05)    =  K0 of ArrayF
C--   ww(06)    =  K1 of ArrayF
C--   ww(07)    =  K2 of ArrayF
C--   ww(08)    =  Number of interpolation points npt
C--   ww(09)    =  K0 of ArrayL
C--   ww(10)    =  K1 of ArrayL
C--   ww(11)    =  K2 of ArrayL
C--
C--   ww(12)    =  First word of ArrayL
C--   ..           ..
C--   ww(nused) =  Last  word of ArrayF
C--
C--   Layout of ArrayL(i,j)
C--   ---------------------
C--   ArrayL(01,j) = iy of mesh of interpolation point j
C--   ArrayL(02,j) = iz of mesh
C--   ArrayL(03,j) = ia relative address of (iy,iz) in pdf table
C--   ArrayL(04,j) = ny mesh width = interpolation order in y
C--   ArrayL(05,j) = nz mesh width = interpolation order in z
C--   ArrayL(06,j) = wy(1) first word of interpolation weights wy(6)
C--   ..
C--   ArrayL(12,j) = wz(1) first word of interpolation weights wz(6)
C--   ..
C--
C--   Layout of ArrayF(i,j)
C--   ---------------------
C--   ArrayF(01,j) = iy of data point j
C--   ArrayF(02,j) = iz of data point j
C--   ArrayF(03,j) = ia relative address in pdf table
C--
C--   The max size of ww is (with 9 mesh points / interpolation point):
C--
C--   nused = 11 + 17*n + 3*9*n = 11 + 44*n (or less)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension yy(*), tt(*), ww(*)
      dimension imin(2), imax(2), karrL(0:2), karrF(0:2)

C--   Address functions in ww
      iaL(i,j) = int(ww( 9)) + int(ww(10))*i + int(ww(11))*j
      iaF(i,j) = int(ww( 5)) + int(ww( 6))*i + int(ww( 7))*j

C--   Initialize
      ierr = 0

C--   Book ArrayL in ww
      imin(1) = 1
      imax(1) = 17
      imin(2) = 1
      imax(2) = npt
      istart  = 12
      call smb_bkmat(imin,imax,karrL,2,istart,ilast)

C--   Book arrayF in ww (assuming 9 mesh points / interpolation point)
      imin(1) = 1
      imax(1) = 3
      imin(2) = 1
      imax(2) = 9*npt
      call smb_bkmat(imin,imax,karrF,2,ilast+1,nused)

C--   Check size of ww
      if(nused.gt.ndim) then
        ierr = 1
        return
      endif

C--   Store partition parameters in ww
      ww( 1) = 0.D0              !initialisation flag
      ww( 9) = dble(karrL(0))
      ww(10) = dble(karrL(1))
      ww(11) = dble(karrL(2))
      ww( 5) = dble(karrF(0))
      ww( 6) = dble(karrF(1))
      ww( 7) = dble(karrF(2))

C--   Reserve a scratch buffer and calculate base address
      idbuf = iqcGimmeScratch()
      if(idbuf.eq.0) then
        ierr = 2
        return
      endif
      ibase = iqcG5ijk(stor7,1,1,idbuf)

C--   Loop over the interpolation points
      nff = 0
      do i = 1,npt
C--     Find interpolation mesh
        call sqcZmesh(yy(i),tt(i),0,iy1,iy2,iz1,iz2,it1)
        ny = iy2-iy1+1
        nz = iz2-iz1+1
C--     Address of first mesh point
        ia = iqcG5ijk(stor7,iy1,iz1,idbuf)
C--     Store the mesh of each interpolation point
        jj        = iaL(1,i)-1
        ww(jj+ 1) = dble(iy1)
        ww(jj+ 2) = dble(iz1)
        ww(jj+ 3) = dble(ia-ibase)
        ww(jj+ 4) = dble(ny)
        ww(jj+ 5) = dble(nz)
C--     Store interpolation weights
        call sqcIntWgt(iy1,ny,it1,nz,yy(i),tt(i),ww(jj+6),ww(jj+12))
C--     Update the list of all mesh points to be filled with data
        iaz = ia-inciz7
        do iz = iz1,iz2
          iaz = iaz+inciz7
          iay = iaz-1
          do iy = iy1,iy2
            iay = iay+1
            if(int(stor7(iay)).ne.1) then !point not yet stored
              stor7(iay) = dble(1)
              nff        = nff+1
              kk         = iaF(1,nff)-1
              ww(kk+1)   = dble(iy)
              ww(kk+2)   = dble(iz)
              ww(kk+3)   = dble(iay-ibase)
            endif
          enddo
        enddo
      enddo
C--   End of loop over the interpolation points

C--   Flag ww initialized and store number of points
      ww(1) = dble(123456)   !initflag
      ww(2) = dble(ipver6)   !parameter key
      ww(3) = dble(0)        !no buffer assigned
      ww(4) = dble(nff)
      ww(8) = dble(npt)

C--   Release scratch buffer
      call sqcReleaseScratch(idbuf)

      return
      end

C     ====================================================
      subroutine sqcFillBuffer(fun,wtb,par,n,ww,ierr,jerr)
C     ====================================================

C--   Fill scratch buffer with function values.
C--
C--   fun    (in) : external function
C--   wtb    (in) : toolbox workspace passed to fun
C--   par    (in) : parameter array passed to fun
C--   n      (in) : number of parameters passed to fun
C--   ww     (in) : workspace previously filled by sqcLstIni
C--   ierr  (out) : 0 = all OK
C--                 1 = ww not initialised
C--                 2 = ww not set-up with current parameters
C--                 3 = no scratch buffer available in memory
C--                 4 = encountered error from fun at loop entry
C--   jerr  (out) : error code from fun
C--
C--   This routine assigns a scratch buffer, if not assigned aready.
C--   The buffer will be released in the call to sqcLstFun.
C--
C--   Argument list of fun
C--
C--   val = fun(wtb,par,n,ix,iq,nf,ithr,ia,val,first,ierr)
C--
C--   wtb    (in) : toolbox workspace (dummy if workspace not used)
C--   par    (in) : parameter list (dummy if no pars)
C--   n      (in) : number of parameters (0 = no pars)
C--   ix,iq  (in) : grid point
C--   nf     (in) : number of flavours at (ix,iq)
C--   ithr   (in) : +-1 = threshold with upper (lower) nf, 0 otherwise
C--   ia     (in) : address of (ix,iq) in pdf table, relative to (1,1)
C--   val    (in) : current value stored at (ix,iq)
C--   first  (in) : true = first call in the loop over (ix,iq)
C--   ierr  (out) : error code set by fun (0 = OK)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      external fun

      dimension wtb(*),par(*),ww(*)

      logical first

C--   Inline address function in ww
      iaF(i,j) = int(ww(5)) + int(ww(6))*i + int(ww(7))*j

      ierr = 0
      jerr = 0

C--   Check ww initialised
      if(int(ww(1)).ne.123456 .and. int(ww(1)).ne.654321) then
        ierr = 1
        return
      endif
C--   Check for evolution parameter change
      if(int(ww(2)).ne.ipver6) then
        ierr = 2
        return
      endif
C--   Check if scratch buffer assigned
      idbuf = int(ww(3))
      if(idbuf.eq.0) then
C--     No buffer assigned, so get one
        idbuf = iqcGimmeScratch()
C--     No more buffers available in stor7
        if(idbuf.eq.0) then
          ierr = 3
          return
        endif
        ww(3) = dble(idbuf)
      endif

C--   Buffer base address
      ibase = iqcG5ijk(stor7,1,1,idbuf)

C--   Loop over mesh points
      nff = int(ww(4))

      do i = 1,nff
        if(i.eq.1) then
          first = .true.
        else
          first = .false.
        endif
        jj   = iaF(1,i)-1
        iy   = int(ww(jj+1))
        iz   = int(ww(jj+2))
        ia   = int(ww(jj+3))
        ix   = nyy2(0)-iy+1
        it   = itfiz5( iz)
        nf   = itfiz5(-iz)
        ithr = 0
        if(iz.ne.nzz5) then
          nfp1 = itfiz5(-(iz+1))
          if(nfp1.eq.nf+1) ithr = -1
        endif
        if(iz.ne.1   ) then
          nfm1 = itfiz5(-(iz-1))
          if(nfm1.eq.nf-1) ithr =  1
        endif
        ff              = stor7(ia+ibase)
        stor7(ia+ibase) = fun(wtb,par,n,ix,it,nf,ithr,ia,ff,first,jerr)
        if(first .and. jerr.ne.0) then
          ierr = 4
          return
        endif
      enddo

*mb
*mb      call sqcReleaseScratch(idbuf)
*mb

      return
      end

C     =====================================
      subroutine sqcFillBuffij(fun,ww,ierr)
C     =====================================

C--   Fill scratch buffer with function values
C--   As sqcFillBuffer but fun is a function of ix and iq only
C--
C--   fun    (in) : external function
C--   ww     (in) : workspace previously filled by sqcLstIni
C--   ierr  (out) : 0 = all OK
C--                 1 = ww not initialised
C--                 2 = ww not set-up with current parameters
C--                 3 = no scratch buffer available in memory
C--
C--   This routine assigns a scratch buffer, if not assigned aready.
C--   The buffer will be released in the call to sqcLstFun.
C--

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      external fun

      dimension ww(*)

C--   Inline address function in ww
      iaF(i,j) = int(ww(5)) + int(ww(6))*i + int(ww(7))*j

      ierr = 0

C--   Check ww initialised
      if(int(ww(1)).ne.123456 .and. int(ww(1)).ne.654321) then
        ierr = 1
        return
      endif
C--   Check for evolution parameter change
      if(int(ww(2)).ne.ipver6) then
        ierr = 2
        return
      endif
C--   Check if scratch buffer assigned
      idbuf = int(ww(3))
      if(idbuf.eq.0) then
C--     No buffer assigned, so get one
        idbuf = iqcGimmeScratch()
C--     No more buffers available in stor7
        if(idbuf.eq.0) then
          ierr = 3
          return
        endif
        ww(3) = dble(idbuf)
      endif

C--   Buffer base address
      ibase = iqcG5ijk(stor7,1,1,idbuf)

C--   Loop over mesh points
      nff = int(ww(4))

      do i = 1,nff
        jj   = iaF(1,i)-1
        iy   = int(ww(jj+1))
        iz   = int(ww(jj+2))
        ia   = int(ww(jj+3))
        ix   = nyy2(0)-iy+1
        it   = itfiz5( iz)
        nf   = itfiz5(-iz)
        ithr = 0
        if(iz.ne.nzz5) then
          nfp1 = itfiz5(-(iz+1))
          if(nfp1.eq.nf+1) ithr = -1
        endif
        if(iz.ne.1   ) then
          nfm1 = itfiz5(-(iz-1))
          if(nfm1.eq.nf-1) ithr =  1
        endif
        if(ithr.ne.-1) then
          stor7(ia+ibase) = fun(ix,it)
        else
          stor7(ia+ibase) = fun(ix,-it)
        endif
      enddo

*mb
*mb      call sqcReleaseScratch(idbuf)
*mb

      return
      end

C     =======================================
      subroutine sqcLstFun(ww,ff,n,nout,ierr)
C     =======================================

C--   Interpolate list of function values
C--
C--   ww     (in) : workspace previously filled by sqcLstIni
C--   ff    (out) : list of interpolated function values
C--   n      (in) : dimension of ff declared in the calling routine
C--   nout  (out) : number of values returned in ff
C--   ierr  (out) : 0 = all OK
C--                 1 = ww not initialized
C--                 2 = ww not set-up with current parameters
C--                 3 = no buffer with function values available
C--
C--   The function values reside in a buffer with id stored in ww(3)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension ww(*), ff(*)

C--   Address functions in ww
      iaL(i,j) = int(ww(9)) + int(ww(10))*i + int(ww(11))*j

      ierr = 0

C--   Check ww initialized
      if(int(ww(1)).ne.123456) then
        ierr = 1
        return
      endif
C--   Check for evolution parameter change
      if(int(ww(2)).ne.ipver6) then
        ierr = 2
        return
      endif
C--   Check if scratch buffer assigned
      idbuf = int(ww(3))
      if(idbuf.eq.0) then
        ierr = 3
        return
      endif

C--   Interpolation loop
      npt   = int(ww(8))
      nout  = min(npt,n)
      ibase = iqcG5ijk(stor7,1,1,idbuf)
      do i = 1,nout
        jj    = iaL(1,i)-1
        ia    = int(ww(jj+ 3)) + ibase
        ny    = int(ww(jj+ 4))
        nz    = int(ww(jj+ 5))
        ff(i) = dqcPdfPol(stor7,ia,ny,nz,ww(jj+6),ww(jj+12))
      enddo

C--   Done with the buffer
      call sqcReleaseScratch(idbuf)

      return
      end

C     ====================================================
      subroutine sqcTabIni(yy,ny,tt,nt,ww,ndim,nused,ierr)
C     ====================================================

C--   Initialise table-interpolation by setting up three arrays in ww:
C--
C--   ArrayY(9,ny) list of meshes and weights for all y-points
C--   ArrayT(9,nt) list of meshes and weights for all z-points
C--   ArrayF(3, m) list of gridpoints to be filled with data
C--
C--   yy,ny  (in) : y points (all inside grid or cuts)
C--   tt,nt  (in) : t points (all inside grid or cuts)
C--   ww    (out) : working array
C--   ndim   (in) : dimension of ww declared in the calling routine
C--   nused (out) : number of words used in ww
C--   ierr  (out) : 0 = OK
C--                 1 = Not enough space in ww (fatal error)
C--                 2 = No temporary buffer available
C--
C--   Layout of ww
C--   ------------
C--   ww(01)    =  Initialisation flag (654321)
C--   ww(02)    =  Evolution parameter key
C--   ww(03)    =  Scratch buffer id (0 = no buffer assigned)
C--   ww(04)    =  Number of points ArrayF
C--   ww(05)    =  K0 of ArrayF
C--   ww(06)    =  K1 of ArrayF
C--   ww(07)    =  K2 of ArrayF
C--   ww(08)    =  Number of y points ny
C--   ww(09)    =  K0 of ArrayY
C--   ww(10)    =  K1 of ArrayY
C--   ww(11)    =  K2 of ArrayY
C--   ww(12)    =  Number of t points nt
C--   ww(13)    =  K0 of ArrayT
C--   ww(14)    =  K1 of ArrayT
C--   ww(15)    =  K2 of ArrayT
C--
C--   ww(16)    =  First word of ArrayY
C--   ..           ..
C--   ww(nused) =  Last  word of ArrayF
C--
C--   Layout of ArrayY(i,j)
C--   ---------------------
C--   ArrayY(01,j) = iy of mesh of interpolation point jy
C--   ArrayY(02,j) = ia base address of (iy,1) in pdf table
C--   ArrayY(03,j) = ny mesh width = interpolation order in y
C--   ArrayY(04,j) = wy(1) first word of interpolation weights wy(6)
C--   ..
C--
C--   Layout of ArrayT(i,j)
C--   ---------------------
C--   ArrayT(01,j) = iz of mesh of interpolation point jt
C--   ArrayT(02,j) = ia relative address of (1,iz) in pdf table
C--   ArrayT(03,j) = nz mesh width = interpolation order in t
C--   ArrayT(04,j) = wz(1) first word of interpolation weights wz(6)
C--   ..
C--
C--   Layout of ArrayF(i,j)
C--   ---------------------
C--   ArrayF(01,j) = iy of data point j
C--   ArrayF(02,j) = iz of data point j
C--   ArrayF(03,j) = ia relative address in pdf table
C--
C--   The max size of ww is (with 9 mesh points / interpolation point):
C--
C--   nused = 15 + 9*(ny+nt) + 3*9*ny*nt

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension yy(*), tt(*), ww(*)
      dimension imin(2), imax(2), karrY(0:2), karrT(0:2), karrF(0:2)

C--   Address functions in ww
      iaY(i,j) = int(ww( 9)) + int(ww(10))*i + int(ww(11))*j
      iaT(i,j) = int(ww(13)) + int(ww(14))*i + int(ww(15))*j
      iaF(i,j) = int(ww( 5)) + int(ww( 6))*i + int(ww( 7))*j

C--   Initialize
      ierr = 0

C--   Book ArrayY in ww
      imin(1) = 1
      imax(1) = 9
      imin(2) = 1
      imax(2) = ny
      istart  = 16
      call smb_bkmat(imin,imax,karrY,2,istart,ilasty)

C--   Book ArrayT in ww
      imin(1) = 1
      imax(1) = 9
      imin(2) = 1
      imax(2) = nt
      call smb_bkmat(imin,imax,karrT,2,ilasty+1,ilastt)

C--   Book arrayF in ww (assuming 9 mesh points / interpolation point)
      imin(1) = 1
      imax(1) = 3
      imin(2) = 1
      imax(2) = 9*ny*nt
      call smb_bkmat(imin,imax,karrF,2,ilastt+1,nused)

C--   Check size of ww
      if(nused.gt.ndim) then
        ierr = 1
        return
      endif

C--   Store partition parameters in ww
      ww( 1) = 0.D0              !initialisation flag
      ww( 9) = dble(karrY(0))
      ww(10) = dble(karrY(1))
      ww(11) = dble(karrY(2))
      ww(13) = dble(karrT(0))
      ww(14) = dble(karrT(1))
      ww(15) = dble(karrT(2))
      ww( 5) = dble(karrF(0))
      ww( 6) = dble(karrF(1))
      ww( 7) = dble(karrF(2))

C--   Loop over the interpolation points in y
      do i = 1,ny
C--     Find interpolation mesh in y
        call sqcZmeshy(yy(i),iy1,iy2)
        my = iy2-iy1+1
C--     Store the mesh of each interpolation point
        jj        = iaY(1,i)-1
        ww(jj+ 1) = dble(iy1)
        ww(jj+ 2) = dble(iy1-1)
        ww(jj+ 3) = dble(my)
C--     Store interpolation weights
        call smb_polwgt(yy(i),ygrid2(iy1),my,ww(jj+4))
      enddo

C--   Loop over the interpolation points in t
      do i = 1,nt
C--     Find interpolation mesh in t
        call sqcZmesht(tt(i),0,iz1,iz2,it1)
        mz = iz2-iz1+1
C--     Store the mesh of each interpolation point
        jj        = iaT(1,i)-1
        ww(jj+ 1) = dble(iz1)
        ww(jj+ 2) = dble((iz1-1)*inciz7)
        ww(jj+ 3) = dble(mz)
C--     Store interpolation weights
        call smb_polwgt(tt(i),tgrid2(it1),mz,ww(jj+4))
      enddo

C--   Reserve a scratch buffer and calculate base address
      idbuf = iqcGimmeScratch()
      if(idbuf.eq.0) then
        ierr = 2
        return
      endif
      ibase = iqcG5ijk(stor7,1,1,idbuf)

C--   Setup list of mesh points
      nff = 0
      do it = 1,nt
        jj  = iaT(1,it)-1
        mz1 = int(ww(jj+1))
        mz  = int(ww(jj+3))
        mz2 = mz+mz1-1
        do iy = 1,ny
          kk  = iaY(1,iy)-1
          my1 = int(ww(kk+1))
          my  = int(ww(kk+3))
          my2 = my+my1-1
          maz = int(ww(jj+2))+ibase+my1-1-inciz7
          do mz = mz1,mz2
            maz = maz+inciz7
            may = maz-1
            do my = my1,my2
              may = may+1
              if(int(stor7(may)).ne.1) then !point not yet stored
                stor7(may) = dble(1)
                nff        = nff+1
                kk         = iaF(1,nff)-1
                ww(kk+1)   = dble(my)
                ww(kk+2)   = dble(mz)
                ww(kk+3)   = dble(may-ibase)
              endif
            enddo
          enddo
        enddo
      enddo

C--   Flag ww initialized and store number of points
      ww(1)  = dble(654321)   !initflag
      ww(2)  = dble(ipver6)   !parameter key
      ww(3)  = dble(0)        !no buffer assigned
      ww(4)  = dble(nff)
      ww(8)  = dble(ny)
      ww(12) = dble(nt)

C--   Release scratch buffer
      call sqcReleaseScratch(idbuf)

      return
      end

C     ================================
      subroutine sqcTabFun(ww,ff,ierr)
C     ================================

C--   Interpolate list of function values
C--
C--   ww     (in) : workspace previously filled by sqcLstIni
C--   ff    (out) : table of interpolated function values
C--   ierr  (out) : 0 = all OK
C--                 1 = ww not initialized
C--                 2 = ww not set-up with current parameters
C--                 3 = no buffer with function values available

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension ww(*), ff(*)

C--   Inline address function for 2-dim array
      iaff(i,j,n) = i + n*(j-1)

C--   Address functions in ww
      iaY(i,j) = int(ww( 9)) + int(ww(10))*i + int(ww(11))*j
      iaT(i,j) = int(ww(13)) + int(ww(14))*i + int(ww(15))*j

      ierr = 0

C--   Check ww initialized
      if(int(ww(1)).ne.654321) then
        ierr = 1
        return
      endif
C--   Check for evolution parameter change
      if(int(ww(2)).ne.ipver6) then
        ierr = 2
        return
      endif
C--   Check if scratch buffer assigned
      idbuf = int(ww(3))
      if(idbuf.eq.0) then
        ierr = 3
        return
      endif

C--   Interpolation loop
      npy   = int(ww( 8))
      npz   = int(ww(12))
      ibase = iqcG5ijk(stor7,1,1,idbuf)
      do iz = 1,npz
        jj  = iaT(1,iz)-1
        jaz = int(ww(jj+2))
        mz  = int(ww(jj+3))
        do iy = 1,npy
          kk        = iaY(1,iy)-1
          jay       = int(ww(kk+2))
          my        = int(ww(kk+3))
          ia        = ibase+jaz+jay
          iaf       = iaff(iy,iz,npy)
          ff(iaf)   = dqcPdfPol(stor7,ia,my,mz,ww(kk+4),ww(jj+4))
        enddo
      enddo

C--   Done with the buffer
      call sqcReleaseScratch(idbuf)

      return
      end

C=======================================================================
C==== Precooked functions to interpolate  ==============================
C=======================================================================

C     =======================================================
      double precision function
     +     dqcPdfSum(w,par,np,ix,iq,nf,ithr,ia,ff,first,jerr)
C     =======================================================

C--   Weighted sum of pdfs
C--
C--   w       (in) : workspace with pdf tables
C--   par     (in) : array of input parameters
C--   np      (in) : number of words in par              (not used)
C--   ix,iq   (in) : grid point                          (not used)
C--   nf      (in) : number of flavours
C--   ithr    (in) : threshold indicator                 (not used)
C--   ia      (in) : relative address of (ix,iq)
C--   ff      (in) : current buffer value                (not used)
C--   first   (in) : set to .true. at loop entry         (not used)
C--   jerr   (out) : error code at loop entry (0 = OK)

C--   Layout of par
C--   -------------
C--   par(1)  =  Karr(0)
C--   par(2)  =  Karr(1)
C--   par(3)  =  Karr(2)
C--   par(4)  =  Karr(3)
C--   par(5)  =  # pdfs for nf = 3
C--   par(6)  =  # pdfs for nf = 4
C--   par(7)  =  # pdfs for nf = 5
C--   par(8)  =  # pdfs for nf = 6
C--
C--   par(9)  =  first word of Parr(i,j,k)
C--   ..          ..
C--
C--   Layout of Parr(2,13+m,3:6)
C--   --------------------------
C--   Parr(1,j,nf) = base address of pdf j at nf
C--   Parr(2,j,nf) = coefficient  of pdf j at nf

      implicit double precision (a-h,o-z)

      dimension w(*), par(*)

      logical first,Ldum

C--   Address function in par
      iaP(i,j,k) = int(par(1))+int(par(2))*i+int(par(3))*j+int(par(4))*k

C--   Avoid compiler warnings for unused arguments
      idum = np
      idum = ix
      idum = iq
      idum = ithr
      Ldum = first
      ddum = ff

      jerr = 0
      pdf  = 0.D0
      nid  = int(par(nf+2))  !nid might be zero
      do j = 1,nid
        i   = iaP(1,j,nf)
        ia0 = int(par(i))
        wgt = par(i+1)
        pdf = pdf + w(ia0+ia) * wgt
      enddo

      dqcPdfSum = pdf

      return
      end

C=======================================================================
C==== One interpolation point ==========================================
C=======================================================================

C     =============================================
      subroutine sqcIntWgt(iy1,ny,it1,nt,y,t,wy,wz)
C     =============================================

C--   Precalculate weights for polynomial interpolation on a pdf table
C--
C--   iy1  (in)  :  lower limit of interpolation mesh in y
C--   ny   (in)  :  mesh size = interpolation order in y
C--   it1  (in)  :  lower limit of interpolation mesh in t
C--   nt   (in)  :  mesh size = interpolation order in t (or z)
C--   y,t  (in)  :  interpolation point
C--   wy   (out) :  weights for interpolation in y
C--   wz   (out) :  weights for interpolation in z
C--
C--   NB: y,t should be within the interpolation mesh, otherwise
C--       we have extrapolation. The routine sqcZmesh figures out
C--       what the interpolation mesh boundaries are for given y,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dimension wy(*),wz(*)

C--   Calculate interpolation weights
      call smb_polwgt(y,ygrid2(iy1),ny,wy)
      call smb_polwgt(t,tgrid2(it1),nt,wz)

      return
      end

C     ======================================================
      double precision function dqcPdfPol(w,ia0,ny,nz,wy,wz)
C     ======================================================

C--   One polynomial interpolation on a pdf table
C--
C--   w    (in)  :  workspace
C--   ia0  (in)  :  address of first mesh point (iy1,iz1,id)
C--   ny   (in)  :  interpolation order in y = mesh size in y
C--   nz   (in)  :  interpolation order in z = mesh size in z
C--   wy   (in)  :  interpolation weights in y
C--   wz   (in)  :  interpolation weights in z
C--
C--   NB: y,t should be within the interpolation mesh, otherwise
C--       we have extrapolation. The routine sqcZmesh figures out
C--       what the interpolation mesh boundaries are for given y,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      dimension w(*),wy(6),wz(6),work(3)

C--   Base address
      iadr = ia0-inciz7
C--   Do nz interpolations in y
      j  = 0
      do iz = 1,nz
        j       = j+1
        iadr    = iadr+inciz7
        work(j) = dmb_polin1(wy,w(iadr),ny)
      enddo
C--   Do one interpolation in t
      dqcPdfPol = dmb_polin1(wz,work,nz)
      
      return
      end

C=======================================================================
C==== Mesh routines ====================================================
C=======================================================================
      
C     ========================================================
      subroutine sqcZmesh(y,t,margin,iymi,iyma,izmi,izma,itmi)
C     ========================================================

C--   Find an interpolation mesh in (iy,iz) around y and t.
C--   By default the size of the mesh is ny = iosp and nt = 3 unless
C--   y or t are at a grid point in which case ny or nt are set to 1.
C--
C--   y,t          (in) : interpolation point
C--   margin       (in) : 0 or 1 points to stay away from a threshold
C--   iymi,iyma   (out) : limits of the interpolation mesh in y
C--   izmi,izma   (out) : limits of the interpolation mesh in z
C--   itmi        (out) : it-value corresponding to izmi

      implicit double precision (a-h,o-z)

      call sqcZmeshy(y,iymi,iyma)
      call sqcZmesht(t,margin,izmi,izma,itmi)

      return
      end

C     =================================
      subroutine sqcZmeshy(y,iymi,iyma)
C     =================================

C--   Find the boundaries of an interpolation mesh in iy around y.
C--   This is not entirely trivial because the mesh should not cross
C--   the grid boundary. The interpolation point y must be
C--   inside the boundaries of the qcdnum y grid.
C--
C--   By default the size of the mesh is ny = iosp unless y is at a grid
C--   point in which case ny is set to 1.
C--
C--   y            (in) : interpolation point
C--   iymi,iyma   (out) : limits of the interpolation mesh in y

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

C--   Get bin in y
      iy = iqcFindIy(y)
      if(iy.eq.-1) stop 'sqcZmeshy: y out of range ---> STOP'

C--   Check if on grid point
      if(iqcYhitIy(y,iy).eq.1) then
C--     On grid point iy
        iyma = iy
        iymi = iy
      else
C--     Mesh width is equal to the spline order ioy2
        iyma = min(iy+ioy2-1,nyy2(0))
        iymi = max(iyma-ioy2+1,0)
      endif

      return
      end

C     =============================================
      subroutine sqcZmesht(t,margin,izmi,izma,itmi)
C     =============================================

C--   Find the boundaries of an interpolation mesh in iz around t.
C--   This is not entirely trivial because the mesh
C--   should not  cross the grid boundaries or thresholds. It may
C--   happen that the mesh is too wide to fit between two thresholds or
C--   between a threshold and a boundary in which case we must adjust
C--   the width of the mesh. The interpolation point t must be
C--   inside the boundaries of the qcdnum t grid.
C--   The margin parameter defines the distance between the upper edge
C--   of the mesh and a threshold. Margin = 0 the upper edge can
C--   touch the threshold, margin = 1 the edge stays one point away
C--   from the threshold. Since the thresholds are at least two points
C--   apart, there should be no problem with margin = 1.
C--
C--   By default the size of the mesh is nt = 3 unless t is at a grid
C--   point in which case nt is set to 1.
C--
C--   t            (in) : interpolation point
C--   margin       (in) : 0 or 1 points to stay away from a threshold
C--   izmi,izma   (out) : limits of the interpolation mesh in z
C--   itmi        (out) : it-value corresponding to izmi

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      if(margin.ne.0 .and. margin.ne.1)
     +   stop  'sqcZmesht: invalid margin'

C--   Get bin in t
      it = iqcItfrmt(t)
      if(it.eq.0)  stop 'sqcZmesht: t out of range ---> STOP'
      iz  = izfit5( it)
      nf  = itfiz5(-iz)
      iz1 = iz15(nf)
      iz2 = iz25(nf)

C--   Check if on grid point
      if(iqcThitIt(t,it).eq.1) then
C--     On grid point it
        izma = iz
        izmi = iz
      else
C--     Mesh width is 3 for interpolation in t
        izma = min(iz+2,iz2-margin)
        izmi = max(izma-2,iz1)
        if(izma.le.izmi)
     +    stop 'sqcZmesht: zero or negative mesh width in t ---> STOP'
      endif

      itmi = itfiz5(izmi)
      itma = itfiz5(izma)
      if(itma-itmi .ne. izma-izmi)
     +    stop 'sqcZmesht: problem with mesh width in t ---> STOP'

      return
      end

C     =================================================================
      subroutine sqcMarkyt(mark,ylst,tlst,margin,iy1,iy2,iz1,iz2,it1,n)
C     =================================================================

C--   Process list of interpolation points and mark gridpoints
C--
C--   mark    (out)  : logical table containing the marks in y and z
C--   ylst(n) (in)   : list of y points
C--   tlst(n) (in)   : list of t points
C--   margin  (in)   : distance away from thresholds [0,1]
C--   iy1(n)  (out)  : list of lower mesh limits in y
C--   iy2(n)  (out)  : list of upper mesh limits in y
C--   iz1(n)  (out)  : list of lower mesh limits in z
C--   iz2(n)  (out)  : list of upper mesh limits in z
C--   it1(n)  (out)  : it-value corresponding to iz1
C--   n       (in)   : number of points in the list

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension ylst(*),tlst(*),iy1(*),iy2(*),iz1(*),iz2(*),it1(*)
      logical   mark(0:mxx0,0:mqq0+7)

C--   Initialize
      do j = 0,mqq0+7
        do i = 0,mxx0
          mark(i,j) = .false.
        enddo
      enddo
C--   Loop over interpolation points
      do i = 1,n
        call sqcZmesh(ylst(i),tlst(i),margin,iya,iyb,iza,izb,ita)
        iy1(i) = iya
        iy2(i) = iyb
        iz1(i) = iza
        iz2(i) = izb
        it1(i) = ita
        do iz = iza,izb
          do iy = iya,iyb
            mark(iy,iz) = .true.
          enddo
        enddo
      enddo

      return 
      end

C     ==========================================
      subroutine sqcMarkyy(marky,ylst,iy1,iy2,n)
C     ==========================================

C--   Process list of interpolation points and mark gridpoints in y
C--
C--   marky   (out)  : logical table containing the marks in y
C--   ylst(n) (in)   : list of y points
C--   iy1(n)  (out)  : list of lower mesh limits in y
C--   iy2(n)  (out)  : list of upper mesh limits in y
C--   n       (in)   : number of points in the list

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension ylst(*),iy1(*),iy2(*)
      logical   marky(0:mxx0)

C--   Initialize
      do i = 0,mxx0
        marky(i) = .false.
      enddo

C--   Loop over interpolation points
      do i = 1,n
        call sqcZmeshy(ylst(i),iy1(i),iy2(i))
        do iy = iy1(i),iy2(i)
          marky(iy) = .true.
        enddo
      enddo

      return 
      end

C     =====================================================
      subroutine sqcMarktt(markz,tlst,margin,iz1,iz2,it1,n)
C     =====================================================

C--   Process list of interpolation points and mark gridpoints in t
C--
C--   markz   (out)  : logical table containing the marks in z
C--   tlst(n) (in)   : list of t points
C--   margin  (in)   : distance away from thresholds [0,1]
C--   iz1(n)  (out)  : list of lower mesh limits in z
C--   iz2(n)  (out)  : list of upper mesh limits in z
C--   it1(n)  (out)  : it-value corresponding to iz1
C--   n       (in)   : number of points in the list

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension tlst(*),iz1(*),iz2(*),it1(*)
      logical   markz(0:mqq0+7)

C--   Initialize
      do i = 0,mqq0+7
        markz(i) = .false.
      enddo
C--   Loop over interpolation points
      do i = 1,n
        call sqcZmesht(tlst(i),margin,iz1(i),iz2(i),it1(i))
        do iz = iz1(i),iz2(i)
          markz(iz) = .true.
        enddo
      enddo

      return 
      end

C=======================================================================
C==== Repository =======================================================
C=======================================================================

C     ====================================
      subroutine sqcPolint(xa,ya,n,x,y,dy)
C     ====================================

C--   Not used anymore

C--   Polynomial tru n points using Neville's algorithm
C--   From Numerical Recipes, chapter 3.1
C--
C--   Input:  xa    table of n abscissa
C--           ya    table of n function values
C--           n     number of (x,y) points = order of polynomial
C--           x     interpolation point
C--   Output: y     interpolated polynomial P(x)
C--           dy    estimated error on y 

      implicit double precision(a-h,o-z)

      dimension xa(*),ya(*),c(10),d(10)

      if(n.gt.10) stop 'sqcPolint: degree n too large --> STOP'
      
      if(n.eq.2) then
        t  = (x-xa(1))/(xa(2)-xa(1))
        y  = (1-t)*ya(1) + t*ya(2)
        dy = 0.D0
        return
      endif  

      ns  = 1
      dif = abs(x-xa(1))
      do i = 1,n
        dift = abs(x-xa(i))
        if(dift.lt.dif) then
          ns  = i
          dif = dift
        endif
        c(i) = ya(i)
        d(i) = ya(i)
      enddo
      y = ya(ns)
      ns = ns-1
      do m = 1,n-1
        do i = 1,n-m
          ho   = xa(i)-x
          hp   = xa(i+m)-x
          w    = c(i+1)-d(i)
          den  = ho-hp
          if(den.eq.0) stop 'sqcPolint: equal abscissa --> STOP'
          den  = w/den
          d(i) = hp*den
          c(i) = ho*den
        enddo
        if(2*ns.lt.n-m) then
          dy = c(ns+1)
        else
          dy = d(ns)
          ns = ns-1
        endif
        y = y+dy
      enddo

      return
      end

C     ==========================================
      subroutine sqcPolin2(xa,nx,ya,ny,za,x,y,z)
C     ==========================================

C--   Not used anymore

C--   Two-dim polynomial interpolation using sqcPolint
C--
C--   xa, nx  (in)    table of nx x-abscissa
C--   ya, ny  (in)    table of ny y-abscissa  
C--   za      (in)    2-dim array za(nx,ny) with function values
C--   x,y     (in)    interpolation point
C--   z       (out)   interpolated polynomial P(x,y)

      implicit double precision(a-h,o-z)

      dimension xa(*),ya(*),za(nx,ny),work(10)

C--   Do ny interpolations in x
      do i = 1,ny
        call sqcPolint(xa,za(1,i),nx,x,work(i),epsi)
      enddo
C--   Do one interpolation in y
      call sqcPolint(ya,work,ny,y,z,epsi)

      return
      end

C     =======================================
      subroutine sqcGetABC(u,v,w,delta,a,b,c)
C     =======================================

C--   Not used anymore

C--   Caculate coefficients for quadratic interpolation in a bin of width
C--   delta. Let the interpolation formula be written as
C--
C--                    f(x) = A(x-x0)^2 + B(x-x0) + C
C--
C--   with x0 the lower edge of the bin and x0+delta the upper edge. Then
C--   this routine calculates the coefficients A, B and C, given delta and
C--   three sample points U = f(x0), V = f(x0+delta/2) and W = f(x0+delta).

      implicit double precision (a-h,o-z)

      a = ( 2*u - 4*v + 2*w)/(delta*delta)
      b = (-3*u + 4*v -   w)/ delta
      c =     u

      return
      end

