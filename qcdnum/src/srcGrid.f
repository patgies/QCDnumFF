
C--  x and Q2 grid routines on the file srcGrid.f
C--  --------------------------------------------
C--
C--   subroutine sqcGryDef(yi,iw,ng,nt,ntot,io)
C--   subroutine sqcGryMat(io)
C--   subroutine sqcGyMake(yy,iw,ng,nt,ntot,io)
C--
C--   double precision function dqcYfrmIy(iy,del)      
C--   integer function iqcIyfrmY(y,del,ny)
C--
C--   integer function iqcFindIg(y)
C--   integer function iqcFindIy(y)      
C--   integer function iqcYhitIy(y,iy)
C--
C--   double precision function dqcYjDiv(j,idiv,jmax)
C--   
C--   subroutine sqcGrTdef(tt,ww,nw,nt,ierr)
C--   subroutine sqcGtMake(ti,wi,nn,tt,nt,nd,ierr)
C--
C--   double precision function dqcTfrmIt(it)
C--   integer function iqcItfrmT(t)
C--   integer function iqcThitIt(t,it)
C--
C--   subroutine sqcSetIzCut
C--   logical function lqcInGrid(i,j)
C--   logical function lqcInside(x,qmu2)
C--   logical function lqcInsidex(x)
C--   logical function lqcInsideq(qmu2)
C--
C--   integer function iqcIyMaxG(iymax0,ig)

C===================================================================
C==  y-Grid routines ===============================================
C===================================================================

C     =========================================
      subroutine sqcGryDef(yi,iw,ng,nt,ntot,io)
C     =========================================

C--   Define equidistant logarithmic y-grid and setup Bsplines.
C--
C--   yi(i)  (in)    upper limits in y
C--   iw(i)  (in)    relative point densities
C--   ng     (in)    number of subgrids (# entries in yi and iw)
C--   nt     (in)    requested # of grid points not counting y = 0
C--   ntot   (out)   generated # of grid points not counting y = 0
C--   io     (in)    order of the Bspline interpolation
C--   
C--   Output in qgrid2.inc
C--
C--   ygrid2(0:nyy2) interpolation ygrid
C--   ymin2(0:nyg2)  lower limit of subgrids (1 - nyg2)
C--                  ymax2(0) lower limit of interpolation grid
C--   ymax2(0:nyg2)  upper limit of subgrids (1 - nyg2)
C--                  ymax2(0) upper limit of interpolation grid
C--   dely2(0:ngy2)  grid spacing for subgrids (1 - nyg2)
C--                  dely2(0) not defined
C--   nyy2(0:ngy2)   # grid points not counting y=0 for subgrids (1 - nyg2)
C--                  nyy2(0) # grid points not counting y=0 of 
C--                  interpolation grid = ntot above
C--   ioy2           spline order in y
C--   smaty2         transformation matrix h = S*a     on sub-grid
C--   sinvy2         transformation matrix a = Sinv*h  on sub-grid
C--   nmaty2         bandwidth of S
C--   ninvy2         bandwidth of Sinv

C--   NB: y-bin numbering is  0, 1, 2, ..., nyy2 with y0 = 0

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension yi(*),iw(*)
      dimension zz(mxx0),mz(mxx0)

C--   Check dimension (keep 10 words safety margin)
      if(nt.lt.2)         stop 'sqcGryDef: nt too small ---> STOP'
      if((nt+10).gt.mxx0) stop 'sqcGryDef: nt too large ---> STOP'
C--   Check number of subgrids
      if(ng.le.0 .or. ng.gt.mxg0) 
     +      stop 'sqcGryDef: invalid number of y-grids ---> STOP' 

C--   Generate grid, subgrids and fill (part of) common block qgrid2.inc
      call sqcGyMake(yi,iw,ng,nt,ntot,io)
C--   Copy to common block qgrid2.inc
      ioy2 = io
C--   Get from common block qgrid2.inc
      nyy  = nyy2(0)
C--   Temporary z-grid with 4 points extra to use for spline definition
C--   Bin numbering from 1 - nyy2+5 instead of 0 - nyy2+4 (extra
C--   points at upper end serve to push differently shaped end-point 
C--   Bsplines beyond the range of the y-grid)
      nztmp = nyy+5
      do i  = 1,nztmp
        zz(i) = i-1
        mz(i) = 1
      enddo
C--   Setup spline definition in y, ns is the number of generated splines
      do jo = 2,io
        call sqcSpyIni(jo,zz,mz,nztmp,ns,nc)
      enddo
C--   Setup spline transformation matrix S
      do jo = 2,3
        call sqcGryMat(jo)
      enddo
C--   Y-grid available
      Lygrid2 = .true.      
      ascut2  = 0.D0
      iymac2  = nyy2(0)
      
      return
      end

C     ========================
      subroutine sqcGryMat(io)
C     ========================

C--   Hardwire the transformation matrix S and its inverse
C--   Store as lower diagonal band matrix with equal values on the diagonals

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      if(io.eq.2) then
        smaty2(1,io) = 1.D0
        nmaty2(io)   = 1
        sinvy2(1,io) = 1.D0
        ninvy2(io)   = 1
      elseif(io.eq.3) then
        smaty2(1,io) = 0.5D0
        smaty2(2,io) = 0.5D0
        nmaty2(io)   = 2
        isign        = 1
        do i = 1,nyy2(0)
          sinvy2(i,io) = isign*2.D0
          isign        = -isign
        enddo
        ninvy2(io) = nyy2(0)
      else
        stop 'sqcGryMat: invalid spline order ---> STOP'
      endif

      return
      end

C     =========================================
      subroutine sqcGyMake(yy,iw,ng,nt,ntot,io)
C     =========================================

C--   Setup y-grid (G0) and equidistant sub-grids (G1,G2,G3 in this example)
C--
C--   yy(1--ng)  (in) Requested subgrid upper limits e.g. y(1),y(2),y(3)
C--   iw(1--ng)  (in) Relative point densities       e.g.   4 ,  2 ,  1
C--   ng         (in) Requested number of sub-grids  e.g.   3
C--   nt         (in) Requested number of gridpoints e.g.  17
C--   ntot      (out) Generated number of gridpoints e.g.  17
C--   io         (in) Interpolation order
C-- 
C--   y0              y1                      y2                      y3
C--    0 1 2 3 4 5 6 7 8   9  10  11  12  13  14      15      16      17
C--    *-*-*-*-*-*-*-*-*---*---*---*---*---*---*-------*-------*-------*
C--    
C--    0 1 2 3 4 5 6 7 8                                 
C--    *-*-*-*-*-*-*-*-*
C--      1 2 3 4 5 6 7 8   iy
C--      1 2 3 4 5 6 7 8   ia(lin)
C--      1 2 3 4 5 6 7 8   ia(quad)
C--                                             
C--    0   1   2   3   4   5   6   7   8   9  10
C--    *---*---*---*---*---*---*---*---*---*---*         
C--                        9  10  11  12  13  14  iy
C--                    9  10  11  12  13  14  15  ia(lin)
C--                9  10  11  12  13  14  15  16  ia(quad)       
C--
C--    0       1       2       3       4       5       6       7       8
C--    *-------*-------*-------*-------*-------*-------*-------*-------*
C--                                                   15      16      17
C--                                           16      17      18      19
C--                                   17      18      19      20      21
C--
C--   Output in common block qgrid2.inc
C--
C--   ygrid2(iy) interpolation ygrid
C--   nyg2       number of sub-grids                 0     1     2     3
C--   ymin2(ig)  lower subgrid limit                y0    y0    y1    y2
C--   ymax2(ig)  upper subgrid limit                y3    y1    y2    y3
C--   dely2(ig)  grid spacing                        -    d1    d2    d3
C--   iwgt2(ig)  relative point density              -     4     2     1
C--   iyma2(ig)  G0 index of ymax2                   -     8    14    17
C--   nyy2(ig)   highest grid point                 17     8    10     8
C--   iaoff2(ig) ia = iy(ig) + iaoff2(ig) iord = 2   -     0     5    11
C--                                       iord = 3   -     0     6    13
C--   nstory2    # of storage words       iord = 2  19
C--                                       iord = 3  21
  
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dimension yy(*), iw(*)
      dimension limy1(mxg0),limy2(mxg0)

C--   Initialise

     
C--   Check weights are integer multiples
      do i = ng-1,1,-1
        if(iw(i+1).le.0) stop 
     +     'sqcGyMake: zero or negative weight encountered ---> STOP'
        irat = iw(i)/iw(i+1)
        if(iw(i).ne.irat*iw(i+1)) stop 
     +     'sqcGyMake: weights not integer multiples ---> STOP'
      enddo
C--   Check yy in ascending order
      ylast = 0.D0
      do i = 1,ng
        if(yy(i).le.ylast) stop 
     +     'sqcGyMake: ygrid not in ascending order ---> STOP'
        ylast = yy(i)
      enddo
C--   Check requested number of grid points take safetymargin of 10
      if(nt.lt.ng) stop 
     +     'sqcGyMake: too little grid points requested (nt) ---> STOP'
      if((nt+10).gt.mxx0) stop 
     +     'sqcGyMake: too many grid points requested (nt) ---> STOP'
C--   Average grid spacing
      del  = 0.D0
      yim1 = 0.D0
      do i = 1,ng
        del  = del + iw(i)*(yy(i)-yim1)
        yim1 = yy(i)
      enddo
      del = del/nt
C--   Adjust del
      nyy2(ng) = int(iw(ng)*yy(ng)/del + 0.5D0)
      del      = iw(ng)*yy(ng)/nyy2(ng)
C--   Grid spacing of each grid
      do i = 1,ng
        dely2(i) = del/iw(i)
      enddo
C--   Adjust grid boundaries
      ymax2(ng) = yy(ng)
      do i = ng-1,1,-1
        iy       = int(yy(i)/dely2(i+1))
        nyy2(i)  = iy*iw(i)/iw(i+1)
        ymax2(i) = nyy2(i)*dely2(i)
      enddo
C--   Number of grid points
      ntot     = nyy2(1)
      limy1(1) = 1
      limy2(1) = ntot
      do i = 2,ng
        limy1(i) = nyy2(i-1)*iw(i)/iw(i-1) + 1
        limy2(i) = nyy2(i)
        ntot     = ntot + limy2(i) - limy1(i) + 1
      enddo
C--   Check number of grid points
      if((ntot+10).gt.mxx0) stop 
     +      'sqcGyMake: too many grid points generated (ntot) ---> STOP'
C--   Generate grid weedout subgrids with equal boundaries
      iy    = 0
      ig    = 0
      do i = 1,ng
        deli = dely2(i)
C--     If limy1 > limy2 we have to reject the subgrid
        ihit = 0
        do npt = limy1(i),limy2(i)
          iy         = iy+1
          ygrid2(iy) = npt*deli
          ihit       = 1
        enddo
C--     Weedout invalid subgrids
        if(ihit.eq.1) then 
          ig = ig+1
          iyma2(ig)  = iy
          dely2(ig)  = dely2(i)
          ymax2(ig)  = ymax2(i)
          nyy2(ig)   = nyy2(i)
          iwgt2(ig)  = iw(i)
        endif
      enddo
      if(iy.ne.ntot) stop 
     +      'sqcGyMake: error generating number of gridpoints ---> STOP'
C--   Fill common block
      nyg2      = ig            !number of subgrids
      nyy2(0)   = ntot          !number of points in interpolation grid
      dely2(0)  = 0.D0          !grid spacing
      ymax2(0)  = ygrid2(ntot)  !upper limit of interpolation grid
      ygrid2(0) = 0.D0          !lower limit of interpolation grid
      ymin2(0)  = 0.D0          !lower limit of interpolation grid
      ymin2(1)  = 0.D0          !lower limit of first subgrid
      do i = 2,nyg2             !lower limit of remaining subgrids
        ymin2(i) = ymax2(i-1)
      enddo
C--   Address offsets
      iaoff2(1) = 0
      do ig = 2,nyg2
        iymi = iqcIyfrmY(ymax2(ig-1),dely2(ig),nyy2(ig))
        iaoff2(ig) = nyy2(ig-1)+iaoff2(ig-1)-iymi+io-1
      enddo
      nstory2 = nyy2(nyg2)+iaoff2(nyg2)
      
C--   NEW CODE -------------------------------------------------------

C--   y0              y1                      y2                      y3
C--    0 1 2 3 4 5 6 7 8   9  10  11  12  13  14      15      16      17
C--    *-*-*-*-*-*-*-*-*---*---*---*---*---*---*-------*-------*-------*
C--    
C--    0 1 2 3 4 5 6 7 8                                 
C--    *-*-*-*-*-*-*-*-*
C--    0 1 2 3 4 5 6 7 8                iy0fiyg2(i=0-nyy2(ig),ig=1-nyg2)
C--                                             
C--    0   1   2   3   4   5   6   7   8   9  10
C--    *---*---*---*---*---*---*---*---*---*---*         
C--    0   2   4   6   8   9  10  11  12  13  14
C--
C--    0       1       2       3       4       5       6       7       8
C--    *-------*-------*-------*-------*-------*-------*-------*-------*
C--    0       4       8      10      12      14      15      16      17
C--    
C--
C--   ygrid2(iy) interpolation ygrid
C--   nyg2       number of sub-grids                 0     1     2     3
C--   ymin2(ig)  lower subgrid limit                y0    y0    y1    y2
C--   ymax2(ig)  upper subgrid limit                y3    y1    y2    y3
C--   iyma2(ig)  G0 index of ymax2                   -     8    14    17
C--   dely2(ig)  grid spacing                        -    d1    d2    d3
C--   iwgt2(ig)  relative point density              -     4     2     1
C--   jymi2(ig)  first point to copy to G0           -     1     5     6
C--   nyy2(ig)   number of grid points              17     8    10     8

C--   Setup subgrid index tables
      do ig = 1,nyg2
        iy0fiyg2(0,ig) = 0
        do iy = 1,nyy2(ig)
          yi              = iy*dely2(ig)
          iy0fiyg2(iy,ig) = iqcFindIy(yi)
        enddo        
      enddo
C--   Set lowest subgrid index to copy
      jymi2(1) = 1
      do ig = 2,nyg2
        jymi2(ig) = iqcIyfrmY(ymin2(ig),dely2(ig),nyy2(ig))+1
      enddo

C--   Debug printout      
*      write(6,*) 'gridpoints ',(nyy2(ig),ig=1,nyg2)
*      write(6,*) 'firstcopy  ',(jymi2(ig),ig=1,nyg2)
*      write(6,*) 'iymax      ',(iyma2(ig),ig=1,nyg2)
*      do iy = 0,nyy2(0)
*         write(6,'(I5,E13.5,5I5)')
*     +   iy,ygrid2(iy),(iy0fiyg2(iy,ig),ig=1,nyg2)
*      enddo
      
C--   END OF NEW CODE ------------------------------------------------      

C--   Debug printout
*      do i = 1,nyg2
*        write(6,*) 'ig,nyy,ymi=',i,nyy2(i)
*      enddo
*      do i = 1,nyg2
*        write(6,*) 'Grid ',i
*        nn = min(nyy2(i),20)
*        write(6,'(20I3)') (iy,iy=nn,0,-1)
*      enddo
*      write(6,*) 'npt,nstor=',nyy2(0),nstory2

      return
      end

C=======================================================================
C==  Search in equidistant subgrid =====================================
C=======================================================================

C     ===========================================
      double precision function dqcYfrmIy(iy,del)
C     ===========================================

C--   Get y from binnumber iy in equidistant subgrid with spacing del.
C--   Warning: iy should run from 0-nyy2
C--
C--   iy  (in)  grid index in subgrid
C--   del (in)  gridspacing of subgrid

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dqcYfrmIy = iy*del

      return
      end

C     ====================================
      integer function iqcIyfrmY(y,del,ny)
C     ====================================

C--   Get binnumber iy from value of y in equidistant subgrid
C--   iy = -1 if y outside grid
C--
C--   y   (in) value of y
C--   del (in) gridspacing of subgrid
C--   ny  (in) total number of grid points not including y = 0 

      implicit double precision (a-h,o-z)
      logical lmb_eq

      include 'qcdnum.inc'
      include 'qpars6.inc'

      iy = int(y/del)
      if((iy.lt.0).or.(iy.gt.ny)) then
        iqcIyfrmY = -1
      elseif(iy.lt.ny .and. lmb_eq(y,(iy+1)*del,aepsi6)) then
        iqcIyfrmY = iy+1     !snap to upper edge
      else
        iqcIyfrmY = iy
      endif

      return
      end

C=======================================================================
C==  Search in y-Grid G0 ===============================================
C=======================================================================
    
C     =============================
      integer function iqcFindIg(y)
C     =============================

C--   Find sub-grid index
 
      implicit double precision (a-h,o-z)
      logical lmb_eq

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

C--   Check if y at lowest grid point
      if(lmb_eq(y,0.D0,aepsi6))                 then
        iqcfindig = -1
        return
      endif
C--   -1 if outside grid
      iqcfindig = -1
      do i = 1,nyg2
        if(lmb_eq(y,ymin2(i),aepsi6))           then
          iqcfindig = max(i-1,1)
          return
        elseif(lmb_eq(y,ymax2(i),aepsi6))       then
          iqcfindig = i
          return
        elseif(y.gt.ymin2(i) .and. y.le.ymax2(i)) then
          iqcfindig = i
          return
        endif
      enddo

      return
      end
    
C     =============================
      integer function iqcFindIy(y)
C     =============================

C--   Find bin number iy corresponding to y
 
      implicit double precision (a-h,o-z)
      logical lmb_eq

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

C--   Accept if y at highest grid point
      if(lmb_eq(y,ygrid2(nyy2(0)),aepsi6))    then
        iqcfindiy = nyy2(0)
        return
      endif
      
C--   Accept y = 0 
      if(lmb_eq(y,0.D0,aepsi6))    then
        iqcfindiy = 0
        return
      endif     
      
C--   -1 if outside grid
      if(y.le.0.D0 .or. y.gt.ygrid2(nyy2(0)) )  then
        iqcfindiy = -1 
        return
      endif

C--   Which sub-grid?
      igrid = iqcFindIg(y)

C--   Sub-grid not found
      if(igrid.eq.-1) stop 'iqcFindIy: cannot find subgrid ---> STOP'

      ylow = 0.D0
      ilow = 0
      if(igrid.gt.1) then 
        ylow = ymax2(igrid-1)
        ilow = iyma2(igrid-1)
      endif

C--   Make sure we dont slip below the lower boundary of the current grid
      if(igrid.gt.1 .and. lmb_eq(y,ymax2(igrid-1),aepsi6)) then
        iy = iyma2(igrid-1)
      else
C--     Add position w.r.t lower limit of current grid
        iy = int(ilow + (y-ylow)/dely2(igrid))
C--     Check for slip
        if(lmb_eq(y,ygrid2(iy+1),aepsi6)) iy = iy+1
      endif

      iqcfindiy = iy

      return
      end

C     ================================
      integer function iqcYhitIy(y,iy)
C     ================================

C--   Returns  1 if y is on (i.e. very close to) gridpoint iy.
C--   Returns -1 otherwise.
C--
C--   del = gridspacing
C--   ny  = total number of grid points not including y = 0 

      implicit double precision (a-h,o-z)
      logical lmb_eq

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      iqcYhitIy = -1
      if((iy.le.0).or.(iy.gt.nyy2(0)))  return
      if(lmb_eq(y,ygrid2(iy),aepsi6)) iqcYhitIy = 1

      return
      end
      
C=======================================================================
C==  Subdivide y-grid ==================================================
C=======================================================================

C--   This subdivided grid is nowhere used in QCDNUM itself but
C--   comes in handy in some of the author's testjobs

C     ===============================================
      double precision function dqcYjDiv(j,idiv,jmax)
C     ===============================================

C--   Sometimes we want a finer sampling than given by the y-grid 
C--   itself (e.g. for plotting). This is done by dividing each bin
C--   in idiv equal parts as shown below for a 2-bin y-grid
C--
C--   0-----------1-----------2    idiv = 1   (original y-grid)
C--   |           |           |
C--   0-----1-----2-----3-----4    idiv = 2   (jmax = 4)  
C--   |           |           |
C--   0---1---2---3---4---5---6    idiv = 3   (jmax = 6)
C--   |           |           |
C--   0--1--2--3--4--5--6--7--8    idiv = 4   (jmax = 8)
C--
C--   The function dqcYjDiv returns the y-value of a gridpoint in the
C--   in the sub-divided grid
C--
C--   Input:   j    [0-jmax] grid point index (see below for jmax)
C--            idiv [>= 1]   number of divisions in each y-grid bin
C--   Output:  jmax          number of grid points of the subdivided grid
C--            dqcYjDiv      y-value of grid point j
C--
C--   Remark:  (1) The y-grid (stored in qgrid2.inc) must have been
C--                defined before; the grid does not need to be
C--                equidistant for dqcYjDiv to work
C--            (2) The function abends when idiv < 1
C--            (3) dqcYjDiv returns zero in case j is not in range         

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

C--   Gotcha..
      if(idiv.le.0) stop 'dqcYjDiv: idiv .le. 0 ---> STOP'
C--   Number of grid points
      jmax = nyy2(0)*idiv
C--   Get y-value
      dqcYjDiv = 0.D0
      if(j.le.0 .or. j.gt.jmax) return
      iy       = (j-1)/idiv
      del      = ( ygrid2(iy+1)-ygrid2(iy) ) / idiv
      ii       = j - iy*idiv
      dqcYjDiv = ygrid2(iy) + ii*del

      return
      end

C=======================================================================
C==  Q2-Grid routines ==================================================
C=======================================================================

C     =============================================
      subroutine sqcGrTdef(tt,ww,nw,nt,loglog,ierr)
C     =============================================

C--   Setup grid in t = ln(Q^2).
C--
C--   tt     (in)    array with nw values of t in strictly ascending order
C--                  tt(1) and tt(nw) are the limits of the grid
C--   ww     (in)    array with nw-1 weight values. The point density in
C--                  the region tt(i) and tt(i+1) is proportional to w(i)
C--   nw     (in)    number of items in tt and ww (NB: ww(nw) is not used)
C--   nt     (inout) requested number of grid points (>=2) on input
C--                  generated number of grid points on exit
C--   loglog (in)    generate loglog grid, if true
C--   ierr   (out)   0  all OK
C--                  1  wrong nw
C--                  2  tt not in ascending order
C--                  3  zero or negative weights
C--                  4  not enough space in tgrid2
C--
C--   The output grid is stored in /qgrid2/tgrid2(mqq0) with ntt2 points

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      dimension tt(*),ww(*)
      
      logical loglog

C--   Copy tt to tgrid2 if requested (i.e. nt <= nw)
      if(nt.le.nw) then
        if(nw.le.mqq0) then
          do i = 1,nw
            tgrid2(i) = tt(i)
          enddo
          nt   = nw
          ntt2 = nw
          ierr = 0
          return
        else
          ierr = 1
          return
        endif
      endif
C--   So we dont want to copy ---> generate grid
C--   log or loglog, thats the question
      if(.not.loglog) then
        call sqcGtMake(tt,ww,nw,tgrid2,nt,mqq0,ierr)
      else
        do i = 1,nw
          tt(i) = log(tt(i))
        enddo
        call sqcGtMake(tt,ww,nw,tgrid2,nt,mqq0,ierr)
        do i = 1,nw
          tt(i) = exp(tt(i))
        enddo
        do i = 1,nt
          tgrid2(i) = exp(tgrid2(i))
        enddo
      endif          
C--   Error condition
      if(ierr.ne.0) return
C--   Copy number of grid points to common block
      ntt2 = nt
C--   T-grid available
      Ltgrid2 = .true.      
      ascut2  = 0.D0

      return
      end

C     ============================================
      subroutine sqcGtMake(ti,wi,nn,tt,nt,nd,ierr)
C     ============================================

C--   Generate grid: ti defines regions, wi the relative point density
C--   in each region and nt the total number of points requested for
C--   the output grid. The routine generates tt with approximately
C--   nt points in total, with point densities as specified by wi.
C--
C--   ti   (in)    list of points in ascending order
C--   wi   (in)    point densities in region (i), up to a common factor
C--   nn   (in)    number of points in ti and wi 
C--   tt   (out)   output grid
C--   nt   (inout) # gridpoints: as requested on input, as generated on exit
C--   nd   (in)    dimension of tt in the calling routine
C--   ierr (out)   0  all OK
C--                1  wrong nn or nd
C--                2  ti not in ascending order
C--                3  zero or negative weights
C--                4  not enough space in tt

      implicit double precision (a-h,o-z)

      dimension ti(*),wi(*),tt(*)

      ierr = 0
      if((nn.lt.2).or.(nd.lt.2)) then
        ierr = 1
        stop 'sqcGrMake: nn or nd lesser than 2 ---> STOP'
      endif

C--   Just copy input to output when nt <= nn
      if(nt.le.nn) then
        if(nt.gt.nd-1) then
          ierr = 4
          stop 'sqcGrMake: too many grid points requested ---> STOP'
        endif
        do i = 1,nn-1
          if(ti(i+1).le.ti(i)) then
            ierr = 2
            stop 'sqcGrMake: ti not in ascending order ---> STOP'
          endif
          tt(i) = ti(i)
        enddo
        tt(nn) = ti(nn)
        nt     = nn
        return
      endif

C--   Get common scaling factor
      sum = 0.D0
      do i = 1,nn-1
        if(ti(i+1).le.ti(i)) then
          ierr = 2
          stop 'sqcGtMake: ti not in ascending order ---> STOP'
        endif
        if(wi(i).le.0.D0) then
          ierr = 3 
          stop 'sqcGtMake: zero or negative weight ---> STOP'
        endif
        sum = sum + wi(i)*(ti(i+1)-ti(i))
      enddo
      fac = (nt-1)/sum

C--   Divide each region in np bins and generate grid points
      nt   = 0
      do i = 1,nn-1
        del = ti(i+1)-ti(i)
        np  = int(max(fac*del*wi(i)+0.5,1.D0))
C--     Generate at least one intermediate point 
        np  = max(np,2)
        bw  = del/np
        t0  = ti(i)
        do j = 1,np
          nt = nt+1
          if(nt.gt.nd-1) then
            ierr = 4
            stop 'sqcGtMake: too many grid points ---> STOP'
          endif
          tt(nt) = t0 + (j-1)*bw
        enddo
      enddo
C--   Upper end of the grid
      nt     = nt+1
      tt(nt) = ti(nn)

      return
      end

C     =======================================
      double precision function dqcTfrmIt(it)
C     =======================================

C--   Get t (= ln Q2) from binnumber it.
C--   Warning: there is no check that it runs from 1-ntt2.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dqcTfrmIt = tgrid2(it) 

      return
      end

C     =============================
      integer function iqcItfrmT(t)
C     =============================

C--   Get binnumber it from value of t (= ln Q2).
C--   Returns it = 0 if t is out of range.

      implicit double precision (a-h,o-z)
      logical lmb_eq,lmb_lt,lmb_gt,lmb_ge

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      save ilast
      data ilast /1/

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(ilast)

C--   t < t1
      if(lmb_lt(t,tgrid2(1),aepsi6))       then
        iqcItfrmt = 0
        ilast     = 1
        return
      endif
C--   t > tmax
      if(lmb_gt(t,tgrid2(ntt2),aepsi6))    then
        iqcItfrmt = 0
        ilast     = 1
      endif

C--   Goto binary search if t < tlast
      if(lmb_lt(t,tgrid2(ilast),aepsi6))   go to 10
C--   t >= tlast, now check t < t(last+1)
      if(lmb_lt(t,tgrid2(ilast+1),aepsi6)) then
        iqcItfrmt = ilast
        return
      endif
C--   Check if t at highest grid point
      if(lmb_eq(t,tgrid2(ntt2),aepsi6))    then
        iqcItfrmt = ntt2
        ilast     = ntt2-1
        return
      endif
       
C--   Binary search
   10 i = 1
      j = ntt2+1
   20 k = (i+j)/2
      if (lmb_lt(t,tgrid2(k),aepsi6)) j = k
      if (lmb_ge(t,tgrid2(k),aepsi6)) i = k
      if (j.gt.i+1)                     go to 20

      iqcItfrmt = i
      ilast     = i
      
      return
      end

C     ================================
      integer function iqcThitIt(t,it)
C     ================================

C--   Returns  1 if t is on (i.e. very close to) gridpoint it.
C--   Returns -1 otherwise. 

      implicit double precision (a-h,o-z)
      logical lmb_eq

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      iqcThitIt = -1
      if((it.le.0).or.(it.gt.ntt2))     return
      if(lmb_eq(t,tgrid2(it),aepsi6)) iqcThitIt = 1

      return

      end

C=======================================================================
C==  Cut managing routines =============================================
C=======================================================================

C     ======================
      subroutine sqcSetIzCut
C     ======================

C--   Put cuts on the z-index
C--   Needs a valid nfmap for izfitL2 and izfitU2 functions to work

      implicit double precision(a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      izmic2 = izfitL2(itmic2)
      izmac2 = izfitU2(itmac2)

      return
      end

C     ===============================
      logical function lqcInGrid(i,j)
C     ===============================

C--   Check if point (i,j) is inside grid

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      if(i.ge.1 .and. i.le.nyy2(0) .and. j.ge.1 .and. j.le.ntt2) then
        lqcInGrid = .true.
      else
        lqcInGrid = .false.
      endif

      return
      end

C     ==================================
      logical function lqcInside(x,qmu2)
C     ==================================

C--   Check if x and qmu2 are inside grid or cuts

      implicit double precision (a-h,o-z)
      logical lqcInsidex, lqcInsideq, Lx, Lq

      Lx = lqcInsidex(x)
      Lq = lqcInsideq(qmu2)
      
      lqcInside = Lx .and. Lq

      return
      end

C     ==============================
      logical function lqcInsidex(x)
C     ==============================

C--   Check if x is inside grid or cuts

      implicit double precision (a-h,o-z)
      logical lmb_le

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

      if(lmb_le(xmic5,x,-aepsi6) .and. lmb_le(x,1.D0,-aepsi6)) then
        lqcInsidex = .true.
      else
        lqcInsidex = .false.
      endif

      return
      end

C     ==============================
      logical function lqcInsideq(q)
C     ==============================

C--   Check if q is inside grid or cuts

      implicit double precision (a-h,o-z)
      logical lmb_le

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

      if(lmb_le(qmic5,q,-aepsi6) .and. lmb_le(q,qmac5,-aepsi6)) then
         lqcInsideq = .true.
      else
         lqcInsideq = .false.
      endif

      return
      end

C     =====================================
      integer function iqcIyMaxG(iymax0,ig)
C     =====================================

C--   Translate y-max index in G0 to an y-max index in Gi
C--
C--   iymax0     (in)    [1,nyy2(0)]    y-max index in G0 
C--   ig         (in)    [1,nyg2]       subgrid index 
C--   iqcIyMaxG  (out)   [1,nyy2(ig)]   y-max index in subgrid

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      
      ymax = ygrid2(iymax0)
      iy   = iqcIyfrmY(ymax,dely2(ig),nyy2(ig)) !find index below ymax  
      if(iy.eq.-1) then 
        iy = nyy2(ig)               !index out of subgrid range
      else
        iy = min(iy+1,nyy2(ig))     !put index above ymax
      endif
      
      iqcIyMaxG = iy
      
      return
      end


