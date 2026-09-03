
C--   This is the file srcEvolNN.f containing the nxn evolution routines

C--   subroutine sqcEvFillA(ww,id,fun)
C--   double precision function sqcEvGetAA(ww,id,it,nf,ithrs)
C--
C--   subroutine sqcEvDfill(w,idf,start,m,n,iz1,iz2)
C--   subroutine sqcEvDglap(w,idw,ida,idf,start,m,n,iz1,iz2,eps)
C--   subroutine sqcNNallg(
C--        ww,iww,ec,iec,ff,iff,idim,nopt,nf,iz1,iz2,nnn,eps)
C--   subroutine sqcNNupdn(
C--      ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,iz1,iz2,nnn,eps)
C--   subroutine sqcNNsubg(
C--        ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,nyg,iz1,iz2,nnn,eps)
C--
C--   subroutine sqcNNmult(w,iw,a,ia,b,ib,n,m,nbnd,idim)
C--   subroutine sqcNNeqs(w,iw,a,ia,b,ib,n,m,idim,ierr)
C--
C--   double precision function dqcNNGetEps(aa,i,ny)
C--   subroutine sqcNNgetVj(ww,id,iz,ig,nyg,buf)
C--   subroutine sqcNNputVj(ww,id,iz,ig,nyg,buf)
C--   subroutine sqcNNFjtoAj(ffin,aout,ny)
C--   subroutine sqcNNAjtoFj(aain,fout,ny)

C==   ==================================================================
C==   Perturbative coefficients ========================================
C==   ==================================================================

C     ================================
      subroutine sqcEvFillA(ww,id,fun)
C     ================================

C--   ww      (in)  store with type-6 tables
C--   id      (in)  global identifier of a type-6 table
C--   fun     (in)  function returning the coefficient value
C--
C--   The syntax of fun is:
C--
C--   double precision function fun(iq,nf,ithresh)
C--
C--   iq      (in)  mu2 gridpoint index (not t!)
C--   nf      (in)  number of lavours at iq
C--   ithresh (in)  0 = not at thr, +-1 at thr with larger/smaller nf

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      external fun

      dimension ww(*)

      nflast = itfiz5(-1)
      ia     = iqcG6ij(ww,1,id)-1
      do iz = 1,nzz2
        it      = itfiz5( iz)
        nf      = itfiz5(-iz)
        iq      = it
        ithresh = 0
        if(it.eq.itchm2 .or. it.eq.itbot2 .or. it.eq.ittop2) then
          if(nf.eq.nflast) then
            ithresh = -1
          elseif(nf.eq.nflast+1) then
            ithresh =  1
          else
            stop 'sqcEVFILLA inconsistent nf --> STOP'
          endif
        endif
        nflast    = nf
        ww(ia+iz) = fun(iq,nf,ithresh)
      enddo
C--   Set itmin and flag table as filled
      call sqcSetMin6(ww,id,1)
      call sqcValidate(ww,id)

      return
      end

C     =======================================================
      double precision function sqcEvGetAA(ww,id,it,nf,ithrs)
C     =======================================================

C--   ww      (in)  store with type-6 tables
C--   id      (in)  global identifier of a type-6 table
C--   it      (in)  t-grid index (same as iq)
C--   nf     (out)  number of flavours
C--   ithrs  (out)  threshold indicator
C--
C--   If it < 0 use nf = (3,4,5) instead of (4,5,6) at iqc,b,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension ww(*)

C--   Find iz index
      ithrs = 0
      iz    = izfit5( it)
      nf    = itfiz5(-iz)
      if(it.gt.0) then
        if( it.eq.itchm2 .or. it.eq.itbot2 .or. it.eq.ittop2) ithrs =  1
      elseif(it.lt.0) then
        if(-it.eq.itchm2 .or.-it.eq.itbot2 .or.-it.eq.ittop2) ithrs = -1
      else
        stop 'sqcEvGetAA encounter it = 0'
      endif
C--   Return alfa (and also nf and ithrs)
      ia         = iqcG6ij(ww,iz,id)
      sqcEvGetAA = ww(ia)

      return
      end

C==   ==================================================================
C==   Evolution routines ===============================================
C==   ==================================================================

C     ==============================================
      subroutine sqcEvDfill(w,idf,start,m,n,iz1,iz2)
C     ==============================================

C--   Same as sqcEvDglap but fill from iz1 to iz2 instead of evolve
C--   sqcEVPLAN should have been called before
C--
C--   w          (in)    Workspace
C--   idf(i)     (in)    Pdf identifiers
C--   start(i,j) (in)    either starting pdf_i(yj) or discontinuity at iz1
C--   m          (in)    first dim of start
C--   n          (in)    abs(n) = number of pdfs to fill simultaneously
C--                      n > 0  = select input mode
C--                      n < 0  = select transfer mode
C--   iz1,2    (inout)   filling from iz1 to iz2
C--
C--   n > 0 input    mode: start = input pdf on entry, untouched on exit
C--   n < 0 transfer mode: start = discontinuity on entry, zero on exit

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

      dimension w(*),idf(*)
      dimension start(m,*)

C--   Number of coupled evolutions
      nabs = abs(n)
      istep = 1
      if(iz2 .lt. iz1) istep = -1

C--   Set start values in iz = -ig of all pdfs
      do i = 1,nabs
        id = idf(i)
        do ig = 0,nyg2
          ia = iqcG5ijk(w,1,-ig,id)-1
          if(n.gt.0) then
            do iy = 1,nyy2(0)
              ia    = ia+1
              w(ia) = start(i,iy) !input mode: set start
            enddo
          else
            do iy = 1,nyy2(0)
              ia    = ia+1
              w(ia) = w(ia)+start(i,iy) !transfer mode: add discontinuity
            enddo
          endif
        enddo
      enddo

C--   Do the work
      do i = 1,nabs
        id  = idf(i)
        ia0 = iqcG5ijk(w,1,0,id)-1
        do iz = iz1,iz2,istep
          ia1 = iqcG5ijk(w,1,iz,id)-1
          do iy = 1,nyy2(0)
            w(ia1+iy) = w(ia0+iy) !copy start to iz
          enddo
        enddo
      enddo

C--   Transfer mode return zero in start array
      if(n.lt.0) then
        do i = 1,nabs
          do iy = 1,nyy2(0)
            start(i,iy) = 0.D0
          enddo
        enddo
      endif

      return
      end

C     ==========================================================
      subroutine sqcEvDglap(w,idw,ida,idf,start,m,n,iz1,iz2,eps)
C     ==========================================================

C--   N-fold evolution at fixed nf
C--   sqcEVPLAN should have been called before
C--
C--   w          (in)    Workspace
C--   idw(i,j,k) (in)    Pij weight table identifiers of order k
C--   ida(i,j,k) (in)    aij coefficient table identifiers of order k
C--   idf(i)     (in)    Pdf identifiers
C--   start(i,j) (in)    either starting pdf_i(yj) or discontinuity at iz1
C--   m          (in)    first 2 dims of idw and ida and first dim of start
C--   n          (in)    abs(n) = number of pdfs to evolve simultaneously
C--                      n > 0  = select input mode
C--                      n < 0  = select transfer mode
C--   iz1,2    (inout)   evolution from iz1 to iz2
C--   eps        (out)   max deviation at midpoint of xbins
C--
C--   n > 0 input    mode: start = input pdf on entry, evolved pdf on exit
C--   n < 0 transfer mode: start = discontinuity on entry, zero on exit

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension w(*),idw(m,m,*),ida(m,m,*),idf(*)
      dimension start(m,*)

C--   Number of coupled evolutions
      nabs = abs(n)

*mb   Print identifieres
*      write(6,'(/'' sqcEvDglap Eij'')')
*      write(6,'('' iord 1 Eij = '',4I10)')
*     + ((ida(i,j,1),j=1,nabs),i=1,nabs)
*      write(6,'('' iord 2 Eij = '',4I10)')
*     + ((ida(i,j,2),j=1,nabs),i=1,nabs)
*      write(6,'('' iord 3 Eij = '',4I10)')
*     + ((ida(i,j,3),j=1,nabs),i=1,nabs)
*      write(6,'(''        Fi  = '',2I10)')
*     + (idf(i),i=1,nabs)

C--   Number of perturbative terms
      nopt = iord6

C--   Number of flavours
      nf = itfiz5(-iz1)

C--   Set start values in iz = -ig of all pdfs
      do i = 1,nabs
        id = idf(i)
        do ig = 0,nyg2
          ia = iqcG5ijk(w,1,-ig,id)-1
          if(n.gt.0) then
            do iy = 1,nyy2(0)
              ia    = ia+1
              w(ia) = start(i,iy) !input mode: set start
            enddo
          else
            do iy = 1,nyy2(0)
              ia    = ia+1
              w(ia) = w(ia)+start(i,iy) !transfer mode: add discontinuity
            enddo
          endif
        enddo
      enddo

C--   Do the work
      call sqcNNallg(w,idw,w,ida,w,idf,m,nopt,nf,iz1,iz2,nabs,eps)

C--   Store endvalues in the start array
      do i = 1,nabs
        id = idf(i)
        ia = iqcG5ijk(w,1,iz2,id)-1
        if(n.gt.0) then
          do iy = 1,nyy2(0)
            ia          = ia+1
            start(i,iy) = w(ia) !input mode: return endvalue
          enddo
        else
          do iy = 1,nyy2(0)
            start(i,iy) = 0.D0  !transfer mode: return zero
          enddo
        endif

      enddo

      return
      end

C     =======================================================
      subroutine sqcNNallg(
     +     ww,iww,ec,iec,ff,iff,idim,nopt,nf,iz1,iz2,nnn,eps)
C     =======================================================

C--   Evolution from iz1 to iz2 with fixed nf over all subgrids
C--
C--   ww           (in)    weight table store
C--   iww(i,j,k)   (in)    weight table identifiers of Pij for each order k
C--   ec           (in)    expansion coefficient store
C--   iec(i,j,k)   (in)    coefficient table id for corresponding Pij
C--   ff           (in)    pdf store
C--   iff(i)       (in)    pdf table identifiers for i = 1,nnn
C--   idim         (in)    first 2 dims of iww and iec
C--   nopt         (in)    number of perturbative terms
C--   nf           (in)    number of flavours
C--   iz1          (in)    start of evolution
C--   iz2          (in)    end of evolution
C--   nnn          (in)    number of coupled equations
C--   eps          (out)   max |lin-quad| mid-between grid points
C--
C--   Input :  Start values stored at iz = -ig for all subgrids ig
C--   Output:  Evolved pdf values for iz1 to iz2
C--            End   values stored at iz = -ig for all subgrids ig

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension ww(*), ec(*), ff(*)
      dimension iww(idim,idim,*), iec(idim,idim,*), iff(*)

      eps = 0.D0

C--   Loop over subgrids (from coarse to dense)
      do ig = nyg2,1,-1
C--     Evolve
        call sqcNNupdn(
     +  ww,iww,ec,iec,idim,ff,iff,ig,nopt,nf,iz1,iz2,nnn,epsi)
C--     Max deviation
        eps = max(eps,epsi)
      enddo

C--   Validate the pdfs
      do i = 1,nnn
        call sqcValidate(ff,iff(i))
      enddo

      return
      end

C     =======================================================
      subroutine sqcNNupdn(
     +   ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,iz1,iz2,nnn,eps)
C     =======================================================

C--   Evolution from iz1 to iz2 with fixed nf in a single subgrid
C--   Does iterative downward evolution
C--
C--   ww         (in)  weight table store
C--   iw(i,j,k)  (in)  weight table identifiers of Pij for each order k
C--   ec         (in)  expansion coefficient store
C--   iec(i,j,k) (in)  coefficient table id corresponding to Pij
C--   idim       (in)  first 2 dims of iw and iec
C--   ff         (in)  pdf store
C--   iff(i)     (in)  pdf table identifiers for i = 1,nnn
C--   ig         (in)  subgrid index
C--   nopt       (in)  number of perturbative terms
C--   nf         (in)  number of flavours
C--   iz1        (in)  start of evolution
C--   iz2        (in)  end of evolution
C--   nnn        (in)  number of coupled equations
C--   eps        (out) max |lin-quad| mid-between grid points
C--
C--   Input :  Subgrid start values stored at iz = -ig in pdf table
C--   Output:  Subgrid pdf   values for iz1 to iz2     in pdf table
C--            Subgrid end   values stored at iz = -ig in pdf table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

      dimension ww(*),ff(*),ec(*)
      dimension iw(idim,idim,*),iec(idim,idim,*),iff(*)

*mb
      real time1, time2
      common/timing/spent(2)
*mb

*mb
      call cpu_time(time1)
*mb

C--   Number of subgrid points after cuts
      nyg = iqcIyMaxG(iymac2,ig)

C--   Now evolve up or down
      if(iz2.ge.iz1 .or.    !upward evolution  .or.
     +   ioy2.eq.2  .or.    !lin interpolation .or.
     +   niter6.lt.0) then  !always use current interpolation scheme

C--     Upward evolution: just evolve
        call sqcNNsubg(
     +  ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,nyg,iz1,iz2,nnn,eps)

      else

C--     Proceed with iterative downward evolution
        if0rem = -mxg0-1
        islast = -mxg0-2
        ioyrem =  ioy2
C--     Linear downward evolution
        ioy2 = 2
        call sqcNNsubg(
     +  ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,nyg,iz1,iz2,nnn,eps)
C--     Iterate?
        if(niter6.eq.0) then
C--       We want linear downward evolution w/o iteration
           ioy2 = ioyrem
          return
        endif
C--     Now proceed with iteration
C--     Copy startvalues from iz1 to the bins if0rem and islast
        do i = 1,nnn
          call sqcTcopyType5(ff,iff(i),iz1,if0rem)
          call sqcTcopyType5(ff,iff(i),iz1,islast)
        enddo
C--     Quad upward evolution
        ioy2 = 3
        call sqcNNsubg(
     +  ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,nyg,iz2,iz1,nnn,eps)
        niter = 0
C--     Now iterate
        do while(niter.lt.niter6)
C--       Set new startvalue
          do i = 1,nnn
            if0 = iqcG5ijk(ff,1,if0rem,iff(i))-1
            if1 = iqcG5ijk(ff,1,   -ig,iff(i))-1
            isl = iqcG5ijk(ff,1,islast,iff(i))-1
            do iy = 1,nyy2(0)
              ff(if1+iy) = ff(if0+iy) + ff(isl+iy)-ff(if1+iy)
            enddo
          enddo
C--       Linear downward evolution
          ioy2 = 2
          call sqcNNsubg(
     +    ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,nyg,iz1,iz2,nnn,eps)
C--       Remember start value
          do i = 1,nnn
            call sqcTcopyType5(ff,iff(i),iz1,islast)
          enddo
C--       Quad upward evolution
          ioy2 = 3
          call sqcNNsubg(
     +    ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,nyg,iz2,iz1,nnn,eps)
          niter = niter+1
        enddo
C--     End of iteration

C--     Restore interpolation degree
        ioy2 = ioyrem
C--     Loop over evolved pdfs
        do i = 1,nnn
C--       Store endpoint of evolution
          call sqcTcopyType5(ff,iff(i),iz2,-ig)
C--       Restore the original starting pdf
          call sqcTcopyType5(ff,iff(i),if0rem,iz1)
        enddo

      endif

*mb
      call cpu_time(time2)
      spent(2) = spent(2)+time2-time1
*mb

      return
      end

C     ==============================================================
      subroutine sqcNNsubg(
     +      ww,iw,ec,iec,idim,ff,iff,ig,nopt,nf,nyg,iz1,iz2,nnn,eps)
C     ==============================================================

C--   Evolution from iz1 to iz2 with fixed nf in a single subgrid
C--
C--   ww         (in)  weight table store
C--   iw(i,j,k)  (in)  weight table identifiers of Pij for each order k
C--   ec         (in)  expansion coefficient store
C--   iec(i,j,k) (in)  coefficient table id corresponding to Pij
C--   idim       (in)  first 2 dims of iw and iec
C--   ff         (in)  pdf store
C--   iff(i)     (in)  pdf table identifiers for i = 1,nnn
C--   ig         (in)  subgrid index
C--   nopt       (in)  number of perturbative terms
C--   nf         (in)  number of flavours
C--   nyg        (in)  upper index y-subgrid after cuts
C--   iz1        (in)  start of evolution
C--   iz2        (in)  end of evolution
C--   nnn        (in)  number of coupled equations
C--   eps        (out) max |lin-quad| mid-between grid points
C--
C--   Input :  Subgrid start values stored at iz = -ig in pdf table
C--   Output:  Subgrid pdf   values for iz1 to iz2     in pdf table
C--            Subgrid end   values stored at iz = -ig in pdf table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension ww(*),ff(*),ec(*)
      dimension iw(idim,idim,*),iec(idim,idim,*),iff(*)

      dimension vvv(mce0*mce0*mxx0),aaa(mce0*mxx0),bbb(mce0*mxx0)
      dimension fff(mxx0)
      dimension ivv(mce0,mce0),iaa(mce0),ibb(mce0)
      dimension sss(mce0*mce0*mxx0),hhh(mce0*mxx0)
      dimension iss(mce0,mce0),ihh(mce0)
      dimension umat(mxx0),usum(mxx0)

*mb
      real time1, time2
      common/timing/spent(2)
*mb

*mb
      call cpu_time(time1)
*mb

C--   Initialization
C--   --------------
      eps = 0.D0
      iy1 = 1
      do i = 1,nyg
        umat(i) = 0.D0
        usum(i) = 0.D0
      enddo
      do i = 1,nnn*nnn*nyg
        sss(i) = 0.D0
      enddo
      do i = 1,nnn*nyg
        hhh(i) = 0.D0
      enddo
C--   Adresses
      ih0 = 0
      is0 = 0
      do i = 1,nnn
        ibb(i) = ih0+1
        ih0    = ih0+nyg
        do j = 1,nnn
          iss(i,j) = is0+1
          is0      = is0+nyg
        enddo
      enddo

C--   Direction of evolution (isign) and first point after it1 (next)
      isign = 1
      next  = iz1+1
      if(iz2.lt.iz1) then 
        isign = -1
        next  = iz1-1
      endif

C--   Calculate vector b at input scale iz1 and store start value
C--   -----------------------------------------------------------
C--   Grid spacing delta = z(next)-z (divided by 2)
      delt = 0.5*abs(zgrid2(next)-zgrid2(iz1))
C--   Transformation matrix divided by delta (U)        
      do i = 1,nmaty2(ioy2)
        umat(i)   = smaty2(i,ioy2)/delt
      enddo
C--   Find t-index
      it1   = itfiz5(iz1)
C--   Weight matrix and  V = U+W at t1
      ityp = 0
      iw0  = 0
      if0  = 0
      do i = 1,nnn
        do j = 1,nnn
          ityp = ityp+1
          if(i.eq.j) then
            do iy = 1,nyg
              vvv(iw0+iy) = umat(iy)
            enddo
          else
            do iy = 1,nyg
              vvv(iw0+iy) = 0.D0
            enddo
          endif
          do k = 1,nopt
            if(iw(i,j,k).ne.0) then
              as = ec(iqcG6ij(ec,iz1,iec(i,j,k)))
              ia = iqcGaddr(ww,1,it1,nf,ig,iw(i,j,k))-1
              do iy = 1,nyg
                ia = ia+1
                vvv(iw0+iy) = vvv(iw0+iy) + isign*as*ww(ia)
              enddo
            endif
          enddo
          ivv(i,j) = iw0+1
          iw0      = iw0+nyg
        enddo
        iaa(i) = if0+1
        if0    = if0+nyg
C--     Get the starting F values from iz = -ig and convert to A values
        call sqcNNgetVj(ff,iff(i),-ig,ig,nyg,fff)
        call sqcNNFjtoAj(fff,aaa(iaa(i)),nyg)
C--     Put F values at iz1 in the pdf table
        call sqcNNputVj(ff,iff(i),iz1,ig,nyg,fff)
C--     Check for oscillation: max deviation quad from lin spline at iz1
        eps = max(eps,dqcNNGetEps(aaa,iaa(i),nyg))
      enddo

C--   Now make the matrix multiplication Va = b
      if(nnn.eq.1) then
        call sqcNSmult(vvv,nyg,aaa,bbb,nyg)
      elseif(nnn.eq.2) then
        iv1 = ivv(1,1)  !FF
        iv2 = ivv(1,2)  !FG
        iv3 = ivv(2,1)  !GF
        iv4 = ivv(2,2)  !GG
        call sqcSGmult(
     +       vvv(iv1),vvv(iv2),vvv(iv3),vvv(iv4),nyg,
     +       aaa(iaa(1)),aaa(iaa(2)),bbb(ibb(1)),bbb(ibb(2)),nyg)
      else
        call sqcNNmult(vvv,ivv,aaa,iaa,bbb,ibb,nnn,nyg,nyg,mce0)
      endif

C--   Evolution loop over t
C--   ---------------------
      do iz = next,iz2,isign

C--     Find t-index
        it  = itfiz5(iz)

C--     Weight matrix and V matrix at t
        iw0  = 0
        if0  = 0
        ityp = 0
        do i = 1,nnn
          do j = 1,nnn
            ityp = ityp+1
            if(i.eq.j) then
              do iy = 1,nyg
                vvv(iw0+iy) = umat(iy)
              enddo
            else
              do iy = 1,nyg
                vvv(iw0+iy) = 0.D0
              enddo
            endif
            do k = 1,nopt
              if(iw(i,j,k).ne.0) then
                as = ec(iqcG6ij(ec,iz,iec(i,j,k)))
                ia = iqcGaddr(ww,1,it,nf,ig,iw(i,j,k))-1
                do iy = 1,nyg
                  ia = ia+1
                  vvv(iw0+iy) = vvv(iw0+iy) - isign*as*ww(ia)
                enddo
              endif
            enddo
            ivv(i,j) = iw0+1
            iw0      = iw0+nyg
          enddo
        enddo

C--     Solve V a = b
        if(nnn.eq.1) then
          call sqcNSeqs(vvv,nyg,aaa,bbb,nyg)
        elseif(nnn.eq.2) then
          iv1 = ivv(1,1)  !FF
          iv2 = ivv(1,2)  !FG
          iv3 = ivv(2,1)  !GF
          iv4 = ivv(2,2)  !GG
          call sqcSGeqs(vvv(iv1),vvv(iv2),vvv(iv3),vvv(iv4),
     +         aaa(iaa(1)),aaa(iaa(2)),bbb(ibb(1)),bbb(ibb(2)),nyg)
        else
          call sqcNNeqs(vvv,ivv,aaa,iaa,bbb,ibb,nnn,nyg,mce0,ierr)
        endif

C--     Convert A to F and store in pdf table
        do i = 1,nnn
          call sqcNNAjtoFj(aaa(iaa(i)),fff,nyg)
          call sqcNNputVj(ff,iff(i),iz,ig,nyg,fff)
        enddo

C--     Update b for the next iteration: not at last iteration
        if(iz.ne.iz2) then
C--       Grid spacing delta = z(next)-z (divided by 2)
          delt = 0.5*abs(zgrid2(iz+isign)-zgrid2(iz))
C--       Calculate U(next) and S = U(next) + U(now); band only       
          do i = 1,nmaty2(ioy2)
            sbnext    = smaty2(i,ioy2)/delt
            usum(i)   = umat(i)+sbnext
            umat(i)   = sbnext
          enddo
C--       Make nnn*nnn copies of S
          jss = 0
          do i = 1,nnn
            ihh(i) = (i-1)*nyg+1
            do j = 1,nnn
              iss(i,j) = jss+1
              do k = 1,nmaty2(ioy2)
                jss = jss+1
                if(i.eq.j) then
                  sss(jss) = usum(k)
                else
                  sss(jss) = 0.D0
                endif
              enddo
            enddo
          enddo
C--       Now make matrix multiplication h = Sa
          if(nnn.eq.1) then
            call sqcNSmult(sss,nmaty2(ioy2),aaa,hhh,nyg)
          elseif(nnn.eq.2) then
            is1 = iss(1,1)  !FF
            is2 = iss(1,2)  !FG
            is3 = iss(2,1)  !GF
            is4 = iss(2,2)  !GG
            call sqcSGmult(
     +           sss(is1),sss(is2),sss(is3),sss(is4),nmaty2(ioy2),
     +           aaa(iaa(1)),aaa(iaa(2)),hhh(ihh(1)),hhh(ihh(2)),nyg)
          else
            call sqcNNmult(sss,iss,aaa,iaa,hhh,ihh,nnn,nyg,nmaty2(ioy2),
     +                     mce0)
          endif
C--       Update b(next) = h(now)-b(now)
          do i = 1,nnn
            ia = ihh(i)-1
            do iy = 1,nyg
              ia      = ia+1
              bbb(ia) = hhh(ia)-bbb(ia)
            enddo
          enddo

C--     To do at the end point of the evolution...
        else
          do i = 1,nnn
C--         Put end values at iz = -ig in the pdf table
            call sqcTcopyType5(ff,iff(i),iz2,-ig)
C--         Check for oscillation: max deviation quad from lin spline
            eps = max(eps,dqcNNGetEps(aaa,iaa(i),nyg))
          enddo
        endif

C--   End of loop over t
      enddo

C--   Thats it...

*mb
      call cpu_time(time2)
      spent(1) = spent(1)+time2-time1
*mb

      return
      end

C==   ==================================================================
C==   NxN matrix routines ==============================================
C==   ==================================================================

C--   ==================================================
      subroutine sqcNNmult(w,iw,a,ia,b,ib,n,m,nbnd,idim)
C--   ==================================================

C--   Multiply  W a = b
C--
C--   where W is an n*n array of lower band Toeplitz matrices
C--   of dimension m, and bandwith nbnd, a is a vector of n arrays of
C--   dimension m and b is a result vector of n arrays of dimension m.
C--   A Toeplitz is stored as a vector of dimension nbnd
C--
C--   |w1             |
C--   |w2 w1          |
C--   |w3 w2 w1       |  with, in this example, m = 5 and nbnd = 3 
C--   |   w3 w2 w1    |
C--   |      w3 w2 w1 |
C--
C--   (in)  w          Linear array with Toeplitz matrices W_alfbet           
C--   (in)  iw(n,n)    Adresses of first element of each W_alfbet
C--   (out) a          Linear array with input vectors a
C--   (in)  ia(n)      Adresses of first element of each a_alf
C--   (in)  b          Linear array with result vectors b
C--   (in)  ib(n)      Adresses of first element of each b_alf 
C--   (in)  n          Number of matrices  alf,bet = 1,...,n
C--   (in)  m          Dimension of each Toeplitz matrix ix = 1,...,m
C--   (in)  nbnd       Matrix bandwidth
C--   (in)  idim       First dim of iw as declared in the calling routine

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*),a(*),b(*)
      dimension iw(idim,*),ia(*),ib(*)

C--   Loop over all elements ix of the pdfs
      do ix = 1,m
        j1 = max(ix+1-nbnd,1)
C--     Loop over all pdfs ialf
        do ialf = 1,n
          balfix = 0.D0
C--       Loop over all pdfs ibet
          do ibet = 1,n
            iwalfbet = iw(ialf,ibet)
            iabet    = ia(ibet)
C--         Loop over the elements jx
            do jx = j1,ix
              walfbetixjx = w(iwalfbet+ix-jx)
              abetjx      = a(iabet+jx-1)
              balfix      = balfix + walfbetixjx*abetjx
            enddo
          enddo
          b(ib(ialf)+ix-1) = balfix
        enddo
      enddo

      return
      end

C     =================================================
      subroutine sqcNNeqs(w,iw,a,ia,b,ib,n,m,idim,ierr)
C     =================================================

C--   Solves coupled triangular equations
C--
C--            W a = b
C--
C--   where W is an n*n array of lower triangular Toeplitz matrices
C--   of dimension m, a is a solution vector of n arrays of dimension m
C--   and b is a right-hand side vector of n arrays of dimension m.
C--   A Toeplitz matrix T is stored as a vector of dimension m: 
C--
C--                                          | t1  0  0  0  0 |
C--                                          | t2 t1  0  0  0 |
C--   T(i,j<=i) = t(i-j+1)  e.g. for 5 x 5:  | t3 t2 t1  0  0 |
C--                                          | t4 t3 t2 t1  0 |
C--                                          | t5 t4 t3 t2 t1 |
C--
C--   (in)  w          Linear array with Toeplitz matrices W_alfbet           
C--   (in)  iw(n,n)    Adresses of first element of each W_alfbet
C--   (out) a          Linear array with solution vectors a
C--   (in)  ia(n)      Adresses of first element of each a_alf
C--   (in)  b          Linear array with right-hand vectors b
C--   (in)  ib(n)      Adresses of first element of each b_alf 
C--   (in)  n          Number of coupled equations  alf,bet = 1,...,n
C--   (in)  m          Dimension of each Toeplitz matrix ix = 1,...,m
C--   (in)  idim       First dimension of iw declared in the calling routine
C--   (out) ierr       0 = OK, #0 = Singularity encountered

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*), a(*), b(*)
      dimension iw(idim,*),ia(*),ib(*)
      dimension vv(mce0,mce0),bb(mce0)
      dimension iwork(mce0)

C--   Loop over all elements ix of the pdfs
      do ix = 1,m
C--     Loop over all pdfs ialf
        do ialf = 1,n
          sumalfix = 0.D0
C--       Loop over all pdfs ibet
          do ibet = 1,n
            iwalfbet = iw(ialf,ibet)
            iabet    = ia(ibet)
C--         Loop over the elements jx
            do jx = 1,ix-1
              valfbetixjx = w(iwalfbet+ix-jx)
              abetjx      = a(iabet+jx-1)
              sumalfix    = sumalfix + valfbetixjx*abetjx
C--         End of loop over elements jx
            enddo
            vv(ialf,ibet) = w(iwalfbet)
C--       End of loop over all pdfs ibet
          enddo
          balfix   = b(ib(ialf)+ix-1)
          bb(ialf) = balfix-sumalfix
C--     End of loop over all pdfs ialf
        enddo
C--     Now solve the nxn matrix equation for each ix
        call smb_dmeqn(n,vv,mce0,iwork,ierr,1,bb)
        if(ierr.ne.0) stop 'sqcNNeqs singular matrix encountered'
C--     Store the result
        do ialf = 1,n
          a(ia(ialf)+ix-1) = bb(ialf)
        enddo 
C--   End of loop over all elements ix of the pdfs
      enddo

      return
      end

C==   ==================================================================
C==   NxN pdf utilities ================================================
C==   ==================================================================

C     ==============================================
      double precision function dqcNNGetEps(aa,i,ny)
C     ==============================================

C--   Get max deviation quad-lin at midpoints
C--
C--   aa      (in)  : Array with A coefficients
C--   i       (in)  : Position of first coefficient in aa
C--   ny      (in)  : upper loop index in y

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dimension aa(*)
      dimension epsi(mxx0)
      
      dqcNNGetEps = 0.D0
      if(ioy2.ne.3) return                

C--   Now get vector of deviations
      call sqcDHalf(ioy2,aa(i),epsi,ny)
C--   Max deviation
      do iy = 1,ny
        dqcNNGetEps = max(dqcNNGetEps,abs(epsi(iy)))
      enddo

      return
      end

C     ==========================================
      subroutine sqcNNgetVj(ww,id,iz,ig,nyg,buf)
C     ==========================================

C--   Contract the subgrid points of a pdf table into the array buf 

C--   ww     (in)  Pdf store
C--   id     (in)  Global pdf identifier
C--   iz     (in)  Z-bin index
C--   ig     (in)  Subgrid index
C--   nyg    (in)  Number of subgrid points after cuts
C--   buf    (out) Array with values of all subgrid points i = 1,...,nyg

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dimension ww(*),buf(*)

      ia = iqcG5ijk(ww,1,iz,id)-1
      do iyg = 1,nyg
        iy0      = iy0fiyg2(iyg,ig)
        buf(iyg) = ww(ia+iy0)
      enddo

      return
      end

C     ==========================================
      subroutine sqcNNputVj(ww,id,iz,ig,nyg,buf)
C     ==========================================

C--   Disperse the array buf over the subgrid points of a pdf table

C--   ww     (in)  Pdf store
C--   id     (in)  Pdf identifier
C--   iz     (in)  Z-bin index
C--   ig     (in)  Subgrid index
C--   nyg    (in)  Number of subgrid points after cuts
C--   buf    (in)  Array with values of all subgrid points i = 1,...,nyg

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dimension ww(*),buf(*)

      ia = iqcG5ijk(ww,1,iz,id)-1
      do iyg = 1,nyg
        iy0        = iy0fiyg2(iyg,ig)
        ww(ia+iy0) = buf(iyg)
      enddo

      return
      end

C     ====================================
      subroutine sqcNNFjtoAj(ffin,aout,ny)
C     ====================================

C--   Convert an array of F values to an array of A values
C--   The F values must come from an equidistant (sub)grid

C--   ffin   (in)   Array of ny F values: ffin(i) i = 1,2,3,...,ny
C--   aout   (out)  Array of ny A values: aout(i) i = 1,2,3,...,ny
C--   ny     (in)   Number of points in ffin and aout
C--
C--   Remark: this conversion does not depend on the (sub)grid spacing

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      
      dimension ffin(*), aout(*)

C--   Convert pdf values to spline coefficients
      call sqcNSeqs(smaty2(1,ioy2),nmaty2(ioy2),aout,ffin,ny)

      return
      end

C     ====================================
      subroutine sqcNNAjtoFj(aain,fout,ny)
C     ====================================

C--   Convert an array of A values to an array of F values
C--   The A values must come from an equidistant (sub)grid

C--   aain   (in)   Array of ny A values: aain(i) i = 1,2,3,...,ny
C--   fout   (out)  Array of ny F values: fout(i) i = 1,2,3,...,ny
C--   ny     (in)   Number of points in aain and fout
C--
C--   Remark: this conversion does not depend on the (sub)grid spacing

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      
      dimension aain(*), fout(*)
      
C--   Convert spline coefficients to pdf values
      call sqcNSmult(smaty2(1,ioy2),nmaty2(ioy2),aain,fout,ny)

      return
      end
