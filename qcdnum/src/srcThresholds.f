C
C--  Threshold routines on the file srcThresholds.f
C--  ----------------------------------------------

C--   subroutine sqcChkIqh(nq,nfix,iqh,iq1,iq2,nfmin,nfmax,ierr)
C--   subroutine sqcChkRqh(qmin,qmax,rqhin,rqhout,ierr)
C--   subroutine sqcThrFFNS(nf)
C--   subroutine sqcThrMFNS(nf,rc,rb,rt)
C--   subroutine sqcThrVFNS(nfix,iqh,nfmin,nfmax)
C--
C--   subroutine sqcRmass2(fthr,rthr)
C--
C--   subroutine sqcNfTab(it0)
C--
C--   integer function izfit2(it)
C--   integer function izfitL2(it)
C--   integer function izfitU2(it)
C--   integer function itfiz2(iz)
C--   integer function nffiz2(iz)
C--
C--   integer function isfromitL(it,itlims,iz,nf)
C--   integer function isfromitU(it,itlims,iz,nf)

C=======================================================================
C===  Handle threshold input ===========================================
C=======================================================================

C     ==========================================================
      subroutine sqcChkIqh(nq,nfix,iqh,iq1,iq2,nfmin,nfmax,ierr)
C     ==========================================================

C--   Checks FFNS/VFNS input thresholds iqc, iqb, iqt
C--   On input, any of the iqc,b,t may be outside the grid range
C--
C--   nq        (in) : number of q gridpoints
C--   nfix      (in) : [0,1] = VFNS and [3,6] = FFNS else ierr = 1
C--   iqh(4:6)  (in) : input iqc,b,t
C--   iq1(3:6) (out) : lower limit of nf = 3,4,5,6
C--   iq2(3:6) (out) : upper limit of nf = 3,4,5,6
C--   nfmin    (out) : lower value of nf
C--   nfmax    (out) : upper value of nf
C--   ierr     (out) : 0 = OK
C--                    1 = invalid nfix
C--                    2 = no thresholds found inside q-grid
C--                    3 = found non-consecutive thresholds
C--                    4 = found thresholds not ascending or too close
C--
C--  Remark: Note boundary overlap iq2(nf) = iq1(nf+1)
C--
C--  Example:   iqc   iqb   iqt
C--              |     |     |
C--  iq->  1  2  3  4  5  6  7  8  9          nf =   3  4  5  6
C--        3  3  3                           ------------------
C--              4  4  4                     iqh =   x  3  5  7   (in)
C--                    5  5  5               iq1 =   1  3  5  7   (out)
C--                          6  6  6         iq2 =   3  5  7  9   (out)
C--

      implicit double precision (a-h,o-z)

      dimension iqh(4:6),iq1(3:6),iq2(3:6)

C--   Initialise
      ierr = 0
      do nf = 3,6
        iq1(nf) = 0
        iq2(nf) = 0
      enddo

C--   FFNS or VFNS thats the question
      if(nfix.eq.0 .or. nfix.eq.1) then
C--     VFNS
        continue
      elseif(nfix.ge.3 .and. nfix.le.6) then
C--     FFNS
        nfmin = nfix
        nfmax = nfix
        iq1(nfix) = 0
        iq2(nfix) = nq+1
        return
      else
C--     Invalid nfix
        ierr = 1
        return
      endif

C--   Now handle VFNS
      nflast = 0
      itlast = 0
      nfmin  = 0

      do nf = 4,6
        iqt = iqh(nf)
        if(iqt.ge.1 .and. iqt.le.nq) then
C--       Found threshold inside grid
          if(nfmin.eq.0) then
C--         Found first threshold
            iq1(nf-1) = 1
            nfmin     = nf-1
          else
C--         Threshold already found
            if(nflast.ne.nf-1) then
C--           Found non-consecutive
              ierr = 3
              return
            elseif(iqt.lt.itlast+2) then
C--           Found non-ascending or too close
              ierr = 4
              return
            endif
          endif
C--       Set threshold
          iq2(nf-1) = iqt
          iq1(nf)   = iqt
          iq2(nf)   = nq
          nfmax     = nf
          nflast    = nf
          itlast    = iqt
        endif
      enddo

C--   No thresholds found
      if(nfmin.eq.0) ierr = 2

      return
      end

C     =================================================
      subroutine sqcChkRqh(qmin,qmax,rqhin,rqhout,ierr)
C     =================================================

C--   Checks input thresholds on the renormalisation scale
C--   On input, any of the rqh may be outside the grid range
C--
C--   qmin         (in) : lower  limit of the q2 grid
C--   qmax         (in) : upper  limit of the q2 grid
C--   rqhin(4:6)   (in) : input  rqc,b,t
C--   rqhout(4:6) (out) : output rqc,b,t
C--   ierr        (out) : 0 = OK
C--                       1 = no thresholds found inside q-grid
C--                       2 = found non-consecutive thresholds
C--                       3 = thresholds not ascending or too close

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      dimension rqhin(4:6), rqhout(4:6)

C--   Initialise
      ierr   = 0
      nflast = 0
      rqlast = 0.D0
      nfmax  = 0
      nfmin  = 0

      do nf = 4,6
        rqh = rqhin(nf)
        if(rqh.ge.qmin .and. rqh.le.qmax) then
C--       Found threshold inside grid
          if(nfmin.eq.0) then
C--         Found first threshold
            nfmin      = nf-1
          else
C--         Threshold already found
            if(nflast.ne.nf-1) then
C--           Found non-consecutive
              ierr = 2
              return
            elseif(rqh.lt.1.02*rqlast) then
C--           Found non-ascending or too close
              ierr = 3
              return
            endif
          endif
C--       Set threshold
          rqhout(nf) = rqh
          rqlast     = rqh
          nfmax      = nf
          nflast     = nf
        endif
      enddo

C--   No thresholds found
      if(nfmin.eq.0) then
        ierr = 1
        return
      endif

C--   Thresholds below the grid
      do nf = 4,nfmin
        rqhout(nf) = 0.0001*nf*qlimd6
      enddo
C--   Thresholds above the grid
      do nf = nfmax+1,6
        rqhout(nf) = 1000.0*nf*qlimu6
      enddo

      return
      end

C     =========================
      subroutine sqcThrFFNS(nf)
C     =========================

C--   Define thresholds for the fixed flavor number scheme
C--
C--   Input   nf           number of flavors [3-6]
C--   Output  qthrs6(4:6)  thresholds on the mu2 grid
C--           tthrs6(4:6)  thresholds on t = log(mu2) grid
C--           rthrs6(4:6)  thresholds on the renor scale

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      if    (nf.eq.3) then
        qthrs6(4) = 4000.0*qlimu6
        qthrs6(5) = 5000.0*qlimu6
        qthrs6(6) = 6000.0*qlimu6
      elseif(nf.eq.4) then
        qthrs6(4) = 0.0004*qlimd6
        qthrs6(5) = 5000.0*qlimu6
        qthrs6(6) = 6000.0*qlimu6
      elseif(nf.eq.5) then
        qthrs6(4) = 0.0004*qlimd6
        qthrs6(5) = 0.0005*qlimd6
        qthrs6(6) = 6000.0*qlimu6
      elseif(nf.eq.6) then
        qthrs6(4) = 0.0004*qlimd6
        qthrs6(5) = 0.0005*qlimd6
        qthrs6(6) = 0.0006*qlimd6
      else
        stop 'sqcThrFFNS: invalid nf'
      endif

      do i = 4,6
        rthrs6(i) = qthrs6(i)
        tthrs6(i) = log(qthrs6(i))
      enddo

      nfix6  = nf
      nfmin6 = nf
      nfmax6 = nf

      return
      end
      
C     ==================================      
      subroutine sqcThrMFNS(nf,rc,rb,rt)
C     ================================== 

C--   Mixed flavour number scheme thresholds

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      
C--   Set nfix6 and thresholds on the fact scale      
      call sqcThrFFNS(nf)
C--   Set thresholds on the renormalisation scale
      rthrs6(4) = rc
      rthrs6(5) = rb
      rthrs6(6) = rt
      
C--   Flag MVNS
      nfix6 = -nf

      return
      end

C     ===========================================
      subroutine sqcThrVFNS(nfix,iqh,nfmin,nfmax)
C     ===========================================

C--   Set thresholds in the vfns
C--
C--   nfix         (in) : 0=dynamic 1=intrinsic heavy flavours
C--   iqh(4:6)     (in) : thresholds iqc,iqb,iqt
C--   nfmin        (in) : lowest  nf in q2-grid
C--   nfmax        (in) : highest nf in q2-grid
C--
C--   Output in common block /qpars6/
C--
C--   nfix6             : set to 0,1 for vfns
C--   qthrs6(4:6)       : thresholds on the mu2 grid
C--   tthrs6(4:6)       : thresholds on t = log(mu2) grid
C--   rthrs6(4:6)       : thresholds on the renor scale
C--
C--   Can ony be called when t-grid is defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension iqh(4:6)

      if(nfix.ne.0 .and. nfix.ne.1) stop 'sqcThrVFNS: nfix not 0 or 1'

C--   Thresholds below the grid
      do nf = 4,nfmin
        qthrs6(nf) = 0.0001*nf*qlimd6
      enddo
C--   Thresholds inside the grid
      do nf = nfmin+1,nfmax
        qthrs6(nf) = exp(tgrid2(iqh(nf)))
      enddo
C--   Thresholds above the grid
      do nf = nfmax+1,6
        qthrs6(nf) = 1000.0*nf*qlimu6
      enddo

      do i = 4,6
        tthrs6(i) = log(qthrs6(i))
      enddo

      call sqcRmass2(qthrs6,rthrs6)

      nfix6  = nfix
      nfmin6 = nfmin
      nfmax6 = nfmax

      return
      end

C     ===============================
      subroutine sqcRmass2(fthr,rthr)
C     ===============================

C--   Convert threshold on the fact scale to the renor scale,
C--   using the scale factors stored in qpars6
C--
C--   fthr(4:6)  (in) : thresholds fcbt2 defined on the fact  scale
C--   rthr(4:6) (out) : thresholds rcbt2 defined on the renor scale

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      dimension fthr(4:6),rthr(4:6)

      do i = 4,6
        rthr(i) = aar6*fthr(i) + bbr6
      enddo

      return
      end

C=======================================================================
C===  Handle flavour subgrids in mu2 ===================================
C=======================================================================

C     ===============================
      subroutine sqcNfTab(w,kset,it0)
C     ===============================

C--   w      : workspace to put the subgrid tables
C--   kset   : table set in w
C--   it0 = 0: no debug printout
C--   it0 # 0: generate debug printout. it0 is the starting point of a
C--            fake evolution in this printout
C--
C--   The lists created by this routine are explained by the following
C--   example of 3 flavors on an 8-point t-grid: nf = 3 from t1-t3,
C--   nf = 4 from t3-t6 and nf = 5 from t6-t8 (note the one-point overlap)
C--   So we have three flavors (3,4,5) which are put in a list (nfsubg6)
C--   The ranges of these flavors are put in lists it1sub6, it2sub6
C--
C--   it --->    1  2  3  4  5  6  7  8     ntsubg6 = 3 --> 1  2  3  4
C--    t        t1 t2 t3 t4 t5 t6 t7 t8     nfsubg6         3  4  5  0
C--   nf3->      3  3  3  0  0  0  0  0     it1sub6         1  3  6  0
C--   nf4->      0  0  4  4  4  4  0  0     it2sub6         3  6  8  0
C--   nf5->      0  0  0  0  0  5  5  5
C--   nf6->      0  0  0  0  0  0  0  0
C--
C--   Zgrid is a working grid with each threshold doubled
C--
C--   iz --->    1  2  3  4  5  6  7  8  9 10       i --->  1  2  3  4
C--   nf         3  3  3  4  4  4  4  5  5  5      nfsubg6  3  4  5  0
C--   zgrid     t1 t2 t3 t3 t4 t5 t6 t6 t7 t8      iz1sub6  1  4  8  0
C--                                                iz2sub6  3  7 10  0
C--
C--                                                nf --->  3  4  5  6
C--                                                izminf6  1  4  8  0
C--                                                izmanf6  3  7 10  0
C--
C--   The following pointer arrays are used to map the indices:
C--
C--   i --->     1  2  3  4  5  6  7  8  9 10
C--   itfiz2     1  2  3  3  4  5  6  6  7  8             (it from  iz)
C--   izfit2     1  2  4  5  6  8  9 10                   (iz from  it)
C--   izfitU2    1  2  4  5  6  8  9 10                   (iz upper nf)
C--   izfitL2    1  2  3  5  6  7  9 10                   (iz lower nf)
C--   nffiz2     3  3  3  4  4  4  4  5  5  5             (nf from  iz)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      dimension ithr(4:6), iq1(3:6), iq2(3:6)
      dimension iz1f(3:6), iz2f(3:6), iz1b(3:6), iz2b(3:6)

      dimension w(*)

C--   Initialize
      do i = 1,4
        nfsubg6(i) = 0
        it1sub6(i) = 0
        it2sub6(i) = 0
        iz1sub6(i) = 0
        iz2sub6(i) = 0
      enddo
      do i = 3,6
        izminf6(i) = 0
        izmanf6(i) = 0
        iq1(i)     = 0
        iq2(i)     = 0
        iz1f(i)    = 0
        iz2f(i)    = 0
        iz1b(i)    = 0
        iz2b(i)    = 0
      enddo
      do i = 1,10
        itlim6(i) = 0
      enddo
      do i = 1,13
        izlim6(i) = 0
      enddo

C--   Get threshold indices
      do nf = 4,6
        ithr(nf) = iqcItfrmT(tthrs6(nf))
      enddo

C--   Flavour ranges
      iq1(nfmin6) = 1
      do nf = nfmin6+1,nfmax6
        iq1(nf) = ithr(nf)
      enddo
      do nf = nfmin6,nfmax6-1
        iq2(nf) = ithr(nf+1)
      enddo
      iq2(nfmax6) = ntt2

C--   Make zgrid lookup tables
      iz = 0
      do nf = nfmin6,nfmax6
        do iq = iq1(nf),iq2(nf)
          iz          = iz+1
        enddo
      enddo
      nzz2 = iz

C--   Subgrid limits it-grid
      ntsubg6 = 0
      do nf = nfmin6,nfmax6
        ntsubg6          = ntsubg6 + 1
        nfsubg6(ntsubg6) = nf                    !Used by EvolveOld only
        it1sub6(ntsubg6) = iq1(nf)               !Used by EvolveOld only
        it2sub6(ntsubg6) = iq2(nf)               !Used by EvolveOld only
      enddo

C--   Parameters for pointer functions in array itlim6(10)
C--   1          : number of subgrid limits, incl endpoints
C--   2--nlims+1 : subgrid limits (max 5)
C--   7--nlims+6 : nf for each subgrid (max 4)
      itlim6(1) = ntsubg6+1
      itlim6(2) = 1
      do i = 1,ntsubg6
        itlim6(i+1) = it1sub6(i)
        itlim6(i+6) = nfsubg6(i)
      enddo
      itlim6(ntsubg6+2) = ntt2

C--   Get izmin, izmax limits
      nlims      = itlim6(1)
      iz1sub6(1) = itlim6(2)
      do i = 2,nlims-1
        isub = isfromitU(itlim6(i+1),itlim6,iz1sub6(i)  ,nf)
        isub = isfromitL(itlim6(i+1),itlim6,iz2sub6(i-1),nf)
      enddo
      iz2sub6(nlims-1) = itlim6(nlims+1) + nlims - 2

C--   Parameters for pointer functions in array izlim6(13)
C--    1          : number of subgrids
C--    2--nlims+1 : iz1 lower subgrid limits (max 4)
C--    6--nlims+5 : iz2 upper subgrid limits (max 4)
C--   10--nlims+9 : nf for each subgrid (max 4)
      izlim6(1) = ntsubg6
      do i = 1,ntsubg6
        izlim6(i+1) = iz1sub6(i)
        izlim6(i+5) = iz2sub6(i)
        izlim6(i+9) = nfsubg6(i)
      enddo

C--   Subgrid pointer tables
C--   Table 1: iz > 0 --> it(iz)
C--            iz < 0 --> nf(iz)
      ig = 1000*kset + 701
      ia = iqcG7ij(w,0,ig)
      do iz = 1,nzz2
        isdum    = isfromiz(iz,izlim6,it,nf)
        w(ia+iz) = dble(it)
        w(ia-iz) = dble(nf)
      enddo

C--   Table 2: it > 0 --> iz(it) 456 grid
C--            it < 0 --> iz(it) 345 grid
      ig = 1000*kset + 702
      ia = iqcG7ij(w,0,ig)
      do it = 1,ntt2
        isdum    = isfromitU(it,itlim6,izU,nf)
        isdum    = isfromitL(it,itlim6,izL,nf)
        w(ia+it) = dble(izU)
        w(ia-it) = dble(izL)
      enddo

C--   These should be kept because still used
      do iz = 1,nzz2
        jz = nzz2+1-iz
        izmanf6(nffiz2(iz)) = iz
        izminf6(nffiz2(jz)) = jz
      enddo

C--   Fill zgrid
      do iz = 1,nzz2
        zgrid2(iz) = tgrid2(itfiz2(iz))
      enddo

C--   Store for later use
      itchm2 = ithr(4)
      itbot2 = ithr(5)
      ittop2 = ithr(6)

C--   Done...
      if(it0.eq.0) return

C--   Debug printout
      write(6,'(/'' sqcNFTAB evolution pointers:'')')
      write(6,'(/'' ic,b,t  '',3I4)') ithr(4),ithr(5),ithr(6)
      write(6,'( '' nf,mi,ma'',3I4)') ntsubg6,nfmin6,nfmax6
C--   Print ranges
      write(6,'(/'' nfsubg6 '',4I4)') nfsubg6
      write(6,'( '' it1sub6 '',4I4)') it1sub6
      write(6,'( '' it2sub6 '',4I4)') it2sub6
      write(6,'( '' iz1sub6 '',4I4)') iz1sub6
      write(6,'( '' iz2sub6 '',4I4)') iz2sub6
      write(6,'(/'' nf--->     3   4   5   6'')')
      write(6,'( '' izminf6 '',4I4)') izminf6
      write(6,'( '' izmanf6 '',4I4)') izmanf6
C--   Print limits
      write(6,'(/'' itlim6 nlims  '',I4 )') itlim6(1)
      write(6,'( '' itlim6 limits '',5I4)') (itlim6(i+1),i=1,itlim6(1))
      write(6,'( '' itlim6 subgnf '',4I4)')
     +                                    (itlim6(i+6),i=1,itlim6(1)-1)
C--   Pointer arrays
      if(nzz2.le.23) then
        write(6,'(/'' i ---> '',23I3)') (i,i=1,23)
        write(6,'( '' izfit  '',23I3)') (izfit2(i) ,i=1,23)
        write(6,'( '' izfitU '',23I3)') (izfitU2(i),i=1,23)
        write(6,'( '' izfitL '',23I3)') (izfitL2(i),i=1,23)
        write(6,'( '' itfiz  '',23I3)') (itfiz2(i) ,i=1,23)
        write(6,'( '' nffiz  '',23I3)') (nffiz2(i) ,i=1,23)
      endif
C--   Fake evolution loop for testing
      write(6,*) ' '
      write(6,*) ' it0 =',it0
      if(it0.gt.0) then
        iz0 = izfitU2(it0)
      else
        iz0 = izfitL2(-it0)
      endif
      ipr = 1
      call sparParTo5(1)
      call sqcEvPlan(iz0,nf0,nfmi,nfma,iz1f,iz2f,iz1b,iz2b,ipr,ie)

      return
      end

C=======================================================================
C===  Pointer functions ================================================
C=======================================================================

C     ===========================
      integer function izfit2(it)
C     ===========================

C--   Replaces the array /qgrid2.inc/ izfit2(it)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      isub = isfromitU(it,itlim6,izfit2,nf)

      return
      end

C     ============================
      integer function izfitL2(it)
C     ============================

C--   Replaces the array /qgrid2.inc/ izfitL2(it)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      isub = isfromitL(it,itlim6,izfitL2,nf)

      return
      end

C     ============================
      integer function izfitU2(it)
C     ============================

C--   Replaces the array /qgrid2.inc/ izfitU2(it)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      isub = isfromitU(it,itlim6,izfitU2,nf)

      return
      end

C     ===========================
      integer function itfiz2(iz)
C     ===========================

C--   Replaces the array /qgrid2.inc/ itfiz2(iz)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      isub = isfromiz(iz,izlim6,itfiz2,nf)

      return
      end

C     ===========================
      integer function nffiz2(iz)
C     ===========================

C--   Replaces the array /qgrid2.inc/ nffiz2(iz)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      isub = isfromiz(iz,izlim6,it,nffiz2)

      return
      end

C=======================================================================
C===  Subgrid index functions ==========================================
C=======================================================================

C     ===========================================
      integer function isfromitL(it,itlims,iz,nf)
C     ===========================================

C--   Get subgrid index 345 grid
C--
C--   it      (in)  : t-grid index
C--   itlims  (in)  : itlims(1)=nlims, 2--6=limits, 7--10=nf
C--   iz     (out)  : z-grid index
C--   nf     (out)  : number of flavours


      implicit double precision (a-h,o-z)

      dimension itlims(*)

      nlims = itlims(1)
      if(it.lt.itlims(2) .or. it.gt.itlims(nlims+1)) then
        isfromitl = 0
        stop 'ISFROMITL: it-index outside limits'
      elseif(nlims.gt.2) then
        isub = 0
        it1  = itlims(2)
        i    = 1
        do while(isub.eq.0)
          it2 = itlims(i+2)
          if(it.ge.it1 .and. it.le.it2) isub = i
          it1 = it2
          i   = i+1
        enddo
        isfromitL = isub
      elseif(nlims.eq.2) then
        isfromitL = 1
      else
        isfromitL = 0
        stop 'ISFROMITL: nlims < 2'
      endif

      iz = it + isfromitL - 1
      nf = itlims(6+isfromitL)

      end

C     ===========================================
      integer function isfromitU(it,itlims,iz,nf)
C     ===========================================

C--   Get subgrid index 456 grid
C--
C--   it      (in)  : t-grid index
C--   itlims  (in)  : itlims(1)=nlims, 2--6=limits, 7--10=nf
C--   iz     (out)  : z-grid index
C--   nf     (out)  : number of flavours

      implicit double precision (a-h,o-z)

      dimension itlims(*)

      nlims = itlims(1)
      if(it.lt.itlims(2) .or. it.gt.itlims(nlims+1)) then
        isfromitU = 0
        stop 'ISFROMITU: it-index outside limits'
      elseif(nlims.gt.2) then
        isub = 0
        it2  = itlims(nlims+1)
        i    = nlims
        do while(isub.eq.0)
          i   = i-1
          it1 = itlims(i+1)
          if(it.ge.it1 .and. it.le.it2) isub = i
          it2 = it1
        enddo
        isfromitU = isub
      elseif(nlims.eq.2) then
        isfromitU = 1
      else
        isfromitU = 0
        stop 'ISFROMITU: nlims < 2'
      endif

      iz = it + isfromitU - 1
      nf = itlims(6+isfromitU)

      end

C     ==========================================
      integer function isfromiz(iz,izlims,it,nf)
C     ==========================================

C--   Get subgrid index from iz grid
C--
C--   iz      (in)  : z-grid index
C--   izlims  (in)  : 1=nsubgrid, 2--5=iz1, 6--9=iz2, 10-13=nf
C--   it     (out)  : t-grid index
C--   nf     (out)  : number of flavours

      implicit double precision (a-h,o-z)

      dimension izlims(*)

      nsubg = izlims(1)
      if(iz.lt.izlims(2) .or. iz.gt.izlims(nsubg+5)) then
        isfromiz = 0
        stop 'ISFROMIZ: iz-index outside limits'
      elseif(nsubg.gt.1) then
        isub = 0
        i    = 1
        iz1  = izlims(i+1)
        do while(isub.eq.0)
          iz2 = izlims(i+5)
          if(iz.ge.iz1 .and. iz.le.iz2) isub = i
          iz2 = iz1
          i   = i+1
        enddo
        isfromiz = isub
      elseif(nsubg.eq.1) then
        isfromiz = 1
      else
        isfromiz = 0
        stop 'ISFROMIZ: nsubg < 1'
      endif

      it = iz - isfromiz + 1
      nf = izlims(9+isfromiz)

      end

