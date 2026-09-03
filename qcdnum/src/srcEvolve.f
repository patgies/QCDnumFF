
C--   This is the file srcEvolve.f with the standard evolution routines

C--   subroutine sqcEvolFG(ityp,mset,func,def,it0,epsi,nfheavy,ierr)
C--   subroutine sqcEvSGNS(ityp,mset,func,isns,n,it0,epsi,nfheavy,ierr)
C--   subroutine sqcEvFixNf(ityp,mset,nf,mode,start,iz1,iz2,eps)
C--   subroutine sqcSetStart(fun,tmatqf,nf,nfheavy,startu,startd)
C--   subroutine sqcEvPlan(iz0,nf0,nfmi,nfma,iz1f,iz2f,iz1b,iz2b,ip,ie)

C     ==============================================================
      subroutine sqcEvolFG(ityp,mset,func,def,it0,epsi,nfheavy,ierr)
C     ==============================================================

C--   Steering routine

C--   ityp     (in) : 1=unpol 2=pol 3=timelike
C--   mset     (in) : output pdf set
C--   func     (in) : starting parameterisations at it0
C--   def      (in) : flavour composition of input functions
C--   it0      (in) : starting scale
C--   epsi    (out) : spline oscillation indicator
C--   nfheavy (out) : largest heavy flavour
C--   ierr    (out) : 0  OK
C--                   1  Start point not inside grid
C--                   2  At least one evolution limit below alphas cut
C--                   3  Input pdfs not linearly independent
C--                   4  Invalid intrinsic heavy quark input

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      logical lmb_ne, nullrows

      external  func
      dimension def(13,12),tmatfq(13,13),tmatqf(13,13)

      dimension iz1f(3:6), iz2f(3:6), iz1b(3:6), iz2b(3:6)
      dimension startu(2,mxx0,12),startd(2,mxx0,12)
      dimension pdiff(2,mxx0,12),eps(12)

      dimension maskh(6,7:12)
C--                d  u  s  c  b  t
      data maskh / 1, 1, 1, 0, 1, 1,      !j =  7
     +             1, 1, 1, 0, 1, 1,      !j =  8
     +             1, 1, 1, 1, 0, 1,      !j =  9
     +             1, 1, 1, 1, 0, 1,      !j = 10
     +             1, 1, 1, 1, 1, 0,      !j = 11
     +             1, 1, 1, 1, 1, 0/      !j = 12

      if(idbug6.ne.0) write(lunerr1,'(/)')

C-1   Initialise -------------------------------------------------------
      ierr = 0
      if(idbug6.eq.2)      then
        call sqcNfTab(pars8,1,it0)
      endif
      do j = 1,13
        do i = 1,13
          tmatfq(i,j) = 0.D0
        enddo
      enddo
C--   Set all pdfs to zero
      do i = 0,12
        idg = iqcIdPdfLtoG(mset,i)
        call sqcPreset(idg,0.D0)
      enddo

C-2   Setup global pdf table ids ---------------------------------------
      idPdfi7(1,1) = iqcIdPdfLtoG(mset,1)     !Singlet
      idPdfi7(2,1) = iqcIdPdfLtoG(mset,0)     !Gluon
      do id = 2,12
        idPdfi7(1,id) = iqcIdPdfLtoG(mset,id) !Nonsinglets
        idPdfi7(2,id) = 0
      enddo

*mb   Print weight table indices
*      write(6,'(/'' IdSpfun ityp = '',I2)') ityp
*      write(6,'('' PQQ '',3I10)') (idspfun('PQQ',k,ityp),k=1,3)
*      write(6,'('' PQG '',3I10)') (idspfun('PQG',k,ityp),k=1,3)
*      write(6,'('' PGQ '',3I10)') (idspfun('PGQ',k,ityp),k=1,3)
*      write(6,'('' PGG '',3I10)') (idspfun('PGG',k,ityp),k=1,3)
*      write(6,'('' PPL '',3I10)') (idspfun('PPL',k,ityp),k=1,3)
*      write(6,'('' PMI '',3I10)') (idspfun('PMI',k,ityp),k=1,3)
*      write(6,'('' PVA '',3I10)') (idspfun('PVA',k,ityp),k=1,3)

*mb   Print all indices
*      write(6,'(/'' Pij indices ityp = '',I2)') ityp
*      write(6,'( '' Singlet gluon'')')
*      do k = 1,3
*        write(6,'('' iord = '',I3, '' Pij = '',4I10)') k,
*     +  ((idWijk7(i,j,k,1,ityp),j=1,2),i=1,2)
*      enddo
*      write(6,'( '' Nonsinglet NS+'')')
*      do k = 1,3
*        write(6,'('' iord = '',I3, '' Pij = '',4I10)') k,
*     +  ((idWijk7(i,j,k,2,ityp),j=1,2),i=1,2)
*      enddo
*      write(6,'( '' Nonsinglet NS-'')')
*      do k = 1,3
*        write(6,'('' iord = '',I3, '' Pij = '',4I10)') k,
*     +  ((idWijk7(i,j,k,3,ityp),j=1,2),i=1,2)
*      enddo
*      write(6,'( '' Nonsinglet VA'')')
*      do k = 1,3
*        write(6,'('' iord = '',I3, '' Pij = '',4I10)') k,
*     +  ((idWijk7(i,j,k,4,ityp),j=1,2),i=1,2)
*      enddo
*
*      write(6,'(/'' Alfas indices'')')
*      do k = 1,3
*        write(6,'('' iord = '',I3, '' Alf = '',4I10)') k,
*     +  ((idEijk7(i,j,k),j=1,2),i=1,2)
*      enddo
*
*      write(6,'(/'' Pdf indices mset = '',I2)') mset
*      do j = 1,12
*        write(6,'('' pdf  = '',I3, '' idg = '',2I10)') j,
*     +  (idPdfi7(i,j),i=1,2)
*      enddo

C-3   Setup the evolution limits and iz0 -------------------------------
C--   Set cuts in iz make sure to include it0
      it     = min(itmic2,abs(it0))
      izmic2 = izfit5(-it)
      it     = max(itmac2,abs(it0))
      izmac2 = izfit5( it)
C--   Start point in z
      iz0    = izfit5(it0)
C--   Determine how to proceed with the evolution
      call sqcEvPlan(iz0,nf0,nfmi,nfma,iz1f,iz2f,iz1b,iz2b,0,ierr)
C--   At least one limit outside grid or alphas cut
      if(ierr.ne.0) return

C-4   Handle the def array also swap indices ---------------------------
      nf = itfiz5(-iz0)
C--   Gluon entry
      tmatfq(1,7) = 1.D0
C--   Copy the 2nf x 2nf part of def to tmatfq
      do j = 1,2*nf
        do i = 1,nf
          tmatfq(j+1,7+i) = def(7+i,j)
          tmatfq(j+1,7-i) = def(7-i,j)
        enddo
      enddo
      nfheavy = nf
C--   Now check and copy intrinsic heavy flavour coefficients
      if(abs(nfix6).eq.1) then
C--     Check for zero rows
        nullrows = .false.
        do j = nf+1,6
          if(.not.nullrows) then
            nullrows = .true.
            do i = 1,6
              if(def(7-i,2*j-1).ne.0.D0) nullrows = .false.
              if(def(7+i,2*j-1).ne.0.D0) nullrows = .false.
              if(def(7-i,2*j  ).ne.0.D0) nullrows = .false.
              if(def(7+i,2*j  ).ne.0.D0) nullrows = .false.
            enddo
            if(.not.nullrows) nfheavy = j
          endif
        enddo
C--     Check coefficients are h+- only
        do j = 2*nf+1,2*nfheavy
          do i = 1,6
            if(lmb_ne(maskh(i,j)*def(7+i,j),0.D0,aepsi6)) then
              ierr = 4
              return
            endif
            if(lmb_ne(maskh(i,j)*def(7-i,j),0.D0,aepsi6)) then
              ierr = 4
              return
            endif
          enddo
        enddo
C--     Copy h+- coefficients
        do j = nf+1,nfheavy
          tmatfq(2*j  ,7-j) = def(7-j,2*j-1)
          tmatfq(2*j  ,7+j) = def(7+j,2*j-1)
          tmatfq(2*j+1,7-j) = def(7-j,2*j  )
          tmatfq(2*j+1,7+j) = def(7+j,2*j  )
        enddo
      endif
C--   For the rest can set anything linearly independent
      do j = nfheavy+1,6
        tmatfq(2*j  ,7-j) = 1.D0
        tmatfq(2*j+1,7+j) = 1.D0
      enddo

C-5   Set start values -------------------------------------------------
      call sqcGetMatQF(tmatfq, tmatqf, jerr)
      if(jerr.ne.0) then
        ierr = 3
        return
      endif
      call sqcSetStart( func, tmatqf, nf, nfheavy, startu, startd )

C-6   Evolution --------------------------------------------------------
C--   Forward
      if(idbug6.ge.1)
     +  write(lunerr1,'(/1X,17(''-''),'' forward '',17(''-''))')
      call sqcEvFixNf(ityp,mset,nf0,+1,startu,iz1f(nf0),iz2f(nf0),eps)
      do nf = nf0+1,nfma
C--     Calculate jumps
        call sqcDoJumps(ityp,mset,nf,iz1f(nf),+1,pdiff)
        call sqcEvFixNf(ityp,mset,nf,-1,pdiff,iz1f(nf),iz2f(nf),eps)
      enddo

C--   Backward
      if(idbug6.ge.1)
     +   write(lunerr1,'(1X,17(''-''),'' reverse '',17(''-''))')
      call sqcEvFixNf(ityp,mset,nf0,+1,startd,iz1b(nf0),iz2b(nf0),eps)
      do nf = nf0-1,nfmi,-1
C--     Calculate jumps
        call sqcDoJumps(ityp,mset,nf,iz1b(nf),-1,pdiff)
        call sqcEvFixNf(ityp,mset,nf,-1,pdiff,iz1b(nf),iz2b(nf),eps)
      enddo

C-7   Finalise ---------------------------------------------------------
C--   Re-set cuts in iz (which may have been stretched to include it0)
      izmic2 = izfit5(-itmic2)
      izmac2 = izfit5( itmac2)
C--   Validate the pdfs
      do id = 0,12
        idglobal = iqcIdPdfLtoG(mset,id)
        call sqcValidate(stor7,idglobal)
      enddo
C--   Max deviation
      epsi = 0.D0
      do id = 1,12
        epsi = max(epsi,eps(id))
      enddo

      if(idbug6.ne.0) write(lunerr1,'(/)')

      return
      end

C     =================================================================
      subroutine sqcEvSGNS(ityp,mset,func,isns,n,it0,epsi,nfheavy,ierr)
C     =================================================================

C--   Steering routine

C--   ityp   (in) : 1=unpol 2=pol 3=timelike
C--   mset   (in) : output pdf set
C--   func   (in) : starting parameterisations at it0
C--   isns   (in) : si/ns specifiers 1=SG 2=NS+ -1=V -2=NS-
C--   n      (in) : number of pdfs to evolve
C--   it0    (in) : starting scale
C--   epsi  (out) : spline oscillation indicator
C--   ierr  (out) : 0  OK
C--                 1  Start point not inside grid
C--                 2  At least one evolution limit below alphas cut

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      external  func

      dimension iz1f(3:6), iz2f(3:6), iz1b(3:6), iz2b(3:6)
      dimension startu(2,mxx0,12),startd(2,mxx0,12)
      dimension isns(*), eps(12)

C--   Pdf type SG, NS+, NS-, VA; number of coupled evolutions
      dimension iptype(-2:2),nevols(-2:2)
C--                -2  -1   0   1   2
      data iptype / 3,  4,  0,  1,  2 /
      data nevols / 1,  1,  0,  2,  1 /


      if(idbug6.ne.0) write(lunerr1,'(/)')

C-1   Initialise -------------------------------------------------------
      ierr = 0
      if(idbug6.eq.2)      then
        call sqcNfTab(pars8,1,it0)
      endif
C--   Set all pdfs to zero
      do i = 0,n
        idg = iqcIdPdfLtoG(mset,i)
        call sqcPreset(idg,0.D0)
      enddo

C-2   Setup global pdf table ids; singlet always at 1, if present ------
      do id = 1,n
        if(isns(id).eq.1) then
          idPdfi7(1,id) = iqcIdPdfLtoG(mset,1)                  !Singlet
          idPdfi7(2,id) = iqcIdPdfLtoG(mset,0)                    !Gluon
        else
          idPdfi7(1,id) = iqcIdPdfLtoG(mset,id)             !Nonsinglets
          idPdfi7(2,id) = 0
        endif
      enddo

C-3   Setup the evolution limits and iz0 -------------------------------
C--   Set cuts in iz make sure to include it0
      it     = min(itmic2,abs(it0))
      izmic2 = izfit5(-it)
      it     = max(itmac2,abs(it0))
      izmac2 = izfit5( it)
C--   Start point in z
      iz0    = izfit5(it0)
C--   Determine how to proceed with the evolution
      call sqcEvPlan(iz0,nf0,nfmi,nfma,iz1f,iz2f,iz1b,iz2b,0,ierr)
C--   At least one limit outside grid or alphas cut
      if(ierr.ne.0) return

*      nf    = itfiz5(-iz0)

C-4   Set start values -------------------------------------------------
      do iy = 1,iymac2
        y = ygrid2(iy)
        x = exp(-y)
        if(isns(1).eq.1) then
C--       Copy singlet/gluon into start
          do id = 0,1
            startu(2-id,iy,1) = func(id,x)
            startd(2-id,iy,1) = func(id,x)
          enddo
C--       Copy nonsinglets into start
          do id = 2,n
            startu(1,iy,id) = func(id,x)
            startd(1,iy,id) = func(id,x)
            startu(2,iy,id) = 0.D0
            startd(2,iy,id) = 0.D0
          enddo
        else
C--       Copy nonsinglets into start
          do id = 1,n
            startu(1,iy,id) = func(id,x)
            startd(1,iy,id) = func(id,x)
            startu(2,iy,id) = 0.D0
            startd(2,iy,id) = 0.D0
          enddo
        endif
      enddo

C-5   Evolution forward    ---------------------------------------------
      iz1   = iz1f(nf0)
      iz2   = iz2f(nf0)
      if(idbug6.ge.1) then
        iq1 = itfiz5(iz1)
        iq2 = itfiz5(iz2)
        write(lunerr1,'(/1X,17(''-''),'' forward '',17(''-''))')
        write(lunerr1,
     +      '('' EVOLVE iq1,2 = '',2I5,''   nf = '',I3)') iq1,iq2,nf0
      endif

!$OMP   PARALLEL DO PRIVATE(I,ID,IPTYP,NN)
      do id = 1,n
        iptyp = iptype(isns(id))
        nn    = nevols(isns(id))
        call sqcEvDglap(stor7,idWijk7(1,1,1,iptyp,ityp), idEijk7,
     +      idPdfi7(1,id), startu(1,1,id), 2, nn, iz1, iz2, eps(id))
      enddo
!$OMP   END PARALLEL DO
*      call sqcEvFixNf(ityp,mset,nf0,+1,startu,iz1f(nf0),iz2f(nf0),eps)

C-6   Evolution backward   ---------------------------------------------
      iz1   = iz1b(nf0)
      iz2   = iz2b(nf0)
      if(idbug6.ge.1) then
        iq1 = itfiz5(iz1)
        iq2 = itfiz5(iz2)
        write(lunerr1,'(1X,17(''-''),'' reverse '',17(''-''))')
        write(lunerr1,
     +      '('' EVOLVE iq1,2 = '',2I5,''   nf = '',I3)') iq1,iq2,nf0
      endif

!$OMP   PARALLEL DO PRIVATE(I,ID,IPTYP,NN)
      do id = 1,n
        iptyp = iptype(isns(id))
        nn    = nevols(isns(id))
        call sqcEvDglap(stor7,idWijk7(1,1,1,iptyp,ityp), idEijk7,
     +      idPdfi7(1,id), startd(1,1,id), 2, nn, iz1, iz2, eps(id))
      enddo
!$OMP   END PARALLEL DO
*      call sqcEvFixNf(ityp,mset,nf0,+1,startd,iz1b(nf0),iz2b(nf0),eps)

C-7   Finalise ---------------------------------------------------------
C--   Re-set cuts in iz (which may have been stretched to include it0)
      izmic2 = izfit5(-itmic2)
      izmac2 = izfit5( itmac2)
C--   Validate the pdfs
      if(isns(1).eq.1) then
        id1 = 0
      else
        id1 = 1
      endif
      do id = id1,n
        idglobal = iqcIdPdfLtoG(mset,id)
        call sqcValidate(stor7,idglobal)
      enddo
C--   Max deviation
      epsi = 0.D0
      do id = 1,n
        epsi = max(epsi,eps(id))
      enddo
C--   Max nf
      nfheavy = nfma

      if(idbug6.ne.0) write(lunerr1,'(/)')

      return
      end

C     ==========================================================
      subroutine sqcEvFixNf(ityp,mset,nf,mode,start,iz1,iz2,eps)
C     ==========================================================

C--   Evolve all pdfs at fixed nf
C--
C--   ityp            (in) : unpol/pol/timelike [1,3]
C--   mset            (in) : pdf set
C--   nf              (in) : number of active flavours
C--   mode            (in) : +1 input start value; -1 input discontinuity
C--   start(2,iy,id)  (in) : input startvalues or discontinuities
C--   iz1             (in) : start point
C--   iz2             (in) : end point
C--   eps(12)        (out) : max deviation

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension start(2,mxx0,*), eps(*)

C--   Evolve 2nf identifiers, fill the rest
      dimension idevol(12,3:6)
C--                 1   2   3   4   5   6   7   8   9  10  11  12
      data idevol / 1,  2,  3,  7,  8,  9,  4,  5,  6, 10, 11, 12, !nf=3
     +              1,  2,  3,  4,  7,  8,  9, 10,  5,  6, 11, 12, !nf=4
     +              1,  2,  3,  4,  5,  7,  8,  9, 10, 11,  6, 12, !nf=5
     +              1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12/ !nf=6

C--   Pdf type SG, NS+, NS-, VA; number of coupled evolutions
      dimension iptype(12),nevols(12)
C--                 1   2   3   4   5   6   7   8   9  10  11  12
      data iptype / 1,  2,  2,  2,  2,  2,  4,  3,  3,  3,  3,  3/
      data nevols / 2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1/

      if(idbug6.ge.1) then
        iq1 = itfiz5(iz1)
        iq2 = itfiz5(iz2)
        if(mode.eq.1) then
          write(lunerr1,
     +      '('' EVOLVE iq1,2 = '',2I5,''   nf = '',I3,''  start'')')
     +        iq1,iq2,nf
        else
          write(lunerr1,
     +      '('' EVOLVE iq1,2 = '',2I5,''   nf = '',I3)')
     +        iq1,iq2,nf
        endif
      endif

C--   Initialize
      do i = 1,12
        eps(i) = 0.D0
      enddo

C--   Evolution direction
      if(iz1.le.iz2) then
        istep = +1
      else
        istep = -1
      endif

C--   Evolve
!$OMP   PARALLEL DO PRIVATE(I,ID,IPTYP,NN)
      do i = 1,2*nf
        id    = idevol(i,nf)
        iptyp = iptype(id)
        nn    = mode * nevols(id)
*mb        write(6,'(''EVOLVE ID, NN, IZ1, IZ2 = '',4I5)') id,nn,iz1,iz2
        call sqcEvDglap(stor7,idWijk7(1,1,1,iptyp,ityp), idEijk7,
     +      idPdfi7(1,id), start(1,1,id), 2, nn, iz1, iz2, eps(id))
      enddo
!$OMP   END PARALLEL DO

C--   Fill inactive heavy quarks
!$OMP   PARALLEL DO PRIVATE(I,ID,IDG,IZ0,IZ)
      do i = 2*nf+1,12
        id  = idevol(i,nf)
        nn  = mode * nevols(id)
*mb        write(6,'(''FILL   ID, NN, IZ1, IZ2 = '',5I5)') id,nn,iz1,iz2,
*mb     +  idPdfi7(1,id)
        call sqcEvDfill(stor7,idPdfi7(1,id),start(1,1,id),2,nn,iz1,iz2)
      enddo
!$OMP   END PARALLEL DO

      idum = mset ! avoid compiler warning

      return
      end

C     ==================================================================
C     EVOLUTION UTILITY ROUTINES
C     ==================================================================

C     ===========================================================
      subroutine sqcSetStart(fun,tmatqf,nf,nfheavy,startu,startd)
C     ===========================================================

C--   Convert set of arbitrary pdfs |fj> to basis pdfs |ei>

C--   fun(j,x)        (in) : returns fj(x) j = 0,...2nfheavy
C--   tmatqf(13,13)   (in) : tmadqf(i,j) contribution of fj to ei
C--   nf              (in) : number of active flavours
C--   nfheavy         (in) : largest heavy flavour
C--   start(2,iy,id) (out) : start(1,iy,1) = singlet
C--                          start(2,iy,1) = gluon
C--                          start(1,iy,i) = ei  i = 2,..,12
C--                          start(2,iy,i) = not used

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'

      external  fun
      dimension tmatqf(13,13), pdf(13), epm(13)
      dimension startu(2,mxx0,*), startd(2,mxx0,*)

C--   Fill start array
      do iy = 1,iymac2
        y        = ygrid2(iy)
        x        = exp(-y)
C--     Get pdf input
        pdf(1) = fun(0,x)
        do id = 1,2*nfheavy
          pdf(id+1) = fun(id,x)
        enddo
        do id = 2*nfheavy+1,12
          pdf(id+1) = 0.D0
        enddo
C--     Transform to basis pdfs
        call sqcPDFtoEPM(tmatqf,pdf,epm,nf)
C--     Copy singlet/gluon into start
        do id = 0,1
          startu(2-id,iy,1) = epm(id+1)
          startd(2-id,iy,1) = epm(id+1)
        enddo
C--     Copy nonsinglets into start
        do id = 2,12
          startu(1,iy,id) = epm(id+1)
          startd(1,iy,id) = epm(id+1)
          startu(2,iy,id) = 0.D0
          startd(2,iy,id) = 0.D0
        enddo
      enddo

      return
      end

C     =================================================================
      subroutine sqcEvPlan(iz0,nf0,nfmi,nfma,iz1f,iz2f,iz1b,iz2b,ip,ie)
C     =================================================================

C--   Setup evolution limits

C--   iz0       (in) : starting point
C--   nf0      (out) : nf at iz0
C--   nfmi     (out) : lowest nf
C--   nfma     (out) : highest nf
C--   iz1f(nf) (out) : first zbin forward  evolution at fixed nf
C--   iz2f(nf) (out) : last  zbin forward  evolution at fixed nf
C--   iz1b(nf) (out) : first zbin backward evolution at fixed nf
C--   iz2b(nf) (out) : last  zbin backward evolution at fixed nf
C--   ip        (in) : 0/1 = no/yes debug printout
C--   ie       (out) : 0 = all OK
C--                    1 = iz0 not inside current evolution cuts
C--                    2 = at least one limit below alphas cut

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension iz1f(3:6), iz2f(3:6), iz1b(3:6), iz2b(3:6)

C--   Check iz0 is inside evolution cuts
      if(iz0.lt.izmic2 .or. iz0.gt.izmac2) then
        ie = 1
        if(ip.eq.1) write(6,*) 'sqcEvplan: iz0 not inside current cuts'
        return
      else
        ie = 0
      endif
C--   Initialise
      nf0  = itfiz5(-iz0   )
      nfmi = itfiz5(-izmic2)
      nfma = itfiz5(-izmac2)
      do nf = 3,6
        iz1f(nf) = 0
        iz2f(nf) = 0
        iz1b(nf) = 0
        iz2b(nf) = 0
      enddo
C--
C--   Forward: get limits and check if not below alphas cut
      do nf = nf0,nfma
        iz1f(nf) = max(iz0   ,iz15(nf))
        iz2f(nf) = min(izmac2,iz25(nf))
        it1      = itfiz5(iz1f(nf))
        if(it1.lt.itmin6) ie = 2
      enddo

C--   Backward: get limits and check if not below alphas cut
      do nf = nf0,nfmi,-1
        iz1b(nf) = min(iz0   ,iz25(nf))
        iz2b(nf) = max(izmic2,iz15(nf))
        it2      = itfiz5(iz2b(nf))
        if(it2.lt.itmin6) ie = 2
      enddo

      if(ip.eq.1) then
        write(6,'(/'' sqcEVPLAN evolution limits:''/)')
        do i = nf0,nfma
          write(6,'('' Forward   iz1,2 = '',2I3,'' with nf = '',I3)')
     +         iz1f(i),iz2f(i),i
        enddo
        do i = nf0,nfmi,-1
          write(6,'('' Reverse   iz1,2 = '',2I3,'' with nf = '',I3)')
     +         iz1b(i),iz2b(i),i
        enddo
      endif

      return
      end
