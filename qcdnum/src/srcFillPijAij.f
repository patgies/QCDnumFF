
C--   This is the file srcFillPijAij.f with Pij and Aij weight routines

C--   subroutine sqcIniWt

C--   subroutine sqcFilWt(filit,lun,ityp,nwlast,ierr)
C--   subroutine sqcFilWU(w,nw,idum,kset,nwords,idPij,idAijk,mxord,ierr)
C--   double precision function dqcTimesNf(iq,nf)
C--   subroutine sqcFilWP(w,nw,idum,kset,nwords,idPij,idAijk,mxord,ierr)
C--   subroutine sqcFilWF(w,nw,idum,kset,nwords,idPij,idAijk,mxord,ierr)

C--   subroutine sqcDumpWt(lun,ityp,key,ierr)
C--   subroutine sqcDumpPij(w,lun,ityp,idPij,idAijk,mxord,ierr)
C--   subroutine sqcReadWt(lun,key,nwlast,ityp,ierr)
C--   subroutine sqcReadPij(w,nw,lun,kset,lastw,idPij,idAijk,mxord,ierr)

C     ==================================================================
C     Weight filling routines ==========================================
C     ==================================================================

C     ===================      
      subroutine sqcIniWt
C     =================== 

C--   Initializes weight tables

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'

C--   Initialise table bookkeeping arrays
      do m = 1,mstp0
        isetp7(m)  = 0
        mxord7(m)  = 0
        do j = 1,mord0
          do i = 1,mpp0
            idPij7(i,j,m) = 0
          enddo
        enddo
        do k = 1,mord0+1
          do j = 1,maa0
            do i = 1,maa0
              idAijk7(i,j,k,m) = 0
            enddo
          enddo
        enddo
      enddo

C--   Another table bookkeeping array
      do n = 1,mstp0
        do m = 1,4
          do k = 1,mord0
            do j = 1,2
              do i = 1,2
                idWijk7(i,j,k,m,n) = 0
              enddo
            enddo
          enddo
        enddo
      enddo

C--   Flag weight tables initialised      
      Lwtini7 = .true.
      
      return
      end                                        
          
C     ===============================================
      subroutine sqcFilWt(filit,lun,ityp,nwlast,ierr)
C     ===============================================

C--   Fill weight tables and book associated pdf tables

C--   filit    (in)   subroutine declared external in the calling routine
C--   lun      (in)   logical unit number only used when reading
C--   ityp     (in)   1=unpol, 2=pol, 3=frag, 4=custom 
C--   nwlast   (out)  last  word used < 0 not enough space
C--   ierr     (out)   0 OK
C--                   -1 weights already exist do nothing
C--                   -2 not enough space
C--                   -3 iset count exceeded [1-mst0]
C--                   -4 disk read error

C--   Output index arrays in /qstor7/
C--
C--   idPij7(i,j,k) i=QQ,QG,GQ,GG,N+,N-,VA; j=iord; k=ityp
C--   idWijk7(i,j,k,m,n) i,j=Q,G; k=iord; m=SG,N+,N-,VA; n=ityp
C--   idAijk7(i,j,k,m) i,j=G,Q,H; k=iord; m=ityp

C--   The subroutine filit should have the following syntax
C--
C--    subroutine filit(w,nw,lun,kset,nwords,idPij,idAijk,mxord,ierr)
C--
C--    w              (in)  store, declared w(*) in filit
C--    nw             (in)  total size of w in words
C--    lun            (in)  not used except in disk reading routine
C--    kset          (out)  table set identifier in w
C--    nwords        (out)  number of words used < 0 not enough space
C--    idPij(7,3)    (out)  Pij(i,iord) i = QQ,QG,GQ,GG,N+,N-,NV
C--    idAijk(3,3,4) (out)  Aij(i,j,iord) i,j = G,Q,H, iord = 4 scratch
C--    mxord         (out)  max order supported by the tables [1-mord0]
C--    ierr          (out)   0 OK
C--                         -1 empty set (never occurs)
C--                         -2 not enough space
C--                         -3 iset count MST0 exceeded
C--                         -4 disk read error

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension idPij(mpp0,mord0),idAijk(maa0,maa0,mord0+1)

      external filit

C--   Nonzero isetp7 indicates that the job is already done
      if(isetp7(ityp).ne.0) then
        ierr = -1
        return
      endif
      
C--   Fill tables (filit is a generic s/r name passed as an argument)
      call filit(stor7,nwf0,lun,kset,nwlast,idPij,idAijk,mxord,ierr)
C--   Check mxord in range 1-3
      if(mxord.lt.1 .or. mxord.gt.3) then
        stop 'sqcFilWt : maxorder not in range [1-3] ---> STOP'
      elseif(ierr.eq.-1) then
        stop 'sqcFilWt : attempt to book an empty set of tables'
      elseif(ierr.ne.0)  then
        return
      endif

C--   Store table indices
      do k = 1,mord0
        do i = 1,mpp0
            idPij7(i,k,ityp) = idPij(i,k)
        enddo
      enddo
      do k = 1,mord0+1
        do j = 1,maa0
          do i = 1,maa0
            idAijk7(i,j,k,ityp) = idAijk(i,j,k)
          enddo
        enddo
      enddo
C--   Again store table indices
      do k = 1,mord0
        idWijk7(1,1,k,1,ityp) = idPij(1,k) !PQQ
        idWijk7(1,2,k,1,ityp) = idPij(2,k) !PQG
        idWijk7(2,1,k,1,ityp) = idPij(3,k) !PGQ
        idWijk7(2,2,k,1,ityp) = idPij(4,k) !PGG
        idWijk7(1,1,k,2,ityp) = idPij(5,k) !NS+
        idWijk7(1,1,k,3,ityp) = idPij(6,k) !NS-
        idWijk7(1,1,k,4,ityp) = idPij(7,k) !VA
      enddo
C--   Store table set identifier and perturbative order
      isetp7(ityp) = kset
      mxord7(ityp) = mxord

      return
      end

C     ==================================================================
      subroutine sqcFilWU(w,nw,idum,kset,nwords,idPij,idAijk,mxord,ierr)
C     ==================================================================

C--   Fill Pij tables unpolarised

C--   w           (in)   store
C--   nw          (in)   total size of w in words
C--   idum        (in)   not used
C--   kset        (out)  table set identifier in w
C--   nwords      (out)  number of words used < 0 not enough space
C--   idPij       (out)  list of Pij  table identifiers
C--   idAijk      (out)  list of Aijk table identifiers
C--   mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)   0 OK
C--                      -1 empty set
C--                      -2 not enough space
C--                      -3 iset count exceeded
C--                      -4 disk read error

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension idPij(mpp0,mord0),idAijk(maa0,maa0,mord0+1)
      dimension itypes(mtyp0)

      external dqcAchi
      external dqcPQQ0R, dqcPQQ0S, dqcPQQ0D               ! (1,1) = PQQ0
      external dqcPQG0A                                   ! (2,1) = PQG0
      external dqcPGQ0A                                   ! (3,1) = PGQ0
      external dqcPGG0A, dqcPGG0R, dqcPGG0S, dqcPGG0D     ! (4,1) = PGG0
      
      external dqcPPL1A, dqcPPL1B                         ! (5,2) = PPL1
      external dqcPMI1B                                   ! (6,2) = PMI1
      external dqcPQQ1A, dqcPQQ1B                         ! (1,2) = PQQ1
      external dqcPQG1A                                   ! (2,2) = PQG1
      external dqcPGQ1A                                   ! (3,2) = PGQ1
      external dqcPGG1A, dqcPGG1B                         ! (4,2) = PGG1
      
      external dqcPPL2A, dqcPPL2B, dqcPPL2D               ! (5,3) = PPL2
      external dqcPMI2A, dqcPMI2B, dqcPMI2D               ! (6,3) = PMI2
      external dqcPVA2A                                   ! (7,3) = PVA2
      external dqcPQQ2A                                   ! (1,3) = PQQ2
      external dqcPQG2A                                   ! (2,3) = PQG2
      external dqcPGQ2A                                   ! (3,3) = PGQ2
      external dqcPGG2A, dqcPGG2B, dqcPGG2D               ! (4,3) = PGG2

      external dqcA000D                                  !  delta kernel
      external dqcAHH1B                                          !  AHH1
      external dqcAGH1A                                          !  AGH1
      external dqcAGQ2A                                          !  AGQ2
      external dqcAGG2A, dqcAGG2B, dqcAGG2D                      !  AGG2
      external dqcAQQ2A, dqcAQQ2B, dqcAQQ2D                      !  AQQ2
      external dqcAHQ2A                                          !  AHQ2
      external dqcAHG2A, dqcAHG2D                                !  AHG2

      external dqcTimesNf

      jdum = idum !avoid compiler warning

C--   Initialize table indices
      do k = 1,mord0
        do i = 1,mpp0
            idPij(i,k) = 0
        enddo
      enddo
      do k = 1,mord0+1
        do j = 1,maa0
          do i = 1,maa0
            idAijk(i,j,k) = 0
          enddo
        enddo
      enddo
      do i = 1,mtyp0
        itypes(i) = 0
      enddo

C--   Max perturbative order
      mxord = 3

C--   Book weight tables in the internal store
      new       = 0
      npar      = 20
      nusr      = 0
C--   Number of Pij and Aij tables + Aij extra tables
      itypes(1) = 8  + 1
      itypes(2) = 17 + 1
      call sqcMakeTab(w,nw,itypes,npar,nusr,new,kset,nwords)
      if(kset.lt.0) then
        ierr = kset
        return
      else
        ierr     = 0
      endif

C--   Weight table offset
      id0   = 1000*kset
           
C--   LO Pij -----------------------------------------------------------
      write(lunerr1,'('' Pij LO'')')
      idPij(1,1)     = id0+201                                   !PQQ LO
      call sqcUwgtRS(w,idPij(1,1),dqcPQQ0R,dqcPQQ0S,dqcAchi,1,ie)
      call sqcUweitD(w,idPij(1,1),dqcPQQ0D,dqcAchi,ie)
      idPij(2,1)     = id0+202                                   !PQG LO
      call sqcUweitA(w,idPij(2,1),dqcPQG0A,dqcAchi,ie)
      idPij(3,1)     = id0+203                                   !PGQ LO
      call sqcUweitA(w,idPij(3,1),dqcPGQ0A,dqcAchi,ie)
      idPij(4,1)     = id0+204                                   !PGG LO
      call sqcUweitA(w,idPij(4,1),dqcPGG0A,dqcAchi,ie)
      call sqcUwgtRS(w,idPij(4,1),dqcPGG0R,dqcPGG0S,dqcAchi,1,ie)
      call sqcUweitD(w,idPij(4,1),dqcPGG0D,dqcAchi,ie)
      idPij(5,1)     = idPij(1,1)                                !PPL LO
      idPij(6,1)     = idPij(1,1)                                !PMI LO
      idPij(7,1)     = idPij(1,1)                                !PVA LO

C--   NLO Pij ----------------------------------------------------------
      write(lunerr1,'('' Pij NLO'')')
      idPij(5,2)     = id0+205                                  !PPL NLO
      call sqcUweitA(w,idPij(5,2),dqcPPL1A,dqcAchi,ie)
      call sqcUweitB(w,idPij(5,2),dqcPPL1B,dqcAchi,1,ie)        !1=delta
      idPij(6,2)     = id0+206                                  !PMI NLO
      idPij(7,2)     = idPij(6,2)                               !PVA NLO
      call sqcUweitB(w,idPij(7,2),dqcPMI1B,dqcAchi,1,ie)        !1=delta
      idPij(1,2)     = id0+207                                  !PQQ NLO
      call sqcUweitA(w,idPij(1,2),dqcPQQ1A,dqcAchi,ie)
      call sqcUweitB(w,idPij(1,2),dqcPQQ1B,dqcAchi,1,ie)        !1=delta
      idPij(2,2)     = id0+208                                  !PQG NLO
      call sqcUweitA(w,idPij(2,2),dqcPQG1A,dqcAchi,ie)
      idPij(3,2)     = id0+209                                  !PGQ NLO
      call sqcUweitA(w,idPij(3,2),dqcPGQ1A,dqcAchi,ie)
      idPij(4,2)     = id0+210                                  !PGG NLO
      call sqcUweitA(w,idPij(4,2),dqcPGG1A,dqcAchi,ie)
      call sqcUweitB(w,idPij(4,2),dqcPGG1B,dqcAchi,1,ie)        !1=delta

C--   NNLO Pij ---------------------------------------------------------
      write(lunerr1,'('' Pij NNLO'')')
      idPij(5,3)     = id0+211                                 !PPL NNLO
      call sqcUweitA(w,idPij(5,3),dqcPPL2A,dqcAchi,ie)
      call sqcUweitB(w,idPij(5,3),dqcPPL2B,dqcAchi,0,ie)      !0=nodelta
      call sqcUweitD(w,idPij(5,3),dqcPPL2D,dqcAchi,ie)
      idPij(6,3)     = id0+212                                 !PMI NNLO
      call sqcUweitA(w,idPij(6,3),dqcPMI2A,dqcAchi,ie)
      call sqcUweitB(w,idPij(6,3),dqcPMI2B,dqcAchi,0,ie)      !0=nodelta
      call sqcUweitD(w,idPij(6,3),dqcPMI2D,dqcAchi,ie)
      idPij(7,3)     = id0+213                                 !PVA NNLO
      call sqcCopyWt(w,idPij(6,3),w,idPij(7,3),0)
      call sqcUweitA(w,idPij(7,3),dqcPVA2A,dqcAchi,ie)
      idPij(1,3)     = id0+214                                 !PQQ NNLO
      call sqcCopyWt(w,idPij(5,3),w,idPij(1,3),0)
      call sqcUweitA(w,idPij(1,3),dqcPQQ2A,dqcAchi,ie)
      idPij(2,3)     = id0+215                                 !PQG NNLO
      call sqcUweitA(w,idPij(2,3),dqcPQG2A,dqcAchi,ie)
      idPij(3,3)     = id0+216                                 !PGQ NNLO
      call sqcUweitA(w,idPij(3,3),dqcPGQ2A,dqcAchi,ie)
      idPij(4,3)     = id0+217                                 !PGG NNLO
      call sqcUweitA(w,idPij(4,3),dqcPGG2A,dqcAchi,ie)
      call sqcUweitB(w,idPij(4,3),dqcPGG2B,dqcAchi,0,ie)      !0=nodelta
      call sqcUweitD(w,idPij(4,3),dqcPGG2D,dqcAchi,ie)

C--   LO Aij -----------------------------------------------------------
      write(lunerr1,'('' Aij LO'')')
      idAijk(1,1,1)  = id0+101                                   !AGG LO
      call sqcUweitD(w,idAijk(1,1,1),dqcA000D,dqcAchi,ie)
      idAijk(2,2,1)  = id0+101                                   !AQQ LO
      idAijk(3,3,1)  = id0+101                                   !AHH LO

C--   NLO Aij ----------------------------------------------------------
      write(lunerr1,'('' Aij NLO'')')
      idAijk(1,3,2)  = id0+102                                  !AGH NLO
      call sqcUweitA(w,idAijk(1,3,2),dqcAGH1A,dqcAchi,ie)
      idAijk(3,3,2)  = id0+103                                  !AHH NLO
      call sqcUweitB(w,idAijk(3,3,2),dqcAHH1B,dqcAchi,1,ie)     !1=delta

C--   NNLO Aij ---------------------------------------------------------
      write(lunerr1,'('' Aij NNLO'')')
      idAijk(1,2,3)  = id0+104                                 !AGQ NNLO
      call sqcUweitA(w,idAijk(1,2,3),dqcAGQ2A,dqcAchi,ie)
      idAijk(1,1,3)  = id0+105                                 !AGG NNLO
      call sqcUweitA(w,idAijk(1,1,3),dqcAGG2A,dqcAchi,ie)
      call sqcUweitB(w,idAijk(1,1,3),dqcAGG2B,dqcAchi,0,ie)   !0=nodelta
      call sqcUweitD(w,idAijk(1,1,3),dqcAGG2D,dqcAchi,ie)
      idAijk(2,2,3)  = id0+106                                 !AQQ NNLO
      call sqcUweitA(w,idAijk(2,2,3),dqcAQQ2A,dqcAchi,ie)
      call sqcUweitB(w,idAijk(2,2,3),dqcAQQ2B,dqcAchi,0,ie)   !0=nodelta
      call sqcUweitD(w,idAijk(2,2,3),dqcAQQ2D,dqcAchi,ie)
      idAijk(3,2,3)  = id0+107                                 !AHQ NNLO
      call sqcUweitA(w,idAijk(3,2,3),dqcAHQ2A,dqcAchi,ie)
      idAijk(3,1,3)  = id0+108                                 !AHG NNLO
      call sqcUweitA(w,idAijk(3,1,3),dqcAHG2A,dqcAchi,ie)
      call sqcUweitD(w,idAijk(3,1,3),dqcAHG2D,dqcAchi,ie)

C--   NNLO Aij extra tables --------------------------------------------
C--   Fill Aqq^+ = Aqq + Ahq
      k             = mord0+1
      idAijk(2,2,k) = id0+109
      call sqcCopyWt(w,idAijk(2,2,3),w,idAijk(2,2,k),0)
      call sqcCopyWt(w,idAijk(3,2,3),w,idAijk(2,2,k),1)
C--   Fill Ahq^+ = Aqq - nf*Ahq
      idAijk(3,2,k) = id0+218
      call sqcCopyWt(w,idAijk(2,2,3),w,idAijk(3,2,k),0)
      call sqcWtimesF(dqcTimesNf,w,idAijk(3,2,3),w,idAijk(3,2,k),-1)

      return
      end

C     ===========================================
      double precision function dqcTimesNf(iq,nf)
C     ===========================================

      implicit double precision (a-h,o-z)

      jq         = iq   !avoid compiler warning
      dqcTimesNf = nf

      return
      end

C     ==================================================================
      subroutine sqcFilWP(w,nw,idum,kset,nwords,idPij,idAijk,mxord,ierr)
C     ==================================================================

C--   Fill Pij tables polarised

C--   w           (in)   store
C--   nw          (in)   total size of w in words
C--   idum        (in)   not used
C--   kset        (out)  table set identifier in w
C--   nwords      (out)  number of words used < 0 not enough space
C--   idPij       (out)  list of Pij table identifiers
C--   idAijk      (out)  Aijk table identifiers
C--   mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)   0 OK
C--                      -1 empty set
C--                      -2 not enough space
C--                      -3 iset count exceeded
C--                      -4 disk read error

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension idPij(mpp0,mord0),idAijk(maa0,maa0,mord0+1)
      dimension itypes(mtyp0)

      external dqcAchi
      external dqcDPQQ0A, dqcDPQQ0B, dqcDPQQ0D           ! (1,1) = DPQQ0
      external dqcDPQG0A                                 ! (2,1) = DPQG0
      external dqcDPGQ0A                                 ! (3,1) = DPGQ0
      external dqcDPGG0A, dqcDPGG0B, dqcDPGG0D           ! (4,1) = DPGG0
      
      external dqcDPPL1B                                 ! (5,2) = DPPL1
      external dqcDPMI1A, dqcDPMI1B                      ! (6,2) = DPMI1
      external dqcDPQS1A                                 ! (1,2) = DPQS1
      external dqcDPQG1A                                 ! (2,2) = DPQG1
      external dqcDPGQ1A                                 ! (3,2) = DPGQ1
      external dqcDPGG1A,dqcDPGG1R,dqcDPGG1S,dqcDPGG1D   ! (4,2) = DPGG1

      external dqcA000D                                  !  delta kernel

      jdum = idum !avoid compiler warning

C--   Initialize table indices
      do k = 1,mord0
        do i = 1,mpp0
            idPij(i,k) = 0
        enddo
      enddo
      do k = 1,mord0+1
        do j = 1,maa0
          do i = 1,maa0
           idAijk(i,j,k) = 0
          enddo
        enddo
      enddo
      do i = 1,mtyp0
        itypes(i) = 0
      enddo

C--   Max perturbative order
      mxord = 2

C--   Book weight tables in the internal store
      new       = 0
      npar      = 20
      nusr      = 0
C--   Number of Pij and Aij tables
      itypes(1) = 3
      itypes(2) = 8
      call sqcMakeTab(w,nw,itypes,npar,nusr,new,kset,nwords)
      if(kset.lt.0) then
        ierr = kset
        return
      else
        ierr     = 0
      endif

C--   Weigt table offset
      id0   = 1000*kset

C--   LO Pij -----------------------------------------------------------
      write(lunerr1,'('' Pij LO'')')
      idPij(1,1)     = id0+101                                  !DPQQ LO
      call sqcUweitA(w,idPij(1,1),dqcDPQQ0A,dqcAchi,ie)
      call sqcUweitB(w,idPij(1,1),dqcDPQQ0B,dqcAchi,1,ie)       !1=delta
      call sqcUweitD(w,idPij(1,1),dqcDPQQ0D,dqcAchi,ie)
      idPij(2,1)     = id0+201                                  !DPQG LO
      call sqcUweitA(w,idPij(2,1),dqcDPQG0A,dqcAchi,ie)
      idPij(3,1)     = id0+102                                  !DPGQ LO
      call sqcUweitA(w,idPij(3,1),dqcDPGQ0A,dqcAchi,ie)
      idPij(4,1)     = id0+202                                  !DPGG LO
      call sqcUweitA(w,idPij(4,1),dqcDPGG0A,dqcAchi,ie)
      call sqcUweitB(w,idPij(4,1),dqcDPGG0B,dqcAchi,1,ie)       !1=delta
      call sqcUweitD(w,idPij(4,1),dqcDPGG0D,dqcAchi,ie)
      idPij(5,1)     = idPij(1,1)                               !DPPL LO
      idPij(6,1)     = idPij(1,1)                               !DPMI LO
      idPij(7,1)     = idPij(1,1)                               !DPVA LO
      
C--   NLO Pij ----------------------------------------------------------
      write(lunerr1,'('' Pij NLO'')')
      idPij(5,2)     = id0+203                                  !PPL NLO
      call sqcUweitB(w,idPij(5,2),dqcDPPL1B,dqcAchi,1,ie)       !1=delta
      idPij(6,2)     = id0+204                                  !PMI NLO
      call sqcUweitA(w,idPij(6,2),dqcDPMI1A,dqcAchi,ie)
      call sqcUweitB(w,idPij(6,2),dqcDPMI1B,dqcAchi,1,ie)       !1=delta
      idPij(7,2)     = idPij(6,2)                               !PVA NLO
      idPij(1,2)     = id0+205                                  !PQQ NLO
      call sqcCopyWt(w,idPij(5,2),w,idPij(1,2),0)
      call sqcUweitA(w,idPij(1,2),dqcDPQS1A,dqcAchi,ie)
      idPij(2,2)     = id0+206                                  !PQG NLO
      call sqcUweitA(w,idPij(2,2),dqcDPQG1A,dqcAchi,ie)
      idPij(3,2)     = id0+207                                  !PGQ NLO
      call sqcUweitA(w,idPij(3,2),dqcDPGQ1A,dqcAchi,ie)
      idPij(4,2)     = id0+208                                  !PGG NLO
      call sqcUweitA(w,idPij(4,2),dqcDPGG1A,dqcAchi,ie)
      call sqcUwgtRS(w,idPij(4,2),dqcDPGG1R,dqcDPGG1S,dqcAchi,1,ie)
      call sqcUweitD(w,idPij(4,2),dqcDPGG1D,dqcAchi,ie)

C--   LO Aij -----------------------------------------------------------
      write(lunerr1,'('' Aij LO'')')
      idAijk(1,1,1)  = id0+103                                   !AGG LO
      call sqcUweitD(w,idAijk(1,1,1),dqcA000D,dqcAchi,ie)
      idAijk(2,2,1)  = id0+103                                   !AQQ LO
      idAijk(3,3,1)  = id0+103                                   !AHH LO

      return
      end
      
C     ==================================================================
      subroutine sqcFilWF(w,nw,idum,kset,nwords,idPij,idAijk,mxord,ierr)
C     ==================================================================

C--   Fill Pij tables timelike  (fragmentation functions)

C--   w           (in)   store
C--   nw          (in)   total size of w in words
C--   idum        (in)   not used
C--   kset        (out)  table set identifier in w
C--   nwords      (out)  number of words used < 0 not enough space
C--   idPij       (out)  list of Pij table identifiers
C--   idAijk      (out)  Aijk table identifiers
C--   mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)   0 OK
C--                      -1 empty set
C--                      -2 not enough space
C--                      -3 iset count exceeded
C--                      -4 disk read error

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension idPij(mpp0,mord0),idAijk(maa0,maa0,mord0+1)
      dimension itypes(mtyp0)
      
      external dqcAchi
      external dqcPQQ0R, dqcPQQ0S, dqcPQQ0D               ! (1,1) = PQQ0
      external dqcTQG0A                                   ! (2,1) = PQG0
      external dqcTGQ0A                                   ! (3,1) = PGQ0
      external dqcPGG0A, dqcPGG0R, dqcPGG0S, dqcPGG0D     ! (4,1) = PGG0
      
      external dqcTPL1A, dqcTPL1B                         ! (5,2) = PPL1
      external dqcTMI1B                                   ! (6,2) = PMI1
      external dqcTQQ1A, dqcTQQ1B                         ! (1,2) = PQQ1
      external dqcTQG1A                                   ! (2,2) = PQG1
      external dqcTGQ1A                                   ! (3,2) = PGQ1
      external dqcTGG1A, dqcTGG1B                         ! (4,2) = PGG1

      external dqcTHG1A                                        ! AHG NLO
      external dqcA000D                                  !  delta kernel

      jdum = idum !avoid compiler warning

C--   Initialize table indices
      do k = 1,mord0
        do i = 1,mpp0
            idPij(i,k) = 0
        enddo
      enddo
      do k = 1,mord0+1
        do j = 1,maa0
          do i = 1,maa0
            idAijk(i,j,k) = 0
          enddo
        enddo
      enddo
      do i = 1,mtyp0
        itypes(i) = 0
      enddo

C--   Max perturbative order
      mxord = 2

C--   Book weight tables in the internal store
      new       = 0
      npar      = 20
      nusr      = 0
C--   Number of Pij and Aij tables
      itypes(1) = 4
      itypes(2) = 10
      call sqcMakeTab(w,nw,itypes,npar,nusr,new,kset,nwords)
      if(kset.lt.0) then
        ierr = kset
        return
      else
        ierr     = 0
      endif

C--   Weigt table offset
      id0   = 1000*kset

C--   LO Pij  ----------------------------------------------------------
      write(lunerr1,'('' Pij LO'')')
      idPij(1,1)     = id0+201                                   !PQQ LO
      call sqcUwgtRS(w,idPij(1,1),dqcPQQ0R,dqcPQQ0S,dqcAchi,1,ie)
      call sqcUweitD(w,idPij(1,1),dqcPQQ0D,dqcAchi,ie)
      idPij(2,1)     = id0+202                                   !PQG LO
      call sqcUweitA(w,idPij(2,1),dqcTQG0A,dqcAchi,ie)
      idPij(3,1)     = id0+203                                   !PGQ LO
      call sqcUweitA(w,idPij(3,1),dqcTGQ0A,dqcAchi,ie)
      idPij(4,1)     = id0+204                                   !PGG LO
      call sqcUweitA(w,idPij(4,1),dqcPGG0A,dqcAchi,ie)
      call sqcUwgtRS(w,idPij(4,1),dqcPGG0R,dqcPGG0S,dqcAchi,1,ie)
      call sqcUweitD(w,idPij(4,1),dqcPGG0D,dqcAchi,ie)
      idPij(5,1)     = idPij(1,1)                                !PPL LO
      idPij(6,1)     = idPij(1,1)                                !PMI LO
      idPij(7,1)     = idPij(1,1)                                !PVA LO

C--   NLO Pij  ---------------------------------------------------------
      write(lunerr1,'('' Pij NLO'')')
      idPij(5,2)     = id0+205                                  !PPL NLO
      call sqcUweitA(w,idPij(5,2),dqcTPL1A,dqcAchi,ie)
      call sqcUweitB(w,idPij(5,2),dqcTPL1B,dqcAchi,1,ie)        !1=delta
      idPij(6,2)     = id0+206                                  !PMI NLO
      idPij(7,2)     = id0+206                                  !PVA NLO
      call sqcUweitB(w,idPij(6,2),dqcTMI1B,dqcAchi,1,ie)        !1=delta
      idPij(1,2)     = id0+207                                  !PQQ NLO
      call sqcUweitA(w,idPij(1,2),dqcTQQ1A,dqcAchi,ie)
      call sqcUweitB(w,idPij(1,2),dqcTQQ1B,dqcAchi,1,ie)        !1=delta
      idPij(2,2)     = id0+208                                  !PQG NLO
      call sqcUweitA(w,idPij(2,2),dqcTQG1A,dqcAchi,ie)
      idPij(3,2)     = id0+209                                  !PGQ NLO
      call sqcUweitA(w,idPij(3,2),dqcTGQ1A,dqcAchi,ie)
      idPij(4,2)     = id0+210                                  !PGG NLO
      call sqcUweitA(w,idPij(4,2),dqcTGG1A,dqcAchi,ie)
      call sqcUweitB(w,idPij(4,2),dqcTGG1B,dqcAchi,1,ie)        !1=delta

C--   LO Aij -----------------------------------------------------------
      write(lunerr1,'('' Aij LO'')')
      idAijk(1,1,1)  = id0+101                                   !AGG LO
      call sqcUweitD(w,idAijk(1,1,1),dqcA000D,dqcAchi,ie)
      idAijk(2,2,1)  = id0+101                                   !AQQ LO
      idAijk(3,3,1)  = id0+101                                   !AHH LO

C--   NLO Aij ----------------------------------------------------------
      write(lunerr1,'('' Aij NLO'')')
      idAijk(3,1,2)  = id0+102 !AHG NLO
      call sqcUweitA(w,idAijk(3,1,2),dqcTHG1A,dqcAchi,ie)

      return
      end

C     ==================================================================
C     Disk read and write ==============================================
C     ==================================================================

C     =======================================
      subroutine sqcDumpWt(lun,ityp,key,ierr)
C     =======================================

C--   Write header info and a set of Pij tables to disk
C--
C--   lun    (in)  logical unit number
C--   ityp   (in)  1=unpol, 2=pol, 3=timelike, 4=custom
C--   key    (in)  key character string
C--   ierr   (out) 0 = OK
C--                1 = write error
C--                2 = ityp does not exist

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qvers1.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      dimension idPij(mpp0,mord0),idAijk(maa0,maa0,mord0+1)

      character*(*) key
      character*50  keyout

C--   Check if ityp exists
      ierr = 2
      kset = isetp7(ityp)
      if(kset.eq.0) return
C--   Reformat key
      call sqcSetKey(key,keyout)
C--   Dump QCDNUM version
      write(lun,err=500) cvers1, cdate1
C--   Dump key
      write(lun,err=500) keyout
C--   Dump ityp 1=unpol, 2=pol, 3=timelike, 4=custom
      write(lun,err=500) ityp
C--   Dump some array sizes
      write(lun,err=500) mxg0, mxx0, mqq0, mst0
      write(lun,err=500) mord0, mnf0, mbp0, mpp0, maa0, mtyp0, mchk0
C--   Dump relevant grid parameters
      write(lun,err=500) nyy2, nyg2, ioy2, dely2
      write(lun,err=500) ntt2
      write(lun,err=500) (tgrid2(i),i=1,ntt2)

C--   Copy table indices
      do k = 1,mord0
        do i = 1,mpp0
            idPij(i,k) = idPij7(i,k,ityp)
        enddo
      enddo
      do k = 1,mord0+1
        do j = 1,maa0
          do i = 1,maa0
            idAijk(i,j,k) = idAijk7(i,j,k,ityp)
          enddo
        enddo
      enddo
      mxord = mxord7(ityp)
C--   Now flush out the store
      call sqcDumpPij(stor7,lun,kset,idPij,idAijk,mxord,ierr)

      return
      
  500 continue
C--   Write error
      ierr = 1

      return
      end

C     =========================================================
      subroutine sqcDumpPij(w,lun,kset,idPij,idAijk,mxord,ierr)
C     =========================================================

C--   Write one set of Pij tables to disk

C--   w           (in)   store
C--   lun         (in)   logical unit number
C--   kset        (in)   table set number in w
C--   idPij       (in)   list of Pij table identifiers
C--   idAijk      (in)   list of Aij table identifiers
C--   mxord       (in)   maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)  0 = OK
C--                      1 = write error

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension w(*)
      dimension itypes(mtyp0)
      dimension idPij(mpp0,mord0),idAijk(maa0,maa0,mord0+1)

      ierr = 0

C--   Fill itypes array
      do ityp = 1,mtyp0
        itypes(ityp) = iqcSgnNumberOfTables(w,kset,ityp)
      enddo

C--   Number of words and position of kset
      npar   = iqcGetNumberOfParams(w,kset)     !number of params
      ifirst = iqcFirstWordOfSet(w,kset)        !first word of kset
      nwdset = iqcGetNumberOfWords(w(ifirst))   !number of words in kset
      ilast  = ifirst + nwdset - 1              !last word of set

C--   Control word
      icword = 123456

C--   Go...
      write(lun,err=500) icword
      write(lun,err=500) nwdset,itypes,npar,idPij,idAijk,mxord
      write(lun,err=500) (w(i),i=ifirst,ilast)

      return

  500 continue
C--   Write error
      ierr = 1

      return
      end

C     ==============================================
      subroutine sqcReadWt(lun,key,nwlast,ityp,ierr)
C     ==============================================

C--   Read weights from a disk file (unformatted read)
C--  
C--   lun          (in)   input logical unit number
C--   key          (in)   key character string
C--   nwlast       (out)  last word used in the store < 0 no space
C--   ityp         (out)  +-1=unpol,+-2=pol,+-3=timelike,+-4=custom
C--                       minus sign if set already exists (no read)
C--   ierr         (out)  0 = all OK
C--                       1 = read error
C--                       2 = problem with QCDNUM version
C--                       3 = key mismatch
C--                       4 = x-mu2 grid not the same
C--                       5 = not enough space
C--                       6 = iset count exceeded [1,mst0]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qvers1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      character*10  cversr
      character*8   cdater
      character*50  keyred
      character*(*) key
      logical lqcSjekey
      dimension nyyr(0:mxg0),delyr(0:mxg0)
      dimension tgridr(mqq0)

      external sqcReadPij

      nwlast  = 0
      ierr    = 0

      rewind(UNIT=lun)

C--   Do nothing except return the weight type
      if(key.eq.'GIVETYPE') then
        read(lun,err=500,end=500) cversr, cdater
        read(lun,err=500,end=500) keyred
        read(lun,err=500,end=500) ityp
        return
      endif

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
C--   Read ityp 1=unpol, 2=pol, 3=timelike, 4=custom
      read(lun,err=500,end=500) ityp
C--   Check if ityp exists
      kset = isetp7(ityp)
      if(kset.ne.0) then
        ityp = -ityp
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
C--   Now read the Pij set
      call sqcFilWt(sqcReadPij,lun,ityp,nwlast,jerr)
      if(jerr.eq.0) then
        ierr = 0
        return
      elseif(jerr.eq.-1) then
C--     Set already exists
        ierr = 0
        ityp = -ityp
        return
      elseif(jerr.eq.-2) then
C--     Not enough space
        ierr = 5
        return
      elseif(jerr.eq.-3) then
C--     Iset count exceeded
        ierr = 6
        return
      elseif(jerr.eq.-4) then
C--     Read error
        ierr = 1
        return
      else
        stop 'sqcReadWt: unknown error code from sqcFilWt'
      endif

      return
      
C--   Read error      
 500  continue     
      ierr = 1
      
      return
      end

C     ==================================================================
      subroutine sqcReadPij(w,nw,lun,kset,lastw,idPij,idAijk,mxord,ierr)
C     ==================================================================

C--   Fill Pij tables by reading them from disk

C--   w           (in)   store
C--   nw          (in)   total size of w in words
C--   lun         (in)   logical unit number
C--   kset        (out)  table set identifier in w
C--   lastw       (out)  last word used in the store < 0 no space
C--   idPij       (out)  list of Pij table identifiers
C--   idAijk      (out)  Aijk table identifiers
C--   mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)   0 OK
C--                      -1 empty set
C--                      -2 not enough space
C--                      -3 iset count exceeded
C--                      -4 disk read error
C--

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension itypes(mtyp0)
      dimension idPij(mpp0,mord0),idAijk(maa0,maa0,mord0+1)
      dimension idPpp(mpp0,mord0),idAaa(maa0,maa0,mord0+1)

      ierr = 0

C--   Read header information
      read(lun,err=500,end=500) icword
      if(icword.ne.123456) goto 500
      read(lun,err=500,end=500) nwords,itypes,npar,idPpp,idAaa,mxord

C--   Book weight tables in the internal store
      new       = 0
      nusr      = 0
      call sqcMakeTab(w,nw,itypes,npar,nusr,new,kset,lastw)
      if(kset.lt.0) then
        ierr = kset
        return
      endif

C--   Setup Pij and Aij arrays
      do k = 1,mord0
        do i = 1,mpp0
          idPij(i,k) = 1000*kset + iqcGetLocalId(idPpp(i,k))
        enddo
      enddo
      do k = 1,mord0+1
        do j = 1,maa0
          do i = 1,maa0
            idAijk(i,j,k) = 1000*kset + iqcGetLocalId(idAaa(i,j,k))
          enddo
        enddo
      enddo

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
      ierr = -4

      return
      end

