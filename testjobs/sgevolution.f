C     ------------------------------------------------------------------
      program SGevolution
C     ------------------------------------------------------------------

C--   Evolve the singlet-gluon using the nxn evolution toolbox
C--   Compare to the standard EVOLFG evolution

C--   You can do FFNS evolution upto NNLO with this example job
C--   For VFNS you must code the threshold loop yourself, see write-up
C--   Forward evolution; note that EVDGLAP can also do backward

C--   ------------------------------------------------------------------
C--   Declarations
C--   ------------------------------------------------------------------

      implicit double precision (a-h,o-z)

C--   Evolution parameters
      data iord /3/, nfin/3/, iosp/3/               !Order, VFNS, Spline
      data q2c/5.D0/, q2b/25.D0/, q0/2.0/               !thresh and mu20
      data amu2/1.0D0/, bmu2/0.0D0/                         !renor scale
      data as0/0.364/, r20/2.D0/                                 !alphas

C--   Grid parameters
      dimension xmin(5), iwt(5)
      data xmin/1.D-5,0.2D0,0.4D0,0.6D0,0.75D0/                  !x grid
      data iwt / 1, 2, 4, 8, 16 /                                !x grid
      data nxin/100/                                             !x grid
      dimension qq(2),wt(2)                                    !mu2 grid
      data qq/2.D0,1.D4/, wt/1.D0,1.D0/, nqin/50/              !mu2 grid

C--   Input pdf definition
      external func                                  !input parton dists
      dimension def(-6:6,12)                         !flavor composition
      data def  /
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,         !dval
     + 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,         !uval
     + 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,         !sval
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,         !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,         !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !sbar
     + 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,         !cval
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !cbar
     + 52*0.    /
      
C--   Pdf output
      data ichk/1/                                  !yes/no check limits

C--   ------------------------------------------------------------------
C--   Declarations for the nxn evolution toolbox
C--   ------------------------------------------------------------------
      parameter (nstoru = 102219)                         !size of store
      dimension storu(nstoru)                               !local store
      dimension idPiju(7,3), idAlfa(3)                !identifier arrays
      dimension idw(2,2,3),idf(2),ida(2,2,3)          !identifier arrays
      dimension itypes(6)                                   !table types
      dimension iqlim(2)                               !evolution limits
      dimension start(2,200)                                !startvalues
      data itypes/6*0/                                 !initialise types
      data ioweit/3/                                !LO,NLO,NNLO weights

C--   Routines that set alphas values
      external AsVal1, AsVal2, AsVal3

C--   As it stands, this example program cannot do VFNS
      if(nfin.eq.0) stop 'No VFNS in this example job'

C--   ------------------------------------------------------------------
C--   First do the standard evolution
C--   ------------------------------------------------------------------

C--   lun = 6 stdout, -6 stdout w/out banner page
      lun = 6
      call qcinit(lun,' ')

C--   Make x-mu2 grid
      call gxmake(xmin,iwt,5,nxin,nx,iosp)                       !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid

C--   Standard QCDNUM weights for unpolarised evolution
      call wtfile(1,'../weights/unpolarised.wgt')

C--   Define renormalisation scale
      call setabr(amu2,bmu2)                                !renor scale

C--   Order should not exceed that of the nxn weight calculation
      if(iord.gt.ioweit) stop 'Evolution order iord too large'

C--   Set evolution parameters
      call setord(iord)                                   !LO, NLO, NNLO
      call setalf(as0,r20)                                 !input alphas
      iqc  = iqfrmq(q2c)                                !charm threshold
      iqb  = iqfrmq(q2b)                               !bottom threshold
      call setcbt(nfin,iqc,iqb,999)                           !FFNS/VFNS
      iq0  = iqfrmq(q0)                                  !starting scale
      
C--   Do here the standard evolution
      call evolfg(1,func,def,iq0,eps)                  !evolve all pdf's

C--   ------------------------------------------------------------------
C--   Do SG evolution with the evolution toolbox
C--   ------------------------------------------------------------------

C--   Weight tables (may not need them all)
      itypes(1) = 5
      itypes(2) = 17
C--   Put 4 pdf and 3 alphas tables in the store (may not need them all)
      itypes(5) = 4
      itypes(6) = 3
C--   Isetw is the table set identifier assigned by QCDNUM
      new  = 0
      npar = 0
      call MakeTab(storu,nstoru,itypes,npar,new,isetw,nwordsu)
C--   Calculate  evolution weigths upto order ioweit
      call FilWT(storu,abs(lun),isetw,idpiju,ioweit)

C--   Fill tables of alphas values upto (as/2pi)^3
      idAlfa(1) = 1000*isetw+601
      idAlfa(2) = 1000*isetw+602
      idAlfa(3) = 1000*isetw+603
      call EvFillA(storu,idAlfa(1),AsVal1)                           !LO
      call EvFillA(storu,idAlfa(2),AsVal2)                          !NLO
      call EvFillA(storu,idAlfa(3),AsVal3)                         !NNLO

C--   Setup the identifiers for si/gl evolution
C--   Weight table ityp 1=qq, 2=qg, 3=gq, 4=gg, 5=ns+, 6=ns-, 7=nsv
      ityp = 0
      do i = 1,2
        do j = 1,2
          ityp = ityp+1
          do k = 1,iord
            idw(i,j,k) = idPiju(ityp,k)
            ida(i,j,k) = idAlfa(k)
          enddo
        enddo  
      enddo
C--   Singlet and gluon table identifiers
      idf(1) = 1000*isetw+501                                   !singlet
      idf(2) = 1000*isetw+502                                     !gluon

C--   Fill startvalues; take them from the standard evolution
      do ix = 1,nx
        start(1,ix) = bvalij(1,1,ix,iq0,1)               !singlet(ix,iq)
        start(2,ix) = bvalij(1,0,ix,iq0,1)                 !gluon(ix,iq)
      enddo

C--   Evolve
      iqlim(1) = iq0
      iqlim(2) = nq
      call EvDglap(storu,idw,ida,idf,start,2,2,iqlim,nf,eps)

C--   Compare
      iq2 = nq
      call sqcNNcompa(storu,idf(1),1,iq2,difs)
      call sqcNNcompa(storu,idf(2),0,iq2,difg)
      write(lun,'(/'' si/gl      dif = '', 2E15.5)') difs,difg
      
      end

C--   ------------------------------------------------------------------
C--   Compares evolution toolbox result with standard evolution
C--   ------------------------------------------------------------------

C     ========================================
      subroutine sqcNNcompa(ww,idw,idf,iq,dif)
C     ========================================

      implicit double precision (a-h,o-z)

      dimension ww(*)

      call grpars(nx,xmi,xma,nq,qmi,qma,iosp)

      dif = 0.D0
      do ix = 1,nx
        fff = bvalij(1,idf,ix,iq,1)
        uuu = EvPdfij(ww,idw,ix,iq,1)
        ddd = abs(uuu-fff)
        dif = max(dif,ddd)
      enddo
      
      return
      end

C--   ------------------------------------------------------------------
C--   External routines needed by the evolution toolbox
C--   ------------------------------------------------------------------

C     ===============================================
      double precision function AsVal1(iq,nf,ithresh)
C     ===============================================

C--   Get alpha_s/2pi

      implicit double precision (a-h,o-z)

      mf     = nf                !avoid compiler warning
      AsVal1 = 0.D0
      if(ithresh .ge. 0) then
        AsVal1 = AlTabN(0, iq,1,ierr)
      elseif(ithresh .eq. -1) then
        AsVal1 = AlTabN(0,-iq,1,ierr)
      else
        stop 'AsVal1: wrong ithresh'
      endif

      return
      end

C     ===============================================
      double precision function AsVal2(iq,nf,ithresh)
C     ===============================================

C--   Get (alpha_s/2pi)^2

      implicit double precision (a-h,o-z)

      mf     = nf                !avoid compiler warning
      AsVal2 = 0.D0
      if(ithresh .ge. 0) then
        AsVal2 = AlTabN(0, iq,2,ierr)
      elseif(ithresh .eq. -1) then
        AsVal2 = AlTabN(0,-iq,2,ierr)
      else
        stop 'AsVal2: wrong ithresh'
      endif

      return
      end

C     ===============================================
      double precision function AsVal3(iq,nf,ithresh)
C     ===============================================

C--   Get (alpha_s/2pi)^3

      implicit double precision (a-h,o-z)

      mf     = nf                !avoid compiler warning
      AsVal3 = 0.D0
      if(ithresh .ge. 0) then
        AsVal3 = AlTabN(0, iq,3,ierr)
      elseif(ithresh .eq. -1) then
        AsVal3 = AlTabN(0,-iq,3,ierr)
      else
        stop 'AsVal3: wrong ithresh'
      endif

      return
      end

C--   ------------------------------------------------------------------
C--   Weight filling routine
C--   ------------------------------------------------------------------

C     ========================================
      subroutine FilWT(w,lun,jset,idpij,mxord)
C     ========================================

C--   Fill Pij tables unpolarised upto NNLO.
C--   The splitting functions are taken from the QCDNUM library

C--   w           (in)   store previously partitioned by maketab
C--   lun         (in)   output stream for messages
C--   jset        (in)   set identifier
C--   idpij(7,3)  (out)  list of Pij table identifiers (global format)
C--   mxord       (in)   maximum perturbative order LO,NLO,NNLO
C--
C--   First index of idpij(7,3):  1=QQ 2=QG 3=GQ 4=GG 5=NS+ 6=NS- 7=VA
C--
C--   Warning when taking splitting functions from QCDNUM:
C--   P_QG has a factor 2n_f inside and this you may not want

      implicit double precision (a-h,o-z)
      
      dimension w(*)        
      dimension idpij(7,3)

C--   These routies are all in the QCDNUM library
      external dqcAchi
      external dqcPQQ0R, dqcPQQ0S, dqcPQQ0D            ! (1,1) = PQQ0
      external dqcPQG0A                                ! (2,1) = PQG0
      external dqcPGQ0A                                ! (3,1) = PGQ0
      external dqcPGG0A, dqcPGG0R, dqcPGG0S, dqcPGG0D  ! (4,1) = PGG0
      
      external dqcPPL1A, dqcPPL1B                      ! (5,2) = PPL1
      external dqcPMI1B                                ! (6,2) = PMI1
      external dqcPQQ1A, dqcPQQ1B                      ! (1,2) = PQQ1
      external dqcPQG1A                                ! (2,2) = PQG1
      external dqcPGQ1A                                ! (3,2) = PGQ1
      external dqcPGG1A, dqcPGG1B                      ! (4,2) = PGG1 
      
      external dqcPPL2A, dqcPPL2B, dqcPPL2D            ! (5,3) = PPL2
      external dqcPMI2A, dqcPMI2B, dqcPMI2D            ! (6,3) = PMI2
      external dqcPVA2A                                ! (7,3) = PVA2
      external dqcPQQ2A                                ! (1,3) = PQQ2
      external dqcPQG2A                                ! (2,3) = PQG2
      external dqcPGQ2A                                ! (3,3) = PGQ2
      external dqcPGG2A, dqcPGG2B, dqcPGG2D            ! (4,3) = PGG2
      
      external dqcAGQ2A                                !  AGQ2
      external dqcAGG2A, dqcAGG2B, dqcAGG2D            !  AGG2
      external dqcAQQ2A, dqcAQQ2B, dqcAQQ2D            !  AQQ2
      external dqcAHQ2A                                !  AHQ2
      external dqcAHG2A, dqcAHG2D                      !  AHG2

      write(lun,'(/'' FilWT fill user weight tables'')')

C--   Loop over spline order
C--   Pij LO
      write(lun,'('' FilWt Pij LO '')')
      idPij(1,1)   = 1000*jset+201 !PQQ LO
      call MakeWRS(w,1000*jset+201,dqcPQQ0R,dqcPQQ0S,dqcAchi,0) !0=delta
      call MakeWtD(w,1000*jset+201,dqcPQQ0D,dqcAchi)
      idPij(2,1)   = 1000*jset+202 !PQG LO
      call MakeWtA(w,1000*jset+202,dqcPQG0A,dqcAchi)
      idPij(3,1)   = 1000*jset+203 !PGQ LO
      call MakeWtA(w,1000*jset+203,dqcPGQ0A,dqcAchi)
      idPij(4,1)   = 1000*jset+204 !PGG LO
      call MakeWtA(w,1000*jset+204,dqcPGG0A,dqcAchi)
      call MakeWRS(w,1000*jset+204,dqcPGG0R,dqcPGG0S,dqcAchi,0) !0=delta
      call MakeWtD(w,1000*jset+204,dqcPGG0D,dqcAchi)
      idPij(5,1)   = 1000*jset+201 !PPL LO
      idPij(6,1)   = 1000*jset+201 !PMI LO
      idPij(7,1)   = 1000*jset+201 !PVA LO

      if(mxord.le.1) return

C--   Pij NLO
      write(lun,'('' FilWt Pij NLO'')')
      idPij(5,2)   = 1000*jset+205 !PPL NLO
      call MakeWtA(w,1000*jset+205,dqcPPL1A,dqcAchi)
      call MakeWtB(w,1000*jset+205,dqcPPL1B,dqcAchi,0)          !0=delta
      idPij(6,2)   = 1000*jset+206 !PMI NLO
      idPij(7,2)   = 1000*jset+206 !PVA NLO
      call MakeWtB(w,1000*jset+206,dqcPMI1B,dqcAchi,0)          !0=delta
      idPij(1,2)   = 1000*jset+207 !PQQ NLO
      call MakeWtA(w,1000*jset+207,dqcPQQ1A,dqcAchi)
      call MakeWtB(w,1000*jset+207,dqcPQQ1B,dqcAchi,0)          !0=delta
      idPij(2,2)   = 1000*jset+208 !PQG NLO
      call MakeWtA(w,1000*jset+208,dqcPQG1A,dqcAchi)
      idPij(3,2)   = 1000*jset+209 !PGQ NLO
      call MakeWtA(w,1000*jset+209,dqcPGQ1A,dqcAchi)
      idPij(4,2)   = 1000*jset+210 !PGG NLO
      call MakeWtA(w,1000*jset+210,dqcPGG1A,dqcAchi)
      call MakeWtB(w,1000*jset+210,dqcPGG1B,dqcAchi,0)          !0=delta

      if(mxord.le.2) return

C--   Pij NNLO
      write(lun,'('' FilWt Pij NNLO'')')
      idPij(5,3)   = 1000*jset+211 !PPL NNLO
      call MakeWtA(w,1000*jset+211,dqcPPL2A,dqcAchi)
      call MakeWtB(w,1000*jset+211,dqcPPL2B,dqcAchi,1)        !1=nodelta
      call MakeWtD(w,1000*jset+211,dqcPPL2D,dqcAchi)
      idPij(6,3)   = 1000*jset+212 !PMI NNLO
      call MakeWtA(w,1000*jset+212,dqcPMI2A,dqcAchi)
      call MakeWtB(w,1000*jset+212,dqcPMI2B,dqcAchi,1)        !1=nodelta
      call MakeWtD(w,1000*jset+212,dqcPMI2D,dqcAchi)
      idPij(7,3)   = 1000*jset+213 !PVA NNLO
      call CopyWgt(w,1000*jset+212,1000*jset+213,0)
      call MakeWtA(w,1000*jset+213,dqcPVA2A,dqcAchi)
      idPij(1,3)   = 1000*jset+214 !PQQ NNLO
      call CopyWgt(w,1000*jset+211,1000*jset+214,0)
      call MakeWtA(w,1000*jset+214,dqcPQQ2A,dqcAchi)
      idPij(2,3)   = 1000*jset+215 !PQG NNLO
      call MakeWtA(w,1000*jset+215,dqcPQG2A,dqcAchi)
      idPij(3,3)   = 1000*jset+216 !PGQ NNLO            
      call MakeWtA(w,1000*jset+216,dqcPGQ2A,dqcAchi)
      idPij(4,3)   = 1000*jset+217 !PGG NNLO
      call MakeWtA(w,1000*jset+217,dqcPGG2A,dqcAchi)
      call MakeWtB(w,1000*jset+217,dqcPGG2B,dqcAchi,1)        !1=nodelta
      call MakeWtD(w,1000*jset+217,dqcPGG2D,dqcAchi)
C--   Aij NNLO      
      write(lun,'('' FilWt Aij NNLO'')')
      call MakeWtA(w,1000*jset+101,dqcAGQ2A,dqcAchi)
      call MakeWtA(w,1000*jset+102,dqcAGG2A,dqcAchi)
      call MakeWtB(w,1000*jset+102,dqcAGG2B,dqcAchi,1)        !1=nodelta
      call MakeWtD(w,1000*jset+102,dqcAGG2D,dqcAchi)
      call MakeWtA(w,1000*jset+103,dqcAQQ2A,dqcAchi)
      call MakeWtB(w,1000*jset+103,dqcAQQ2B,dqcAchi,1)        !1=nodelta
      call MakeWtD(w,1000*jset+103,dqcAQQ2D,dqcAchi)
      call MakeWtA(w,1000*jset+104,dqcAHQ2A,dqcAchi)
      call MakeWtA(w,1000*jset+105,dqcAHG2A,dqcAchi)
      call MakeWtD(w,1000*jset+105,dqcAHG2D,dqcAchi)

      return
      end

C--   ------------------------------------------------------------------
C--   Pdf input to the standard evolution
C--   ------------------------------------------------------------------

C     ======================================
      double precision function func(ipdf,x)
C     ======================================

      implicit double precision (a-h,o-z)

                     func = 0.D0
      if(ipdf.eq. 0) func = xglu(x)
      if(ipdf.eq. 1) func = xdnv(x)
      if(ipdf.eq. 2) func = xupv(x)
      if(ipdf.eq. 3) func = 0.D0
      if(ipdf.eq. 4) func = xdbar(x)
      if(ipdf.eq. 5) func = xubar(x)
      if(ipdf.eq. 6) func = xsbar(x)
      if(ipdf.eq. 7) func = 0.D0
      if(ipdf.eq. 8) func = 0.D0
      if(ipdf.eq. 9) func = 0.D0
      if(ipdf.eq.10) func = 0.D0
      if(ipdf.eq.11) func = 0.D0
      if(ipdf.eq.12) func = 0.D0

      return
      end
 
C     =================================
      double precision function xupv(x)
C     =================================

      implicit double precision (a-h,o-z)

      data au /5.107200D0/

      xupv = au * x**0.8D0 * (1.D0-x)**3.D0

      return
      end

C     =================================
      double precision function xdnv(x)
C     =================================

      implicit double precision (a-h,o-z)

      data ad /3.064320D0/

      xdnv = ad * x**0.8D0 * (1.D0-x)**4.D0

      return
      end

C     =================================
      double precision function xglu(x)
C     =================================

      implicit double precision (a-h,o-z)

      data ag    /1.7D0 /

      common /msum/ glsum, uvsum, dvsum

      xglu = ag * x**(-0.1D0) * (1.D0-x)**5.D0

      return
      end

C     ==================================
      double precision function xdbar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      data adbar /0.1939875D0/

      xdbar = adbar * x**(-0.1D0) * (1.D0-x)**6.D0

      return
      end

C     ==================================
      double precision function xubar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      xubar = xdbar(x) * (1.D0-x)

      return
      end

C     ==================================
      double precision function xsbar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      xsbar = 0.2D0 * (xdbar(x)+xubar(x))

      return
      end







