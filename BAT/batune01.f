C     ------------------------------------------------------------------
      program batune01
C     ------------------------------------------------------------------
C--   Qcdnum tuning for BAT - 01
C--   - Read input pdfs created by batune00
C--   - Evolve 100-30000 GeV2 on 300x150 grid
C--   - Write high density reference splines
C     ------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      include '../splint/inc/splint.inc'
      include '../splint/inc/spliws.inc'
      
      real tim1, tim2
      
      external func                                          !input pdfs
      dimension def(-6:6,12)                   !input flavor composition
      dimension xmin(5), iwt(5)
      data xmin/1.D-3,0.2D0,0.4D0,0.6D0,0.75D0/                  !x grid
      data iwt / 1, 2, 4, 8, 16 /                                !x grid
      data ngx/5/, nxin/300/, iosp/3/                            !x grid
      dimension qq(2),wt(2)                                    !mu2 grid
      data qq/1.D2,3.D4/,wt/2*1.D0/                            !mu2 grid
      data ngq/2/, nqin/150/                                   !mu2 grid

      dimension ia0(0:10)                              !spline addresses
      dimension pars(13)                            !evolution parmeters
      character*2 pnam(0:10)
      data pnam /'gl','dv','db','uv','ub','sv','sb','cv','cb','bv','bb'/

C--   lun = 6 stdout, -6 stdout w/out banner page
      lun    = 6
      lunout = abs(lun)
      call qcinit(lun,' ')                                  !init QCDNUM
      call ssp_SpInit(13)                                   !init SPLINT
C--   Read input gluon and get parameters
      ia0(0) = isp_SpRead('../splines/'//pnam(0)//'.spl')
      do i = 1,13
        pars(i) = dsp_SpGetVal(ia0(0),i)                !read parameters
      enddo
      qq(1) = dsp_SpGetVal(ia0(0),14)                        !read scale
C--   Read input quarks and get def
      do i = 1,10
        ia0(i) = isp_SpRead('../splines/'//pnam(i)//'.spl')
        do j = -6,6
          def(j,i) =  dsp_SpGetVal(ia0(i),j+7)       !get pdf definition
        enddo
      enddo
C--   Store sdresses of input splines
      do i = 0,10
        call ssp_Uwrite(i+1,dble(ia0(i)))
      enddo
C--   Setup QCDNUM
      call gxmake(xmin,iwt,ngx,nxin,nx,iosp)                     !x-grid
      call gqmake(qq,wt,ngq,nqin,nq)                           !mu2-grid
      call wtfile(1,'../weights/bat01.wgt')           !calculate weights
      call setord(int(pars(1)))                           !LO, NLO, NNLO
      call setalf(pars(2),pars(3))                         !input alphas
      nfin = int(pars(4))                                     !FFNS/VFNS
      iqc  = iqfrmq(pars(5))                            !charm threshold
      iqb  = iqfrmq(pars(6))                           !bottom threshold
      nfin = 5                                !re-set to FFNS 5 flavours
      call setcbt(nfin,iqc,iqb,999)              !thresholds in the VFNS
      iq0  = iqfrmq(qq(1))                                  !start scale
      iq1  = nq                                               !end scale
      qq1  = qfrmiq(iq1)                                      !end scale

C--   Print settings
      call PrSettings(lunout,def)

C--   Evolution
      call cpu_time(tim1)
      ntry = 1
      do i = 1,ntry
        call evolfg(1,func,def,iq0,eps)                          !evolve
      enddo
      call cpu_time(tim2)
      write(lunout,'(/'' Evol  : '',F10.2,'' ms'')')
     +                  1000*(tim2-tim1)/ntry

C--   Structure function weights
      call zmwfile('../weights/zmstf.wgt')

C--   Gluon, singlet and valence 2-dim reference splines
      call cpu_time(tim1)
      call refpdf('../splines/xg300x150.spl',0)
      call refpdf('../splines/xs300x150.spl',1)
      call refpdf('../splines/xv300x150.spl',7)
      call cpu_time(tim2)
      write(lunout,'(/'' RefPdf: '',F10.2,'' ms'')') 1000*(tim2-tim1)
C--   FL, F2 and xF3 2-dim reference splines
      call cpu_time(tim1)
      call refstf('../splines/FL300x150.spl',1)
      call refstf('../splines/F2300x150.spl',2)
      call refstf('../splines/F3300x150.spl',3)
      call cpu_time(tim2)
      write(lunout,'(/'' RefStf: '',F10.2,'' ms'')') 1000*(tim2-tim1)

C--   Plot gl + si vs x
      xmi = xmin(1)
      xma = 1.D0
      call pltvsx('../plots/bt01vsx1.dat',xmi,xma,qq(1),1.D0,1.D0)
      call pltvsx('../plots/bt01vsx2.dat',xmi,xma,qq(2),1.D0,1.D0)

C--   Plot gl + si vs q2
      qmi = qq(1)
      qma = qq(2)
      call pltvsq('../plots/bt01vsq1.dat',qmi,qma,xmin(1),1.D0,1.D0)
      call pltvsq('../plots/bt01vsq2.dat',qmi,qma,0.02D0 ,4.D0,3.D0)
      call pltvsq('../plots/bt01vsq3.dat',qmi,qma,0.50D0 ,1.D3,4.D1)
      call pltvsq('../plots/bt01vsq4.dat',qmi,qma,0.90D0 ,1.D0,1.D0)

      nw = isp_SpSize(0)
      write(lunout,'(/'' Memory size = '',I9)') nw
      nw = isp_SpSize(1)
      write(lunout,'( '' Words  used = '',I9)') nw
      nw = isp_SpSize(ia0(0))
      write(lunout,'( '' Spline size = '',I9)') nw

*      call ssp_Mprint

      end

C     ------------------------------------------------------------------
C     Reference splines ------------------------------------------------
C     ------------------------------------------------------------------

C     ============================
      subroutine refpdf(fnam,ipdf)
C     ============================

C--   Create, fill and dump 2-dim pdf reference spline

      implicit double precision (a-h,o-z)

      character fnam*(*)
      dimension pars(13)

      external sppdf
      
      ia  = isp_S2Make(3,3)
      call ssp_Uwrite(1,dble(ipdf))
      call ssp_S2Fill(ia,sppdf,0.D0)
      call cpypar(pars,13,0)
      do i = 1,13
        call ssp_SpSetVal(ia,i,pars(i))
      enddo
      q2 = qfrmiq(1)
      call ssp_SpSetVal(ia,14,q2)
      call ssp_SpDump(ia,fnam)
      call ssp_Erase(ia)

      write(6,'(/'' REFPDF: write file '',A)') fnam

      return
      end

C     ============================================
      double precision function sppdf(ix,iq,first)
C     ============================================

C--   Function to spline pdf

      implicit double precision(a-h,o-z)
      logical first
      save ipdf

      if(first) ipdf = int(dsp_Uread(1))

      sppdf = BVALIJ(1,ipdf,ix,iq,1)

      return
      end

C     ============================
      subroutine refstf(fnam,istf)
C     ============================

C--   Create, fill and dump 2-dim stf reference spline

      implicit double precision (a-h,o-z)

      character fnam*(*)
      dimension pars(13)

      dimension cpr(-6:6), pro(-6:6), val(-6:6)
      data cpr/ 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4./
      data val/-1.,-1.,-1.,-1.,-1.,-1., 0., 1., 1., 1., 1., 1., 1./

      do j = -6,6
        pro(j) = cpr(j)/9.D0
      enddo

      ia  = isp_S2Make(3,3)
      if(istf.eq.3) then
        call ssp_S2F123(ia,1,val,istf,0.D0)
      else
        call ssp_S2F123(ia,1,pro,istf,0.D0)
      endif
      call cpypar(pars,13,0)
      do i = 1,13
        call ssp_SpSetVal(ia,i,pars(i))
      enddo
      q2 = qfrmiq(1)
      call ssp_SpSetVal(ia,14,q2)
      call ssp_SpDump(ia,fnam)
      call ssp_Erase(ia)

      write(6,'(/'' REFSTF: write file '',A)') fnam

      return
      end

C     ------------------------------------------------------------------
C     Plot routines ----------------------------------------------------
C     ------------------------------------------------------------------

C     ============================================
      subroutine pltvsx(fnam,xmi,xma,q2,facg,facs)
C     ============================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfunx

      call ssp_Uwrite(1,q2)
      call ssp_Uwrite(2,facg)
      call ssp_Uwrite(3,facs)

      call ffplot(fnam,pfunx,2,xmi,xma,-300,' ')

      return
      end

C     ==========================================
      double precision function pfunx(i,z,first)
C     ==========================================

C--   gl/si vs x

      implicit double precision(a-h,o-z)
      logical first

      save qfix, factor
      dimension factor(0:1)

      if(first) then
        qfix      = dsp_Uread(1)
        factor(0) = dsp_Uread(2)
        factor(1) = dsp_Uread(3)
      endif

      pfunx = factor(i-1) * BVALXQ(1, i-1, z, qfix, 1)

      return
      end

C     ============================================
      subroutine pltvsq(fnam,qmi,qma,xx,facg,facs)
C     ============================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfunq

      call ssp_Uwrite(1,xx)
      call ssp_Uwrite(2,facg)
      call ssp_Uwrite(3,facs)

      call ffplot(fnam,pfunq,2,qmi,qma,-300,' ')

      return
      end

C     ==========================================
      double precision function pfunq(i,z,first)
C     ==========================================

C--   gl/si vs q

      implicit double precision(a-h,o-z)
      logical first

      save xfix, factor
      dimension factor(0:1)

      if(first) then
        xfix      = dsp_Uread(1)
        factor(0) = dsp_Uread(2)
        factor(1) = dsp_Uread(3)
      endif

      pfunq = factor(i-1) * BVALXQ(1, i-1, xfix, z, 1)

      return
      end

C     ------------------------------------------------------------------
C     Print QCDNUM settings --------------------------------------------
C     ------------------------------------------------------------------

C     ==================================
      subroutine PrSettings(lunout, def)
C     ==================================

      implicit double precision (a-h,o-z)

      dimension pars(13), def(-6:6,12)

      call CPYPAR(pars, 13, 0)
      call GRPARS(nx, xmi, xma, nq, qmi, qma, jord )
      write(lunout,'(/'' nx, xmi      :'',I5, E13.5)') nx, xmi
      write(lunout,'( '' nq, qmi, qma :'',I5,2E13.5)') nq, qmi, qma
      write(lunout,'( '' iord         :'',I5       )') int(pars(1))
      write(lunout,'( '' alfas, q2alf :'',5X,2E13.5)') pars(2), pars(3)
      write(lunout,'( '' nfix,  q2cbt :'',I5,3E13.5)') int(pars(4)),
     +                                        pars(5), pars(6), pars(7)
      write(lunout,'( '' iq0, q20     :'',I5, E13.5)') iqfrmq(qmi), qmi
      write(lunout,'(/,'' def  '',11(I3,1X))') (i,i=-5,5)
      do i = 1,10
        write(lunout,'(I4,2X,11F4.0)') i,(def(j,i),j=-5,5)
      enddo

      return
      end
      
C     ------------------------------------------------------------------
C     Evolution input pdfs ---------------------------------------------
C     ------------------------------------------------------------------

C     ======================================
      double precision function func(ipdf,x)
C     ======================================

      implicit double precision (a-h,o-z)

      save ia0
      dimension ia0(0:10)                              !spline addresses

      func = 0.D0
      if(ipdf.eq.-1) then
        do i = 0,10
          ia0(i) = int(dsp_Uread(i+1))
        enddo
      endif
      if(ipdf.eq. 0) func = dsp_FunS1(ia0( 0), x, 1)
      if(ipdf.eq. 1) func = dsp_FunS1(ia0( 1), x, 1)
      if(ipdf.eq. 2) func = dsp_FunS1(ia0( 2), x, 1)
      if(ipdf.eq. 3) func = dsp_FunS1(ia0( 3), x, 1)
      if(ipdf.eq. 4) func = dsp_FunS1(ia0( 4), x, 1)
      if(ipdf.eq. 5) func = dsp_FunS1(ia0( 5), x, 1)
      if(ipdf.eq. 6) func = dsp_FunS1(ia0( 6), x, 1)
      if(ipdf.eq. 7) func = dsp_FunS1(ia0( 7), x, 1)
      if(ipdf.eq. 8) func = dsp_FunS1(ia0( 8), x, 1)
      if(ipdf.eq. 9) func = dsp_FunS1(ia0( 9), x, 1)
      if(ipdf.eq.10) func = dsp_FunS1(ia0(10), x, 1)
      if(ipdf.eq.11) func = 0.D0
      if(ipdf.eq.12) func = 0.D0

      return
      end








