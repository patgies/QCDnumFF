C     ------------------------------------------------------------------
      program batune02
C     ------------------------------------------------------------------
C--   Qcdnum tuning for BAT - 02
C--   - Read input pdfs created by BATUNE00
C--   - Read gluon, singlet and valence reference splines (BATUNE01)
C--   - Tune QCDNUM grid with respect to reference splines
C     ------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      include '../splint/inc/splint.inc'
      include '../splint/inc/spliws.inc'
      
      real tim1, tim2
      
      external func                                          !input pdfs
      dimension def(-6:6,12)                   !input flavor composition
      dimension xmin(5), iwt(5)
      data xmin/1.D-3,0.2D0,0.4D0,0.6D0,0.8D0/                   !x grid
      data iwt / 1, 2, 4, 8, 16 /                                !x grid
      data ngx/5/, nxin/100/, iosp/3/                            !x grid
      dimension qq(3),wt(3)                                    !mu2 grid
      data qq/1.D2,1.D3,3.D4/,wt/3*1.D0/                       !mu2 grid
      data ngq/3/, nqin/50/                                    !mu2 grid

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
      ngx = 5
      call gxmake(xmin,iwt,ngx,nxin,nx,iosp)                     !x-grid
      call gqmake(qq,wt,ngq,nqin,nq)                           !mu2-grid
      call wtfile(1,'../weights/bat02.wgt')           !calculate weights
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

C--   Read pdf 2-dim reference splines
      iag = isp_SpRead('../splines/xg300x150.spl')
      ias = isp_SpRead('../splines/xs300x150.spl')
      iav = isp_SpRead('../splines/xv300x150.spl')

C--   Plot comparison gl, si, and va vs x
      xmi  = xmin(1)
*      xmi  = 0.1D0
      xma  = 0.9D0
      qfix = qq(3)
      call pltvsx('../plots/bt02gx100x50.dat',xmi,xma,qfix,0,iag)
      call pltvsx('../plots/bt02sx100x50.dat',xmi,xma,qfix,1,ias)
      call pltvsx('../plots/bt02vx100x50.dat',xmi,xma,qfix,7,iav)

C--   Plot gl + si vs q2
*      qmi  = qq(1)
*      qma  = qq(2)
*      xfix = xmin(1)
*      call pltvsq('../plots/bt02fq200.dat',qmi,qma,-xfix,iag,ias)

      nw = isp_SpSize(0)
      write(lunout,'(/'' Memory size = '',I9)') nw
      nw = isp_SpSize(1)
      write(lunout,'( '' Words  used = '',I9)') nw
      nw = isp_SpSize(ia0(0))
      write(lunout,'( '' Spline size = '',I9)') nw

*      call ssp_Mprint

      end

C     ------------------------------------------------------------------
C     Plot routines ----------------------------------------------------
C     ------------------------------------------------------------------

C     ==========================================
      subroutine pltvsx(fnam,xmi,xma,q2,ipdf,ia)
C     ==========================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfunx

      call ssp_Uwrite(1,q2)
      call ssp_Uwrite(2,dble(ipdf))
      call ssp_Uwrite(3,dble(ia))

      call ffplot(fnam,pfunx,4,xmi,xma,-300,' ')

      return
      end

C     ==========================================
      double precision function pfunx(i,z,first)
C     ==========================================

C--   plot vs x: 1=pdf, 2=ref, 3=pdf-ref, 4=(pdf-ref)/ref

      implicit double precision(a-h,o-z)
      logical first

      save qfix, ipdf, ia

      if(first) then
        qfix = dsp_Uread(1)
        ipdf = int(dsp_Uread(2))
        ia   = int(dsp_Uread(3))
      endif

      pdf   = BVALXQ(1, ipdf, z, qfix, 1)
      ref   = dsp_FunS2(ia, z, qfix, 1)
      if(i.eq.1) then
        pfunx = pdf
      elseif(i.eq.2) then
        pfunx = ref
      elseif(i.eq.3) then
        pfunx = pdf-ref
      elseif(i.eq.4) then
        pfunx = (pdf-ref)/ref
      else
        pfunx = 0.D0
      endif

      return
      end

C     ==========================================
      subroutine pltvsq(fnam,qmi,qma,xx,ipdf,ia)
C     ==========================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfunq

      call ssp_Uwrite(1,xx)
      call ssp_Uwrite(2,dble(ipdf))
      call ssp_Uwrite(3,dble(ia))

      call ffplot(fnam,pfunq,4,qmi,qma,-300,' ')

      return
      end

C     ==========================================
      double precision function pfunq(i,z,first)
C     ==========================================

C--   plot vs q: 1=pdf, 2=ref, 3=pdf-ref, 4=(pdf-ref)/ref

      implicit double precision(a-h,o-z)
      logical first

      save xfix, ipdf, ia

      if(first) then
        xfix = dsp_Uread(1)
        ipdf = int(dsp_Uread(2))
        ia   = int(dsp_Uread(3))
      endif

      pdf   = BVALXQ(1, ipdf, xfix, z, 1)
      ref   = dsp_FunS2(ia, xfix, z, 1)
      if(i.eq.1) then
        pfunq = pdf
      elseif(i.eq.2) then
        pfunq = ref
      elseif(i.eq.3) then
        pfunq = pdf-ref
      elseif(i.eq.4) then
        pfunq = (pdf-ref)/ref
      else
        pfunq = 0.D0
      endif

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








