C     ------------------------------------------------------------------
      program batune00
C     ------------------------------------------------------------------
C--   Qcdnum tuning for BAT - 00
C--   - Evolve Micky Mouse pdfs from 1-100 GeV to get input for next job
C--   - VFNS-(2,25) NNLO on 100x50 grid with xmin = 1.D-3
C--   - Spline all pdfs (no top) and write them to disk (11 files)
C     ------------------------------------------------------------------
      implicit double precision (a-h,o-z)

*      include '../splint/inc/splint.inc'
*      include '../splint/inc/spliws.inc'
      
      real tim1, tim2
      
      data iord/3/, nfin/0/                                  !NNLO, VFNS
      data as0/0.364/, r20/2.D0/                           !input alphas
      external func                                          !input pdfs
      dimension def(-6:6,12)                   !input flavor composition
      data def  /
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,    ! 1 = dval
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,    ! 2 = dbar
     + 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,    ! 3 = uval
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,    ! 4 = ubar
     + 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,    ! 5 = sval
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,    ! 6 = sbar
     + 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,    ! 7 = cval
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,    ! 8 = cbar
     + 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,    ! 9 = bval
     + 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,    !10 = bbar
     + 26*0.    /
      dimension xmin(5), iwt(5)
      data xmin/1.D-3,0.2D0,0.4D0,0.6D0,0.75D0/                  !x grid
      data iwt / 1, 2, 4, 8, 16 /                                !x grid
      data ngx/5/, nxin/100/, iosp/3/                            !x grid
      dimension qq(4),wt(4)                                    !mu2 grid
      data qq/1.D0,2.D0,25.D0,1.D2/, wt/4*1.D0/                !mu2 grid
      data ngq/4/, nqin/50/                                    !mu2 grid

      external sfun, pfun
      dimension ia0(0:10), ia1(0:10)                   !spline addresses
      dimension pars(13)                            !evolution parmeters
      character*2 pnam(0:10)
      data pnam /'gl','dv','db','uv','ub','sv','sb','cv','cb','bv','bb'/

C--   lun = 6 stdout, -6 stdout w/out banner page
      lun    = 6
      lunout = abs(lun)
      call qcinit(lun,' ')                                  !init QCDNUM
      call ssp_SpInit(13)                                   !init SPLINT
      call gxmake(xmin,iwt,ngx,nxin,nx,iosp)                     !x-grid
      call gqmake(qq,wt,ngq,nqin,nq)                           !mu2-grid
      call wtfile(1,'../weights/bat00.wgt')           !calculate weights
      call zmfillw(nusedf)                       !structure function wts
      call setord(iord)                                   !LO, NLO, NNLO
      call setalf(as0,r20)                                 !input alphas
      iqc  = iqfrmq(qq(2))                              !charm threshold
      iqb  = iqfrmq(qq(3))                             !bottom threshold
      call setcbt(nfin,iqc,iqb,999)               !thesholds in the VFNS
      iq0  = iqfrmq(qq(1))                                  !start scale
      iq1  = nq                                               !end scale
      qq1  = qfrmiq(iq1)                                      !end scale

C--   Print settings
      call PrSettings(lunout,def)

C--   Evolution
      call cpu_time(tim1)
      call evolfg(1,func,def,iq0,eps)                            !evolve
      call cpu_time(tim2)
      write(lunout,'(/'' Evol  : '',F10.2,'' ms'')') 1000*(tim2-tim1)

C--   Spline pdfs at iq0 and iq1
      call cpu_time(tim1)
      do ipdf = 0,10
        ia0(ipdf)  = isp_SxMake(5)                         !spline at q0
        ia1(ipdf)  = isp_SxMake(5)                         !spline at q1
        if(ipdf.eq.0) then                                        !gluon
          do i = 1,13
            call ssp_Uwrite(i,0.D0)
          enddo
          call ssp_Uwrite(7,1.D0)
        else                                                     !quarks
          do i = 1,13
            call ssp_Uwrite(i,def(i-7,ipdf))
          enddo
        endif
        call ssp_SxFill(ia0(ipdf),sfun,iq0)
        call ssp_SxFill(ia1(ipdf),sfun,iq1)
      enddo
      call cpu_time(tim2)
      write(lunout,'( '' Splint: '',F10.2,'' ms'')') 1000*(tim2-tim1)

C--   Store evolution parameters, scale, and dump gluon spline at iq1
      call cpypar(pars,13,0)
      do i = 1,13
        call ssp_SpSetVal(ia1(0),i,pars(i))            !store parameters
      enddo
      call ssp_SpSetVal(ia1(0),14,qq1)                      !store scale
      call ssp_SpDump(ia1(0),'../splines/'//pnam(0)//'.spl')
C--   Store definition and dump pdf splines at iq1
      do i = 1,10
        do j = -6,6
          call ssp_SpSetVal(ia1(i),j+7,def(j,i))   !store pdf definition
        enddo
        call ssp_SpDump(ia1(i),'../splines/'//pnam(i)//'.spl')
      enddo

C--   Write spline addresses into the store
      do i = 0,10
        call ssp_Uwrite(i+1,dble(ia1(i)))
      enddo
C--   Get plot limits
      call ssp_SpLims(ia1(0),nu,umi,uma,nv,vmi,vma,nactive)
C--   Plot pdfs and stfs at iq1 (16 columns)
      call ffplot('../plots/batune00.dat',pfun,16,umi,uma,-200,' ')

      nw = isp_SpSize(0)
      write(lunout,'(/'' Memory size = '',I9)') nw
      nw = isp_SpSize(1)
      write(lunout,'( '' Words  used = '',I9)') nw
      nw = isp_SpSize(ia0(0))
      write(lunout,'( '' Spline size = '',I9)') nw
      nw = isp_SpSize(2)
      write(lunout,'( '' JunkAddress = '',I9)') nw

*      call ssp_Erase(1)
*      nw = isp_SpSize(1)
*      write(lunout,'(/'' Words  used = '',I9)') nw
*      call ssp_Mprint

      end

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
C     Spline input function --------------------------------------------
C     ------------------------------------------------------------------

C     ===========================================
      double precision function sfun(ix,iq,first)
C     ===========================================

C--   Function to spline

      implicit double precision(a-h,o-z)
      logical first

      save def
      dimension def(-6:6)

      if(first) then
        do i = -6,6
          def(i) = dsp_Uread(i+7)           !read def array
        enddo
      endif
      iset = 1
      if(def(0).eq.1.D0) then
        sfun = sumfij(iset,def,0,ix,iq,1)   !gluon
      else
        sfun = sumfij(iset,def,1,ix,iq,1)   !quarks
      endif

      return
      end

C     ------------------------------------------------------------------
C     Plot function  ---------------------------------------------------
C     ------------------------------------------------------------------

C     =========================================
      double precision function pfun(i,z,first)
C     =========================================

C--   Plot pdfs and stfs versus x at 100 GeV2
C--   1-11 pdf; 12=proton;   13=valence;
C--             14=FLproton; 15=F2proton; 16=xF3valence

      implicit double precision(a-h,o-z)
      logical first

      common/pass/ ia0(0:10), ia1(0:10)                !spline addresses

      dimension cpr(-6:6), pro(-6:6), val(-6:6)
      data cpr/ 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4./
      data val/-1.,-1.,-1.,-1.,-1.,-1., 0., 1., 1., 1., 1., 1., 1./

      if(first) continue

      do j = -6,6
          pro(j) = cpr(j)/9.D0
      enddo

      if(i.le.11) then
        ia = int(dsp_Uread(i))
        pfun = dsp_FunS1(ia,z,1)
      elseif(i.eq.12) then
        pfun = SUMFXQ( 1, pro, 1, z, 1.D2, 1 )
      elseif(i.eq.13) then
        pfun = BVALXQ( 1, 7, z, 1.D2, 1)
      elseif(i.eq.14) then
        call ZMSTFUN ( 1, pro, z, 1.D2, pfun, 1, 1 )
      elseif(i.eq.15) then
        call ZMSTFUN ( 2, pro, z, 1.D2, pfun, 1, 1 )
      elseif(i.eq.16) then
        call ZMSTFUN ( 3, val, z, 1.D2, pfun, 1, 1 )
      else
        pfun = 0.D0
      endif

      return
      end

C     ------------------------------------------------------------------
C     Evolution input functions ----------------------------------------
C     ------------------------------------------------------------------

C     ======================================
      double precision function func(ipdf,x)
C     ======================================

      implicit double precision (a-h,o-z)

                     func = 0.D0
      if(ipdf.eq. 0) func = xglu(x)
      if(ipdf.eq. 1) func = xdval(x)
      if(ipdf.eq. 2) func = xdbar(x)
      if(ipdf.eq. 3) func = xuval(x)
      if(ipdf.eq. 4) func = xubar(x)
      if(ipdf.eq. 5) func = xsval(x)
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
      double precision function xglu(x)
C     =================================

      implicit double precision (a-h,o-z)

      data ag    /1.7D0 /

      xglu = ag * x**(-0.1D0) * (1.D0-x)**5.D0

      return
      end
 
C     ==================================
      double precision function xdval(x)
C     ==================================

      implicit double precision (a-h,o-z)

      data ad /3.064320D0/

      xdval = ad * x**0.8D0 * (1.D0-x)**4.D0

      return
      end

C     ==================================
      double precision function xuval(x)
C     ==================================

      implicit double precision (a-h,o-z)

      data au /5.107200D0/

      xuval = au * x**0.8D0 * (1.D0-x)**3.D0

      return
      end

C     ==================================
      double precision function xsval(x)
C     ==================================

      implicit double precision (a-h,o-z)

      xsval = 0.D0*x

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







