C     ------------------------------------------------------------------
      program batune03
C     ------------------------------------------------------------------
C--   Qcdnum tuning for BAT - 03
C--   - Read input pdfs created by batune00
C--   - Evolve with 100x50 grid
C--   - Read batune02 reference splines F2, FL, xF3
C--   - Tune spline nodes of F2, FL, xF3
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
      character*2 Fn

      dimension xnd(100), qnd(100)

      dimension cpr(-6:6), pro(-6:6), val(-6:6)
      data cpr/ 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4./
      data val/-1.,-1.,-1.,-1.,-1.,-1., 0., 1., 1., 1., 1., 1., 1./

      external funsg

      do j = -6,6
        pro(j) = cpr(j)/9.D0
      enddo

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
      call wtfile(1,'../weights/bat03.wgt')           !calculate weights
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
      call evolfg(1,func,def,iq0,eps)                            !evolve

C--   Structure function weights
      call zmwfile('../weights/zmstf.wgt')

C--   Read reference splines
      irFL = isp_SpRead('../splines/FL300x150.spl')
      irF2 = isp_SpRead('../splines/F2300x150.spl')
      irF3 = isp_SpRead('../splines/F3300x150.spl')
C--   Make 2-dim structure function splines
      istx = 5
      istq = 10
      scut = 0.D0
C--   Just to get the nodes
      ia     = isp_S2Make(istx,istq)
      call ssp_Unodes(ia, xnd, 100, nu)
      call ssp_Vnodes(ia, qnd, 100, nv)
*      call ssp_Nprint(ia)
      call ssp_Erase(ia)
C--   Edit nodes
      xnd(90) = 0.13D0
      xnd(91) = 0.16D0
      xnd(92) = 0.33D0
C--   Make 2-dim structure function splines
      call cpu_time(tim1)
      iaFL = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaFL,1,pro,1,0.D0)
      iaF2 = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaF2,1,pro,2,0.D0)
      iaF3 = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaF3,1,val,3,scut)
      call cpu_time(tim2)
      write(lunout,
     +        '(/'' Stfsplines  : '',F10.2,'' ms'')') 1000*(tim2-tim1)

C--   Momentum sum rule
      iasg  = isp_S2Make(3,3)
      call ssp_S2Fill(iasg,funsg,0.D0)
      x1    = 1.D-3
      x2    = 1.D0
      q1    = 1.D2
      q2    = 3.D4
      sval  = dsp_IntS2(iasg,x1,x2,q1,q2,0.D0,5)/(q2-q1)
      write(lunout,'(/'' Sum rule : '',F10.3)') sval

C--   Structure function plots -----------------------------------------

C--   Decide here which structure function to plot
      Fn  = 'F3'
C--   Pass spline addresses
      call ssp_Uwrite(1,dble(irF3))
      call ssp_Uwrite(2,dble(iaF3))

C--   2-dim structure function plot
      call fsplot('../plots/bt03'//Fn//'xq.dat',4,0.9D0,0)
      call fnplot('../plots/bt03ndxq.dat',iaFL)
      call fcplot('../plots/bt03rsxq.dat',iaFL,300.D0)
*      call xmplot('../plots/xminvsq.dat',isptwo)
*      call qmplot('../plots/qmaxvsx.dat',isptwo)

C--   1-dim plot fixed values
      xp1 = 1.D-3
      xp2 = 1.D-2
      xp3 = 1.D-1
      qp1 = 1.D2
      qp2 = 2.D3
      qp3 = 3.D4

C--   Plot proton F123 vs x
      xmi  = xmin(1)
      xma  = 0.9D0
      call stfvsx('../plots/bt03'//Fn//'x1.dat',xmi,xma,qp1)
      call stfvsx('../plots/bt03'//Fn//'x2.dat',xmi,xma,qp2)
      call stfvsx('../plots/bt03'//Fn//'x3.dat',xmi,xma,qp3)

C--   Plot proton F123 vs q2
      qmi  = qq(1)
      qma  = qq(3)
      call stfvsq('../plots/bt03'//Fn//'q1.dat',qmi,qma,xp1)
      call stfvsq('../plots/bt03'//Fn//'q2.dat',qmi,qma,xp2)
      call stfvsq('../plots/bt03'//Fn//'q3.dat',qmi,qma,xp3)

      nw = isp_SpSize(0)
      write(lunout,'(/'' Memory size = '',I9)') nw
      nw = isp_SpSize(1)
      write(lunout,'( '' Words  used = '',I9)') nw
      nw = isp_SpSize(iaF2)
      write(lunout,'( '' Spline size = '',I9)') nw

*      call ssp_Mprint

      end

C     ------------------------------------------------------------------
C     Singlet + gluon filling routine ----------------------------------
C     ------------------------------------------------------------------

C     ============================================
      double precision function funsg(ix,iq,first)
C     ============================================

C--   Function to spline gluon + singlet

      implicit double precision(a-h,o-z)
      logical first

      if(first) continue

      funsg = BVALIJ(1,0,ix,iq,1) + BVALIJ(1,1,ix,iq,1)

      return
      end

C     ------------------------------------------------------------------
C     Surface plot routines --------------------------------------------
C     ------------------------------------------------------------------

C     =====================================
      subroutine fsplot(fnam,ipl,xmax,ichk)
C     =====================================

C--   Write file for 2-dim surface plot
C--   ipl: 1=ref, 2=stf, 3=stf-ref, 4=(stf-ref)/ref

      implicit double precision(a-h,o-z)

      character*(*) fnam

      dimension xarr(200), qarr(200)

C--   Open file
      lun = imb_NextL(0)
      open(unit=lun,file=fnam,
     +              form='formatted',status='unknown',err=500)

C--   Spline addresses
      iaR = int(dsp_Uread(1))
      iaF = int(dsp_Uread(2))

C--   Get reference nodes
      call ssp_Unodes(iaR,xarr,200,nx)
      call ssp_Vnodes(iaR,qarr,200,nq)

C--   Write formatted file
      do ix = 1,nx
        xx = xarr(ix)
        if(xx.le.xmax) then
          do iq = 1,nq
            qq = qarr(iq)
            ff = pfun2(ipl,xx,qq,ichk)
            write(lun,'(3E13.5)') xx,qq,ff
          enddo
          write(lun,'('' '')')
        endif
      enddo

      write(6,'(/'' FSPLOT: write file '',A)') fnam

      close(lun)

      return

 500  continue
      stop ' FSPLOT: cannot open file'

      end


C     ===========================================
      double precision function pfun2(i,x,q,ichk)
C     ===========================================

C--   plot stf vs x,q

      implicit double precision(a-h,o-z)

      logical lmb_ne

      iaR = int(dsp_Uread(1))
      iaF = int(dsp_Uread(2))

      stref = dsp_FunS2(iaR, x, q, ichk)
      stfun = dsp_FunS2(iaF, x, q, ichk)

      if(i.eq.1) then
        pfun2 = stref
      elseif(i.eq.2) then
        pfun2 = stfun
      elseif(i.eq.3) then
        pfun2 = abs(stfun - stref)
      elseif(i.eq.4 .and. lmb_ne(stref,0.D0,1.D-12)) then
        pfun2 = abs((stfun - stref) / stref)
      else
        pfun2 = 0.D0
      endif

      return
      end


C     ============================
      subroutine fnplot(fnam, ias)
C     ============================

C--   Write file to create node-grid on surface plot

      implicit double precision(a-h,o-z)

      character*(*) fnam

      dimension xarr(100), qarr(100)

C--   Open file
      lun = imb_NextL(0)
      open(unit=lun,file=fnam,
     +              form='formatted',status='unknown',err=500)

C--   Get ispl nodes
      call ssp_Unodes(ias,xarr,100,nx)
      call ssp_Vnodes(ias,qarr,100,nq)

C--   Write formatted file
      do ix = 1,nx
        do iq = 1,nq
          xx = xarr(ix)
          qq = qarr(iq)
          ff = 0.D0
          write(lun,'(3E13.5)') xx,qq,ff
        enddo
        write(lun,'('' '')')
      enddo

      write(6,'(/'' FNPLOT: write file '',A)') fnam

      close(lun)

      return

 500  continue
      stop ' FNPLOT: cannot open file'

      end

C     ================================
      subroutine fcplot(fnam, ias, rs)
C     ================================

C--   Write file to plot rscut on surface plot

      implicit double precision(a-h,o-z)
      logical lmb_le, lmb_ge

      character*(*) fnam

      dimension xarr(100), qarr(100)

C--   Open file
      lun = imb_NextL(0)
      open(unit=lun,file=fnam,
     +              form='formatted',status='unknown',err=500)

C--   Get ispl nodes
      call ssp_Unodes(ias,xarr,100,nx)
      call ssp_Vnodes(ias,qarr,100,nq)

      if(lmb_ge(rs*rs,qarr(nq)/xarr(1),-1.D-9)) goto 510
      if(lmb_le(rs*rs,qarr(1)/xarr(nx),-1.D-9)) goto 510

C--   Write formatted file
      np = 20
C--   Bracket y
      rsc = log(rs*rs)
      ymi = -log(xarr(nx))
      yma = -log(xarr(1))
      tmi =  log(qarr(1))
      tma =  log(qarr(nq))
      y1  =  max(rsc-tma,ymi)
      y2  =  min(rsc-tmi,yma)
      bw  =  (y2-y1)/(np-1)
C--   Loop over np plot points
      do i = 1,np
        yi = y1 + (i-1)*bw
        ti = rsc-yi
        xx = exp(-yi)
        qq = exp( ti)
        ff = 0.D0
        write(lun,'(3E13.5)') xx,qq,ff
      enddo

      write(6,'(/'' FCPLOT: write file '',A)') fnam

      close(lun)

      return

 500  continue
      stop ' FCPLOT: cannot open file'
 510  continue
      write(lun,'(3E13.5)') xarr(1),qarr(1),0.D0

      end

C     ============================
      subroutine xmplot(fnam, ias)
C     ============================

C--   Write file to plot xmin on surface plot

      implicit double precision(a-h,o-z)

      include '../splint/inc/splint.inc'
      include '../splint/inc/spliws.inc'

      character*(*) fnam

      dimension iatab(4), iatmp(4)

C--   Open file
      lun = imb_NextL(0)

      open(unit=lun,file=fnam,
     +              form='formatted',status='unknown',err=500)

      call sspGetIaTwoD(w,ias,iat,iau,nu,iav,nv,iatab,iatmp)

      do it = 1,nv
        tt  = w(iav+it-1)
        iym = max(int(w(iav+it+nv-1)),1)
        ym  = w(iau+iym-1)
        xx  = exp(-ym)
        qq  = exp( tt)
        ff  = 0.D0
        write(lun,'(3E13.5)') xx,qq,ff
      enddo

      write(6,'(/'' XMPLOT: write file '',A)') fnam

      close(lun)

      return

 500  continue
      stop ' XMPLOT: cannot open file'
* 510  continue
*      stop ' XMPLOT: error writing file'

      end

C     ============================
      subroutine qmplot(fnam, ias)
C     ============================

C--   Write file to plot qmax on surface plot

      implicit double precision(a-h,o-z)

      include '../splint/inc/splint.inc'
      include '../splint/inc/spliws.inc'

      character*(*) fnam

      dimension iatab(4), iatmp(4)

C--   Open file
      lun = imb_NextL(0)

      open(unit=lun,file=fnam,
     +              form='formatted',status='unknown',err=500)

      call sspGetIaTwoD(w,ias,iat,iau,nu,iav,nv,iatab,iatmp)

      do iy = 1,nu
        yy  = w(iau+iy-1)
        itm = max(int(w(iau+iy+nu-1)),1)
        tm  = w(iav+itm-1)
        xx  = exp(-yy)
        qq  = exp( tm)
        ff  = 0.D0
        write(lun,'(3E13.5)') xx,qq,ff
      enddo

      write(6,'(/'' QMPLOT: write file '',A)') fnam

      close(lun)

      return

 500  continue
      stop ' QMPLOT: cannot open file'

      end

C     ------------------------------------------------------------------
C     1-dim plot routines ----------------------------------------------
C     ------------------------------------------------------------------

C     ==================================
      subroutine stfvsx(fnam,xmi,xma,q2)
C     ==================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfunx

      call ssp_Uwrite(3,q2)
      call ssp_Uwrite(4,xmi)
      call ffplot(fnam,pfunx,5,xmi,xma,-300,' ')

      return
      end


C     ==========================================
      double precision function pfunx(i,z,first)
C     ==========================================

C--   plot stf vs x

      implicit double precision(a-h,o-z)
      logical first

      save iaR, iaF, qfix, xcut

      if(first) then
        iaR  = int(dsp_Uread(1))
        iaF  = int(dsp_Uread(2))
        qfix = dsp_Uread(3)
        xcut = dsp_Uread(4)
      endif


      stref = dsp_FunS2(iaR, z, qfix, 1)
      stfun = dsp_FunS2(iaF, z, qfix, 1)

      if(i.eq.1) then
        pfunx = stref
      elseif(i.eq.2 .and. z.ge.xcut) then
        pfunx = stref
      elseif(i.eq.3 .and. z.ge.xcut) then
        pfunx = stfun
      elseif(i.eq.4 .and. z.ge.xcut) then
        pfunx = stfun - stref
      elseif(i.eq.5 .and. z.ge.xcut) then
        pfunx = (stfun - stref) / stref
      else
        pfunx = 0.D0
      endif

      return
      end

C     ==================================
      subroutine stfvsq(fnam,qmi,qma,xx)
C     ==================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfunq

      call ssp_Uwrite(3,xx)
      call ssp_Uwrite(4,qma)
      call ffplot(fnam,pfunq,5,qmi,qma,-300,' ')

      return
      end

C     ==========================================
      double precision function pfunq(i,z,first)
C     ==========================================

C--   plot stf vs q

      implicit double precision(a-h,o-z)
      logical first

      save iaR, iaF, xfix, qcut

      if(first) then
        iaR  = int(dsp_Uread(1))
        iaF  = int(dsp_Uread(2))
        xfix = dsp_Uread(3)
        qcut = dsp_Uread(4)
      endif

      if(first) continue

      stref = dsp_FunS2(iaR, xfix, z, 1)
      stfun = dsp_FunS2(iaF, xfix, z, 1)

      if(i.eq.1) then
        pfunq = stref
      elseif(i.eq.2 .and. z.le.qcut) then
        pfunq = stref
      elseif(i.eq.3 .and. z.le.qcut) then
        pfunq = stfun
      elseif(i.eq.4 .and. z.le.qcut) then
        pfunq = stfun - stref
      elseif(i.eq.5 .and. z.le.qcut) then
        pfunq = (stfun - stref) / stref
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








