C     ------------------------------------------------------------------
      program batune04
C     ------------------------------------------------------------------
C--   Qcdnum tuning for BAT - 04
C--   - Read input pdfs created by batune00
C--   - Evolve with 100x50 grid
C--   - Read batune02 reference splines F2, FL, xF3
C--   - Make 22x7 splines for F2, FL, xF3 up and down
C--   - Tune cross section spline
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

      external xfun

      dimension iaF(7), xnd(100), qnd(100)

      dimension dplus(-6:6)
      data dplus / 0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0,
     +             1.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0/
      dimension dminu(-6:6)
      data dminu / 0.D0,-1.D0, 0.D0,-1.D0, 0.D0,-1.D0, 0.D0,
     +             1.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0/
      dimension uplus(-6:6)
      data uplus / 1.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     +             0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 1.D0/
      dimension uminu(-6:6)
      data uminu /-1.D0, 0.D0,-1.D0, 0.D0,-1.D0, 0.D0, 0.D0,
     +             0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 1.D0/

      dimension cpr(-6:6), pro(-6:6), val(-6:6)
      data cpr/ 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4./
      data val/-1.,-1.,-1.,-1.,-1.,-1., 0., 1., 1., 1., 1., 1., 1./

      dimension bx1(500), bx2(500), bq1(500), bq2(500)
      dimension xint(500), rint(500)

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
      iaF(1) = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaF(1),1,dplus,1,scut)
      iaF(2) = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaF(2),1,uplus,1,scut)
      iaF(3) = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaF(3),1,dplus,2,scut)
      iaF(4) = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaF(4),1,uplus,2,scut)
      iaF(5) = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaF(5),1,dminu,3,scut)
      iaF(6) = isp_S2User(xnd,100,qnd,100)
      call ssp_S2F123(iaF(6),1,uminu,3,scut)
      call cpu_time(tim2)
      write(lunout,
     +        '(/'' Stfsplines  : '',F10.2,'' ms'')') 1000*(tim2-tim1)

C--   Store spline addresses
      do i = 1,6
        call ssp_Uwrite(i,dble(iaF(i)))
      enddo
      call ssp_Uwrite(8,dble(irFL))
      call ssp_Uwrite(9,dble(irF2))

C--   cross-section spline
      istx   = 1
      istq   = 2
      rscut  = 370.D0
      rs     = 300.D0
      call cpu_time(tim1)
      iaF(7) = isp_S2Make(istx,istq)
      call ssp_Uwrite(7,dble(iaF(7)))
      call ssp_SpSetVal(iaF(7),1,rs)
      call ssp_S2Fill(iaF(7),xfun,rscut)
      call cpu_time(tim2)
      write(lunout,
     +        '( '' Xsecspline  : '',F10.2,'' ms'')') 1000*(tim2-tim1)

      rsc = dsp_RsCut(iaF(7))
      rsm = dsp_RsMax(iaF(7),rsc)

      write(6,'(/'' RsCut =     '',F10.0)') rsc
      write(6,'( '' RsMax =     '',F10.0)') rsm

C--   Read integration bins
      call ReadBins('../splines/BatBins.dat',bx1,bx2,bq1,bq2,500,nbins)

C--   Integrate
      call cpu_time(tim1)
      roots = 300.D0
      np    = 5
      ntry  = 100
      ncros = 0
      do j = 1,ntry
      do i = 1,nbins
        icros   = iCrossSc(bx1(i),bx2(i),bq1(i),bq2(i),roots)
        xint(i) = dsp_IntS2(iaf(7),bx1(i),bx2(i),bq1(i),bq2(i),roots,np)
        if(icros.eq.1) ncros = ncros+1
      enddo
      enddo
      call cpu_time(tim2)
      write(lunout,'('' Bins cross  : '',I10)') ncros/ntry
      write(lunout,
     +  '( '' Integration : '',F10.2,'' ms'')') 1000*(tim2-tim1)/ntry

C--   See if integration OK
*      call DumpBint('../splines/IntBins2.dat',xint,nbins)
      call ReadBint('../splines/IntBins2.dat',rint,nbins)
      call CompBint(xint,rint,nbins)

C--   Cross-section plots ----------------------------------------------

C--   2-dim cross-section plot
C--   ipl: 1=ref, 2=stf, 3=stf-ref, 4=(stf-ref)/ref
      ipl    = 4
      call xsplot('../plots/bt04Xsexq.dat',ipl,0.9D0,0)
      call fnplot('../plots/bt04Xndxq.dat',iaF(7))
      call fcplot('../plots/bt04Xrsxq.dat',iaF(7),rs)
      call fcplot('../plots/bt04Xyyxq.dat',iaF(7),0.7D0*rs)
*      call xmplot('../plots/xminvsq.dat',isptwo)
*      call qmplot('../plots/qmaxvsx.dat',isptwo)

C--   Plot xsec vs x
      xmi  = xmin(1)
      xma  = 1.0D0
      call pltvsx('../plots/bt04xvsx1.dat',xmi,xma,qq(1))
      call pltvsx('../plots/bt04xvsx2.dat',1.D-2,xma, 2.D3)
      call pltvsx('../plots/bt04xvsx3.dat',1.D-1,xma,qq(3))

C--   Plot xsec vs q2
      qmi  = qq(1)
      qma  = qq(3)
      call pltvsq('../plots/bt04xvsq1.dat',qmi,1.D3,5.D-3)
      call pltvsq('../plots/bt04xvsq2.dat',qmi,1.D4,5.D-2)
      call pltvsq('../plots/bt04xvsq3.dat',qmi,qma,5.D-1)

C--   Plot xsec along y=1
      call pltvsy('../plots/bt04xaty1.dat',xmi,xma,1.0D0)
      call pltvsy('../plots/bt04xaty2.dat',xmi,xma,0.7D0)

      nw = isp_SpSize(0)
      write(lunout,'(/'' Memory size = '',I9)') nw
      nw = isp_SpSize(1)
      write(lunout,'( '' Words  used = '',I9)') nw
      nw = isp_SpSize(ia0(0))
      write(lunout,'( '' Spline size = '',I9)') nw

*      call ssp_Mprint

      end

C     ------------------------------------------------------------------
C     Integration bins -------------------------------------------------
C     ------------------------------------------------------------------

C     =============================================
      integer function iCrossSc(x1, x2, q1, q2, rs)
C     =============================================

C--   0 = below cut, 1 = crossed by cut, 2 above cut
C--
C--   y1,2    (in): bin limits in y
C--   t1,2    (in): bin limits in t
C--   sc      (in): log roots^2; 0=no cut
C--
C--   Returns 0 if no cut

      implicit double precision(a-h,o-z)
      logical lmb_ge, lmb_le

      deps0 = 1.D-9
      y1 = -log(x2)
      y2 = -log(x1)
      t1 =  log(q1)
      t2 =  log(q2)

      if(lmb_le(rs,0.D0,-deps0)) then
        iCrossSc = 0
      else
        sc = log(rs*rs)
        if(lmb_le(t2+y2,sc,-deps0)) then
          iCrossSc = 0
        elseif(lmb_ge(t1+y1,sc,-deps0)) then
          iCrossSc = 2
        else
          iCrossSc = 1
        endif
      endif
        
      return
      end

C     ===============================================
      subroutine ReadBins(fname,bx1,bx2,bq1,bq2,n,nb)
C     ===============================================

C--   Read cross-section integration bins

      implicit double precision (a-h,o-z)

      character*(*) fname
      dimension bx1(*), bx2(*), bq1(*), bq2(*)
      character*10 txt

      lun = imb_NextL(0)
      open(unit=lun,file=fname,
     +          form='formatted',status='unknown',err=500)
      idum = 0
      nb   = 0
      do while(idum.eq.0)
        nb = nb+1
        if(nb.gt.n) stop ' READBINS: bin array size exceeded'
        read(unit=lun,fmt=*,end=100,err=510)
     +                 bx1(nb),bx2(nb),bq1(nb),bq2(nb)
      enddo

 100  continue
      nb = nb-1
      call smb_itoch(nb,txt,m)
      write(6,'(/'' READBINS: read '',A,'' bins from '',A)')
     +                 txt(1:m),fname

      close(lun)
      return

 500  continue
      stop ' READBINS: cannot open file'
 510  continue
      stop ' READBINS: error reading file'

      end

C     =====================================
      subroutine DumpBint(fname,xint,nbins)
C     =====================================

C--   Dump bin integrals

      implicit double precision (a-h,o-z)

      character*(*) fname
      dimension xint(*)

      lun = imb_NextL(0)
      open(unit=lun,file=fname,
     +          form='unformatted',status='unknown',err=500)

      write(lun,err=510) nbins
      do i = 1,nbins
        write(lun,err=510) xint(i)
      enddo

      close(lun)
      return

 500  continue
      stop ' DUMPBINT: cannot open file'
 510  continue
      stop ' DUMPBINT: error writing file'

      end

C     =====================================
      subroutine ReadBint(fname,rint,nbins)
C     =====================================

C--   Read bin integrals

      implicit double precision (a-h,o-z)

      character*(*) fname
      dimension rint(*)

      lun = imb_NextL(0)
      open(unit=lun,file=fname,
     +          form='unformatted',status='unknown',err=500)

      read(lun,err=510,end=510) nbins
      do i = 1,nbins
        read(lun,err=510,end=510) rint(i)
      enddo

      close(lun)
      return

 500  continue
      stop ' READBINT: cannot open file'
 510  continue
      stop ' READBINT: error reading file'

      end

C     =====================================
      subroutine CompBint(xint,rint,nbins)
C     =====================================

C--   Compare bin integrals

      implicit double precision (a-h,o-z)

      dimension xint(*), rint(*), dint(500)

      if(nbins.gt.500) stop 'CompBint: too many bins'

      call smb_vminv(xint,rint,dint,nbins)
      dev  = dmb_vnorm(1,dint,nbins)
      ref  = dmb_vnorm(1,rint,nbins)
      
      write(6,'(/'' Bin integrals dI/I = '',E13.5)') dev/ref

      return
      end

C     ------------------------------------------------------------------
C     cross-section function
C     ------------------------------------------------------------------

C     ===========================================
      double precision function xfun2(i,x,q,ichk)
C     ===========================================

C--   xseq vs x,q from splines
C--   i: 1=xref; 2=xref below limit; 3=xfun-xref; 4=(xfun-xref)/xref

      implicit double precision(a-h,o-z)

      logical lmb_ne

      iaX  = int(dsp_Uread(7))
      irFL = int(dsp_Uread(8))
      irF2 = int(dsp_Uread(9))

      FL   = dsp_FunS2(irFL, x, q, ichk)
      F2   = dsp_FunS2(irF2, x, q, ichk)
      F3   = 0.D0

      rs   = dsp_SpGetVal(iaX,1)
      s    = rs*rs
      qma  = x*s
      xref = xsec(x, q, s, FL, F2, F3)
      xfun = dsp_FunS2(iaX, x, q, ichk)

      if(i.eq.1) then
        xfun2 = xref
      elseif(i.eq.2 .and. q.le.qma) then
        xfun2 = xref
      elseif(i.eq.3 .and. lmb_ne(xfun,0.D0,-1.D-9)) then
        xfun2 = abs(xfun - xref)
      elseif(i.eq.4 .and. lmb_ne(xfun,0.D0,-1.D-9)) then
        xfun2 = abs((xfun - xref) / xref)
      else
        xfun2 = 0.D0
      endif

      return
      end


C     ===================================================
      double precision function xsec(x, q, s, FL, F2, F3)
C     ===================================================

C--   Compute xsec

      implicit double precision (a-h,o-z)

      data pi/3.141592654/
  
      alf  = 1.D0/137.D0
      fac  = 2.D0*pi*alf*alf
      y    = q/(x*s)
      omy  = 1.D0-y
      ypl  = 1.D0 + omy*omy
      ymi  = 1.D0 - omy*omy
      xsec = fac*(ypl*F2 - y*y*FL + ymi*F3)/(x*q*q)

      return
      end

C     ------------------------------------------------------------------
C     xsec spline input function ---------------------------------------
C     ------------------------------------------------------------------

C     =============================================
      double precision function xfun(ix, iq, first)
C     =============================================

C--   Function to spline xsec

      implicit double precision (a-h,o-z)
      logical first

      save iaF
      dimension iaF(7)

      if(first) then
        do i = 1,7
          iaF(i) = int(dsp_Uread(i))
        enddo
      endif

      x    = xfrmix(ix)
      q    = qfrmiq(iq)
      rs   = dsp_SpGetVal(iaF(7),1)
      s    = rs*rs
      FL   = (1.D0/9.D0)*dsp_FunS2(iaF(1), x, q, 1) +
     +       (4.D0/9.D0)*dsp_FunS2(iaF(2), x, q, 1)
      F2   = (1.D0/9.D0)*dsp_FunS2(iaF(3), x, q, 1) +
     +       (4.D0/9.D0)*dsp_FunS2(iaF(4), x, q, 1)
      F3   = 0.D0

      xfun = xsec(x, q, s, FL, F2, F3)

      return
      end

C     ------------------------------------------------------------------
C     Surface plot routines --------------------------------------------
C     ------------------------------------------------------------------

C     =====================================
      subroutine xsplot(fnam,ipl,xmax,ichk)
C     =====================================

C--   Write file for 2-dim cross-section plot
C--   ipl: 1=xref, 2=xref berlow limit , 3=xsec-xref, 4=(xsec-xref)/xref

      implicit double precision(a-h,o-z)

      character*(*) fnam

      dimension xarr(200), qarr(200)

C--   Open file
      lun = imb_NextL(0)
      open(unit=lun,file=fnam,
     +              form='formatted',status='unknown',err=500)

C--   Refernce spline address
      iaR = int(dsp_Uread(8))

C--   Get reference nodes
      call ssp_Unodes(iaR,xarr,200,nx)
      call ssp_Vnodes(iaR,qarr,200,nq)

C--   Write formatted file
      do ix = 1,nx
        xx = xarr(ix)
        if(xx.le.xmax) then
          do iq = 1,nq
            qq = qarr(iq)
            ff = xfun2(ipl,xx,qq,ichk)
            write(lun,'(3E13.5)') xx,qq,ff
          enddo
          write(lun,'('' '')')
        endif
      enddo

      write(6,'(/'' XSPLOT: write file '',A)') fnam

      close(lun)

      return

 500  continue
      stop ' XSPLOT: cannot open file'

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
      subroutine pltvsx(fnam,xmi,xma,q2)
C     ==================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfunx

      call ssp_Uwrite(10,q2)

      call ffplot(fnam,pfunx,4,xmi,xma,-300,' ')

      return
      end


C     ==========================================
      double precision function pfunx(i,z,first)
C     ==========================================

C--   plot xsec vs x

      implicit double precision(a-h,o-z)
      logical first

      save qfix

      if(first) then
        qfix = dsp_Uread(10)
      endif

      pfunx = xfun2(i,z,qfix,1)

      return
      end

C     ==================================
      subroutine pltvsq(fnam,qmi,qma,xx)
C     ==================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfunq

      call ssp_Uwrite(10,xx)

      call ffplot(fnam,pfunq,4,qmi,qma,-300,' ')

      return
      end

C     ==========================================
      double precision function pfunq(i,z,first)
C     ==========================================

C--   plot xsec vs q

      implicit double precision(a-h,o-z)
      logical first

      save xfix

      if(first) then
        xfix = dsp_Uread(10)
      endif

      pfunq = xfun2(i,xfix,z,1)

      return
      end

C     ====================================
      subroutine pltvsy(fnam,xmi,xma,yval)
C     ====================================

      implicit double precision (a-h,o-z)

      character*(*) fnam

      external pfuny

      iaX = int(dsp_Uread(7))
      rs  = dsp_SpGetVal(iaX,1)
      call ssp_Uwrite(10,rs*rs)
      call ssp_Uwrite(11,yval)

      call ssp_SpLims(iaX,nu,u1,u2,nv,v1,v2,n)

      xmin = max(xmi,v1/(yval*rs*rs))
      xmax = min(xma,v2/(yval*rs*rs))
      call ffplot(fnam,pfuny,4,xmin,xmax,-300,' ')

      return
      end

C     ==========================================
      double precision function pfuny(i,z,first)
C     ==========================================

C--   plot xsec along y = 1

      implicit double precision(a-h,o-z)
      logical first

      save s, y

      if(first) then
        s = dsp_Uread(10)
        y = dsp_Uread(11)
      endif

      pfuny = xfun2(i,z,z*y*s,1)

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
