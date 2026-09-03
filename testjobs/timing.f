C     ----------------------------------------------------------------
      program timing
C     ----------------------------------------------------------------
C--   Enjoy the speed of QCDNUM and ZMSTF
C     ----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      
      real rmb_urand, rval, tim1, tim2
      
      data iord/3/, nfin/0/                                  !NNLO, VFNS
      data as0/0.364/, r20/2.D0/                           !input alphas
      external func                                          !input pdfs
      dimension def(-6:6,12)                   !input flavor composition
      data def  /
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,     !1 = dval
     + 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,     !2 = uval
     + 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,     !3 = sval
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,     !4 = dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,     !5 = ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,     !6 = sbar
     + 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,     !7 = cval
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,     !8 = cbar
     + 52*0.    /
      dimension xmin(5), iwt(5)
      data xmin/1.D-5,0.2D0,0.4D0,0.6D0,0.75D0/                  !x grid
      data iwt / 1, 2, 4, 8, 16 /                                !x grid
      data ngx/5/, nxin/100/, iosp/3/                            !x grid
      dimension qq(2),wt(2)                                    !mu2 grid
      data qq/2.D0,1.D4/, wt/1.D0,1.D0/                        !mu2 grid
      data ngq/2/, nqin/50/                                    !mu2 grid
      data q2c/3.D0/, q2b/25.D0/, q0/2.0/               !thresh and mu20
      dimension xx(1000),q2(1000),ff(1000)                       !output
      dimension proton(-6:6)                                     !proton
      data proton/4.,1.,4.,1.,4.,1.,0.,1.,4.,1.,4.,1.,4./
      
C--   Divide proton by 9
      do i = -6,6
        proton(i) = proton(i)/9.D0
      enddo

C--   Generate 1000 random x-mu2 points in the HERA kinematic range
      iseed = 0
      roots = 300.D0                                     !HERA COM enery
      shera = roots*roots
      xlog1 = log(xmin(1))
      xlog2 = log(0.99D0)
      qlog1 = log(qq(1))
      qlog2 = log(qq(ngq))
      ntot  = 0
      do while(ntot.lt.1000)
        rval  = rmb_urand(iseed)
        xlog  = xlog1 + rval*(xlog2-xlog1)
        xxxx  = exp(xlog)
        rval  = rmb_urand(iseed)
        qlog  = qlog1 + rval*(qlog2-qlog1)
        qqqq  = exp(qlog)
        if(qqqq.le.xxxx*shera) then
          ntot     = ntot+1
          xx(ntot) = xxxx
          q2(ntot) = qqqq
        endif
      enddo

C--   lun = 6 stdout, -6 stdout w/out banner page
      lun    = 6
      lunout = abs(lun)      
      call qcinit(lun,' ')                                   !initialize
      call gxmake(xmin,iwt,ngx,nxin,nx,iosp)                     !x-grid
      call gqmake(qq,wt,ngq,nqin,nq)                           !mu2-grid
      call wtfile(1,'../weights/unpolarised.wgt')     !calculate weights
      call zmfillw(nwords)                                !weights ZMSTF
      call setord(iord)                                   !LO, NLO, NNLO
      call setalf(as0,r20)                                 !input alphas
      iqc  = iqfrmq(q2c)                                !charm threshold
      iqb  = iqfrmq(q2b)                               !bottom threshold
      call setcbt(nfin,iqc,iqb,999)               !thesholds in the VFNS
      iq0  = iqfrmq(q0)                                     !start scale

C--   Evolution and structure function timing loop
      write(lunout,
     +     '(/'' Wait: 1000 evols and 2.10^6 stfs will take ...'')')
      call cpu_time(tim1)
      do iter = 1,1000
        call evolfg(1,func,def,iq0,eps)                          !evolve
        call zmstfun(1,proton,xx,q2,ff,1000,1)                !FL proton
        call zmstfun(2,proton,xx,q2,ff,1000,1)                !F2 proton
      enddo  
      call cpu_time(tim2)
      write(lunout,'(/'' CPU : '',F10.4,'' secs'')') tim2-tim1
      end
      
C     ----------------------------------------------------------------

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







