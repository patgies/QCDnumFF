C     ------------------------------------------------------------------
      program getgluon
C     ------------------------------------------------------------------
C--   Use different interpolation routines to get the gluon
C     ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      
      data as0/0.364/, r20/2.D0/, iord/3/, nfin/0/   !alphas, NNLO, VFNS
      dimension xmin(5), iwt(5)                                  !x grid
      data xmin/1.D-5,0.2D0,0.4D0,0.6D0,0.75D0/                  !x grid
      data iwt / 1, 2, 4, 8, 16 /                                !x grid
      data nxin/100/, iosp/3/                                    !x grid
      dimension qlim(2),wt(2)                                  !mu2 grid
      data qlim/2.D0,1.D4/, wt/1.D0,1.D0/, nqin/60/            !mu2 grid
      data q2c/3.D0/, q2b/25.D0/, q0/2.0/                    !thresh, q0

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

      external func                                          !input pdfs

      data xx /1.D-4/, qq /1.D2/                                 !output
      dimension pdf(-6:6), coef(-6:6)                            !output
      data coef /13*0.D0/

      call QCINIT(6,' ')                                           !init
      call GXMAKE(xmin,iwt,5,nxin,nx,iosp)                        !xgrid
      call GQMAKE(qlim,wt,2,nqin,nq)                              !qgrid
      call WTFILE(1,'../weights/unpolarised.wgt')              !wt unpol
      call SETALF(as0,r20)                                     !input as
      iqc  = IQFRMQ(q2c)                                          !charm
      iqb  = IQFRMQ(q2b)                                         !bottom
      call SETCBT(nfin,iqc,iqb,999)                         !VFNS thresh
      iq0  = IQFRMQ(q0)                                     !input scale
      call SETORD(iord)                                     !LO,NLO,NNLO
      call EVOLFG(1,func,def,iq0,eps)                            !evolve

      call GETINT('lunq',lunout)               !get logicasl unit number

C--   All gluons below must be the same
      gluon = BVALXQ(1,0,xx,qq,1)
      write(lunout,'( '' BVALXQ gluon = '',E15.5)') gluon

      gluon = FVALXQ(1,0,xx,qq,1)
      write(lunout,'( '' FVALXQ gluon = '',E15.5)') gluon

      call ALLFXQ(1,xx,qq,pdf,0,1)
      write(lunout,'( '' ALLFXQ gluon = '',E15.5)') pdf(0)

      gluon = sumfxq(1,coef,0,xx,qq,1)
      write(lunout,'( '' SUMFXQ gluon = '',E15.5)') gluon

      call FFLIST(1,coef,0,xx,qq,gluon,1,1)
      write(lunout,'( '' FFLIST gluon = '',E15.5)') gluon

      call FTABLE(1,coef,0,xx,1,qq,1,gluon,1)
      write(lunout,'( '' FTABLE gluon = '',E15.5)') gluon

      ix = ixfrmx(xx)
      iq = iqfrmq(qq)
      gluon = BVALIJ(1,0,ix,iq,1)
      write(lunout,'(/'' BVALIJ gluon = '',E15.5)') gluon

      gluon = FVALIJ(1,0,ix,iq,1)
      write(lunout,'( '' FVALIJ gluon = '',E15.5)') gluon

      call ALLFIJ(1,ix,iq,pdf,0,1)
      write(lunout,'( '' ALLFIJ gluon = '',E15.5)') pdf(0)

      gluon = SUMFIJ(1,coef,0,ix,iq,1)
      write(lunout,'( '' SUMFIJ gluon = '',E15.5)') gluon

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







