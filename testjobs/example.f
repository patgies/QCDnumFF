C     ----------------------------------------------------------------
      program example
C     ----------------------------------------------------------------
C--   12-03-07  Basic QCDNUM example job
C     ----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      data as0/0.364/, r20/2.D0/, iord/3/, nfin/0/   !alphas, NNLO, VFNS
      data itype/3/, iset/1/                   !unpol/pol/tlike, pdf set
      data idbug/0/                                     !debug printpout
      data xmin/1.D-4/, nxin/100/, iosp/3/               !x grid, splord
      dimension qq(2),wt(2)                                    !mu2 grid
      data qq/2.D0,1.D4/, wt/1.D0,1.D0/, nqin/60/              !mu2 grid
      data q2c/3.D0/, q2b/25.D0/, q0/2.0/               !thresh and mu20
      data x/1.D-3/, q/1.D3/, qmz2/8315.25D0/             !output scales

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

      dimension pdf(-6:6)                                        !pdfout

C--   Weight files
      character*26 fnam(3)
C--               12345678901234567890123456
      data fnam /'../weights/unpolarised.wgt',
     +           '../weights/polarised.wgt  ',
     +           '../weights/timelike.wgt   '/

C--   Set-up -----------------------------------------------------------
      lun    = 6                                   ! -6 supresses banner
      call qcinit(lun,' ')                                   !initialize
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      call wtfile(itype,fnam(itype))                  !calculate weights
      call setord(iord)                                   !LO, NLO, NNLO
      call setalf(as0,r20)                                 !input alphas
      iqc  = iqfrmq(q2c)                                !charm threshold
      iqb  = iqfrmq(q2b)                               !bottom threshold
      call setcbt(nfin,iqc,iqb,999)               !thesholds in the vfns
C--   Evolution --------------------------------------------------------
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                                     !start scale
      call setint('edbg',idbug)                          !debug printout
      call evolfg(jset,func,def,iq0,eps)               !evolve all pdf's
C--   Get results   ----------------------------------------------------
      call allfxq(iset,x,q,pdf,0,1)               !interpolate all pdf's
      csea = 2.D0*pdf(-4)                             !charm sea at x,Q2
      asmz = asfunc(qmz2,nfout,ierr)                        !alphas(mz2)
C--   Print -----------------------------------------------------------
      call getint('lunq',lunout)
      write(lunout,'('' x, q, CharmSea ='',3E13.5)') x,q,csea
      write(lunout,'('' as(mz2)        ='', E13.5)') asmz

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







