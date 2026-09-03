C     ----------------------------------------------------------------
      program kkD0scan
C     ----------------------------------------------------------------
C--   Same physics as kkD0.f, but scans D_c^D0(x,Q2) over x at several
C--   Q values and dumps the result to a data file for plotting.
C     ----------------------------------------------------------------
      implicit double precision (a-h,o-z)

      external funcC, funcB
      dimension def(-6:6,12)

C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
      data def  /
     + 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,        !dval
     + 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,        !uval
     + 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,        !sval
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,        !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,        !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,        !sbar
     + 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,        !cval
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,        !cbar
     + 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,        !bval
     + 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,        !bbar
     + 26*0.    /                                            !tval,tbar

      dimension pdfC(-6:6), pdfB(-6:6)

      data qmc2  /2.25D0/, qmb2 /25.D0/                     !mc=1.5,mb=5 GeV
      data as0  /0.118D0/, r20 /8315.25D0/       !alphas(Mz2)=0.118, LO ref
      data xmin /1.D-5/, nxin /300/, iosp /3/               !x grid, splord
      dimension qq(4), wt(4)
      data qq/1.D0, 2.25D0, 25.D0, 1.D6/
      data wt/1.D0, 1.D0, 2.D0, 1.D0/
      data nqin/100/

C--   scan settings ------------------------------------------------------
      integer nqscan, nxscan
      parameter (nqscan=5, nxscan=41)
      dimension qscan(nqscan)
      data qscan /2.5D0, 5.0D0, 10.0D0, 30.0D0, 91.19D0/     !Q in GeV

C--   Set-up -----------------------------------------------------------
      lun = 6
      call qcinit(lun,' ')
      call gxmake(xmin,1,1,nxin,nx,iosp)
      call gqmake(qq,wt,4,nqin,nq)
      call wtfile(3,'../weights/timelike.wgt')
      call setord(1)
      call setalf(as0,r20)
      iqc = iqfrmq(qmc2)
      iqb = iqfrmq(qmb2)
      call setcbt(0,iqc,iqb,999)

      iq0C = iqfrmq(qmc2)
      call evolfg(13,funcC,def,iq0C,epsC)

      iq0B = iqfrmq(qmb2)
      call evolfg(23,funcB,def,iq0B,epsB)

C--   Scan and dump to file ----------------------------------------------
      open(unit=51,file='kkD0_scan.dat',status='unknown')
      write(51,'(A)') '# Q_GeV  x  xD_D0(x,Q2)  D_D0(x,Q2)'

      do iq = 1, nqscan
        qval = qscan(iq)
        q2   = qval*qval
        do ix = 1, nxscan
C--       log-spaced x from 0.01 to 0.9
          xlog = dlog(0.01D0) +
     +           (dlog(0.9D0)-dlog(0.01D0))*dble(ix-1)/dble(nxscan-1)
          x = dexp(xlog)

          call allfxq(1,x,q2,pdfC,0,1)
          if (q2.ge.qmb2) then
            call allfxq(2,x,q2,pdfB,0,1)
            xDtot = pdfC(4) + pdfB(4)
          else
            xDtot = pdfC(4)
          endif

          write(51,'(F10.4,2X,E13.5,2X,E13.5,2X,E13.5)')
     +          qval, x, xDtot, xDtot/x
        enddo
      enddo

      close(51)
      write(lun,'('' Scan written to kkD0_scan.dat'')')

      end

C     ----------------------------------------------------------------

C     ==========================================
      double precision function funcC(ipdf,x)
C     ==========================================
      implicit double precision (a-h,o-z)
      data xn /0.694D0/, eps /0.101D0/

      funcC = 0.D0
      if (ipdf.eq.8) then
        den   = (1.D0-x)**2 + eps*x
        funcC = xn * x**2 * (1.D0-x)**2 / den**2
      endif

      return
      end

C     ==========================================
      double precision function funcB(ipdf,x)
C     ==========================================
      implicit double precision (a-h,o-z)
      data xn /81.7D0/, alfa /1.81D0/, beta /4.95D0/

      funcB = 0.D0
      if (ipdf.eq.10) then
        funcB = xn * x**(alfa+1.D0) * (1.D0-x)**beta
      endif

      return
      end
