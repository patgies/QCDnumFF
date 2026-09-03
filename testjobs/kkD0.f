C     ----------------------------------------------------------------
      program kkD0
C     ----------------------------------------------------------------
C--   Time-like (fragmentation function) evolution of the Kniehl-Kramer
C--   D0 fragmentation function, LO, mu0=mQ convention.
C--   B.A. Kniehl and G. Kramer, "Charmed-Hadron Fragmentation Functions
C--   from CERN LEP1 Revisited", DESY 06-102 (arXiv:hep-ph/0607306).
C--
C--   The paper has TWO separate non-perturbative inputs:
C--     - charm-quark FF (Peterson form) at mu0 = mc
C--     - bottom-quark FF (power law)    at mu0 = mb
C--   with everything else zero at that flavour's own threshold.
C--
C--   Because DGLAP evolution is linear, the full result is obtained
C--   by evolving each input separately on the same grid and adding:
C--     D_total(x,Q) = D_c-evolved(x,Q) + D_b-evolved(x,Q)   for Q >= mb
C--     D_total(x,Q) = D_c-evolved(x,Q)                      for mc <= Q < mb
C     ----------------------------------------------------------------
      implicit double precision (a-h,o-z)

      external funcC, funcB                          !input FF functions
      dimension def(-6:6,12)                          !flavour composition

C--   Standard QCDNum flavour basis (shared by both evolutions):
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
      data xmin /1.D-5/, nxin /150/, iosp /3/               !x grid, splord
C--   mu2 grid starts below qmc2 so that setcbt has grid points below
C--   the charm threshold (required by QCDNUM)
      dimension qq(4), wt(4)
      data qq/1.D0, 2.25D0, 25.D0, 1.D6/
      data wt/1.D0, 1.D0, 2.D0, 1.D0/
      data nqin/100/

C--   Set-up -----------------------------------------------------------
      lun = 6
      call qcinit(lun,' ')                                    !initialize
      call gxmake(xmin,1,1,nxin,nx,iosp)                          !x-grid
      call gqmake(qq,wt,4,nqin,nq)                             !mu2-grid
      call wtfile(3,'../weights/timelike.wgt')     !itype=3: time-like wgts
      call setord(1)                                        !LO evolution
      call setalf(as0,r20)                                   !input alphas
      iqc = iqfrmq(qmc2)
      iqb = iqfrmq(qmb2)
      call setcbt(0,iqc,iqb,999)                    !dynamic VFNS thresholds

C--   Evolution A: charm-quark input at mu0 = mc -------------------------
      iq0C = iqfrmq(qmc2)
      call evolfg(13,funcC,def,iq0C,epsC)      !itype=3(t-like),iset=1 -> jset=13

C--   Evolution B: bottom-quark input at mu0 = mb -----------------------
      iq0B = iqfrmq(qmb2)
      call evolfg(23,funcB,def,iq0B,epsB)      !itype=3(t-like),iset=2 -> jset=23

C--   Get results  --------------------------------------------------------
      x  = 0.1D0
      q2 = 100.D0                                   !an example scale > qmb2

      call allfxq(1,x,q2,pdfC,0,1)
      call allfxq(2,x,q2,pdfB,0,1)

      if (q2.ge.qmb2) then
        xDtot = pdfC(4) + pdfB(4)          !c-quark component: A + B
      else
        xDtot = pdfC(4)                    !only the charm evolution applies
      endif

      call getint('lunq',lunout)
      write(lunout,'('' x, Q2, x*D_c^D0(x,Q2) ='',2E13.5,E14.5)')
     +      x, q2, xDtot
      write(lunout,'('' -> D_c^D0(x,Q2) (divide by x) ='',E14.5)')
     +      xDtot/x

      end

C     ----------------------------------------------------------------

C     ==========================================
      double precision function funcC(ipdf,x)
C     ==========================================
C--   Charm-quark Peterson input D_c(x,mc^2), D0 meson, LO
C--   (Table I of hep-ph/0607306: N=0.694, epsilon=0.101)
      implicit double precision (a-h,o-z)
      data xn /0.694D0/, eps /0.101D0/

      funcC = 0.D0
C--   only the cbar column (ipdf=8) is populated; def then gives c=cbar=D_c
      if (ipdf.eq.8) then
        den   = (1.D0-x)**2 + eps*x
        funcC = xn * x**2 * (1.D0-x)**2 / den**2      !x * D_c(x,mc^2)
      endif

      return
      end

C     ==========================================
      double precision function funcB(ipdf,x)
C     ==========================================
C--   Bottom-quark power-law input D_b(x,mb^2), D0 meson, LO
C--   (Table I of hep-ph/0607306: N=81.7, alpha=1.81, beta=4.95)
      implicit double precision (a-h,o-z)
      data xn /81.7D0/, alfa /1.81D0/, beta /4.95D0/

      funcB = 0.D0
C--   only the bbar column (ipdf=10) is populated; def then gives b=bbar=D_b
      if (ipdf.eq.10) then
        funcB = xn * x**(alfa+1.D0) * (1.D0-x)**beta   !x * D_b(x,mb^2)
      endif

      return
      end
