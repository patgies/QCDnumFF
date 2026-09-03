C     ----------------------------------------------------------------
      program pdfsets
C     ----------------------------------------------------------------
C--   Play around with pdf sets in memory
C--
C--   set 1 : unpolarised NNLO
C--   set 2 : polarised   NLO
C--   set 5 : unpolarised LO
C--   set 6 : unpolarised NLO
C--   set 7 : unpolarised NNLO
C--   set 8 : unpolarised NNLO + 2 extra pdfs (singlet and proton)
C--
C--   This tests if PDFCPY, PDFEXT, FTABLE, FFLIST are working correctly
C     ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      
      data as0/0.364/, r20/2.D0/, iord/3/, nfin/0/   !alphas, NNLO, VFNS
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
      dimension xmin(5), iwt(5)
      data xmin/1.D-5,0.2D0,0.4D0,0.6D0,0.75D0/                  !x grid
      data iwt / 1, 2, 4, 8, 16 /                                !x grid
      data nxin/100/, iosp/3/                                    !x grid
      dimension qlim(2),wt(2)                                  !mu2 grid
      data qlim/2.D0,1.D4/, wt/1.D0,1.D0/, nqin/60/            !mu2 grid
      data q2c/3.D0/, q2b/25.D0/, q0/2.0/                    !thresh, q0
      dimension proton(-6:6)                                     !proton
      data proton/4.,1.,4.,1.,4.,1.,0.,1.,4.,1.,4.,1.,4./

      dimension xx(5),qq(3)                                      !output
      data xx /1.D-5, 1.D-4, 1.D-3, 1.D-2, 1.D-1/                !output
      data qq /1.D1 , 1.D2 , 1.D3 /                              !output
      dimension xf1(5,3),xf2(5,3),xf8(5,3),gt1(5,3)              !output
      dimension psi(5,3),pns(5,3)                                !output
      dimension xl(15),ql(15),gl1(15)                            !output
      dimension epar(13)                                         !output

      external fimport             !function to read pdfs to be imported
      
C--   Divide proton by 9
      do i = -6,6
        proton(i) = proton(i)/9.D0
      enddo

      call qcinit(6,' ')                                           !init
      do i = -6,6
        call qstore('Write',i+7,proton(i))                 !store proton
      enddo
      call gxmake(xmin,iwt,5,nxin,nx,iosp)                        !xgrid
      call gqmake(qlim,wt,2,nqin,nq)                              !qgrid
      call wtfile(1,'../weights/unpolarised.wgt')         !weights unpol
      call wtfile(2,'../weights/polarised.wgt')           !weights   pol
      call setalf(as0,r20)                                     !input as
      iqc  = iqfrmq(q2c)                                          !charm
      iqb  = iqfrmq(q2b)                                         !bottom
      call setcbt(nfin,iqc,iqb,999)                         !VFNS thresh
      iq0  = iqfrmq(q0)                                     !input scale

      call getint('lunq',lun)                                !output lun

C--   Unpolarised LO, NLO, NNLO in sets 5, 6, 7 ------------------------
      do iord = 1,3
        call setord(iord)
        call evolfg(1,func,def,iq0,eps)
        call pdfcpy(1,4+iord)
      enddo

C--   Print gluon table of set 1 ---------------------------------------
      call ftable(1,proton,0,xx,5,qq,3,gt1,1)    !gluon     set 1
      write(lun,'('' Gluon table  IX = '',I2,4I13)') (i,i=1,5)
      do j = 1,3
        write(lun,'('' IQ ='',I3,5E13.5)') j,(gt1(i,j),i=1,5)
      enddo

C--   Make gluon list of set 1 -----------------------------------------
      k = 0
      do j = 1,3
        do i = 1,5
          k     = k+1
          xl(k) = xx(i)
          ql(k) = qq(j)
        enddo
      enddo
      call fflist(1,proton,0,xl,ql,gl1,15,1)      !gluon    set 1

C--   Check that table = list ------------------------------------------
      diff = -10.D0
      k    = 0
      do j = 1,3
        do i = 1,5
          k    = k+1
          diff = max( diff, abs(gt1(i,j)-gl1(k)) )
        enddo
      enddo
      write(lun,'(/'' This should be small:'',E15.5)') diff

C--   Polarised pdfs NLO in default set 2 ------------------------------
      call setord(2)
      call evolfg(2,func,def,iq0,eps)

C--   Import NNLO unpol (from set 1) into set 8 with 2 extra pdfs ------
C--   Make sure that extpdf uses the params of set 1
      call pushcp                                        !save pars
      call usepar(1)                                     !use pars set 1
      offset = 0.D0
C--   s/r fimport stores proton singlet and nonsinglet in idpdf 7 and 8
      call extpdf(fimport,8,2,offset,eps)
      call pullcp                                          !restore pars

C--   Make proton table of set 1 and 8 ---------------------------------
      call ftable(1,proton,1,xx,5,qq,3,xf1,1)          !proton     set 1
      call ftable(8,proton,1,xx,5,qq,3,xf8,1)          !proton     set 8
      call ftable(8,proton,13,xx,5,qq,3,psi,1)         !singlet    set 8
      call ftable(8,proton,14,xx,5,qq,3,pns,1)         !nonsinglet set 8

C--   See if set 1 is properly imported into set 8  --------------------
      diff = -10.D0
      do j = 1,3
        do i = 1,5
          diff = max( diff, abs(xf8(i,j)-xf1(i,j)) )
        enddo
      enddo
      write(lun,'('' This should be small:'',E15.5)') diff

C--   See if proton singlet and nonsinglet are properly stored in 8 ----
      diff = -10.D0
      do j = 1,3
        do i = 1,5
          diff = max( diff, abs(xf8(i,j)-psi(i,j)-pns(i,j)) )
        enddo
      enddo
      write(lun,'('' This too            :'',E15.5)') diff

C--   Now print proton table of set 2 ----------------------------------
      call ftable(2,proton,1,xx,5,qq,3,xf2,1)              !proton set 2
      write(lun,'(/'' Proton table IX = '',I2,4I13)') (i,i=1,5)
      do j = 1,3
        write(lun,'('' IQ ='',I3,5E13.5)') j,(xf2(i,j),i=1,5)
      enddo

C--   Print some parameter values for all data sets --------------------
      write(lun,'(/''    ISET    IORD    ITYP    NPDF     KEY'')')
      do iset = 1,8
        npdfs  = nptabs(iset)
        if(npdfs.ne.0) then                         !check if set exists
          call cpypar(epar,13,iset)
          iorder = int(epar( 1))
          ipdfs  = int(epar(13))
          key    = keypar(iset)
          write(lun,'(5I8)') iset, iorder, ipdfs, npdfs, key
        endif
      enddo

      end

C     ------------------------------------------------------------------

C     ======================================================
      double precision function fimport( ipdf, x, q, first )
C     ======================================================

      implicit double precision (a-h,o-z)
      dimension proton(-6:6)
      save proton
      logical first

      if(first) then
        do i = -6,6
          call qstore('Read',i+7,proton(i))
        enddo
      endif

      if(ipdf.ge.-6 .and. ipdf.le.6) then
        fimport = fvalxq(1,ipdf,x,q,1)               ! -6:6 = qbar, g, q
      elseif(ipdf.eq.7) then
        fimport = sumfxq(1,proton,2,x,q,1)          ! 7 = proton singlet
      elseif(ipdf.eq.8) then
        fimport = sumfxq(1,proton,3,x,q,1)       ! 8 = proton nonsinglet
      else
        fimport = 0.D0                             ! any other pdf index
      endif

      return
      end

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







