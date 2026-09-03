 
C     =================
      program testzmstf
C     =================

C--   Test the zmstf package

      implicit double precision (a-h,o-z)

      real tim1,tim2,tim3,tim4
      
C--   Coarse grid for output
      common /tabout/ xxtab(99),qqtab(99),
     +                xmi,xma,qmi,qma,nxtab,nqtab      
C--   qcdnum input
      dimension pdfin(-6:6,12)
      data pdfin/  
C--    tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--    -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
     +  0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,        !dval
     +  0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,        !uval
     +  0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,        !sval
     +  0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,        !dbar
     +  0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,        !ubar
     +  0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,        !sbar
     +  0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,        !cval
     +  0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,        !cbar
     +  0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,        !bval
     +  0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,        !bbar
     + -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,        !tval
     +  1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./        !tbar
C--   Holds q-grid definition
      dimension qarr(2),warr(2)
C--   Input pdf functions
      external func
C--   Copy pdfs to set 5      
      external mypdfs      
      
C--   Define output distributions
      dimension dnv(-6:6),upv(-6:6),del(-6:6),uds(-6:6)
      dimension pro(-6:6),duv(-6:6)
C--             tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--             -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
      data dnv / 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0. /
      data upv / 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0. /
      data del / 0., 0., 0., 0.,-1., 1., 0., 0., 0., 0., 0., 0., 0. /
      data uds / 0., 0., 0., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0. /
      data pro / 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. /
      data duv / 0., 0., 0., 0.,-1.,-1., 0., 1., 1., 0., 0., 0., 0. /
      
C--   cbt quark masses
      dimension hqmass(3)
      
      dimension xlist(100),qlist(100),flist(100)
      dimension xx(1),qq(1),F2(1),FL(1),FLp(1),xF3(1)

      lun    = 6
      lunout = abs(lun)
      call qcinit(lun,' ')

C--   LO/NLO/NNLO
      iord = 3
C--   VFNS
      nfix = 0        
C--   Define alpha_s
      alf = 0.35D0
      q2a = 2.0D0
      
C--   Coarse x-Q2 grid for printout (defines xmi,xma,qmi,qma)
      call makepgrid
            
C--   More grid parameters
      iosp  = 3
      n_x   = 100
      n_q   = 60
      q0    = 2.D0

C--   QCDNUM calls to pass the values defined above
      call setord(iord)
      call setalf(alf,q2a)
C--   x-mu2 grid
      call gxmake(xmi,1,1,n_x,nxout,iosp)
      qarr(1) = qmi
      qarr(2) = qma
      warr(1) = 1.D0
      warr(2) = 1.D0
      call gqmake(qarr,warr,2,n_q,nqout)
      iq0 = iqfrmq(q0)
C--   Quark masses      
      hqmass(1) =  2.0D0
      hqmass(2) =  5.0D0
      hqmass(3) = 80.0D0
C--   Thresholds
      iqc = iqfrmq(hqmass(1)*hqmass(1)) 
      iqb = iqfrmq(hqmass(2)*hqmass(2)) 
      iqt = iqfrmq(hqmass(3)*hqmass(3))     
      call setcbt(nfix,iqc,iqb,iqt)

C--   Try to read the weight file and create one if that fails
      call wtfile(1,'../weights/unpolarised.wgt')

C--   Try to read the weight file and create one if that fails
      call zmreadw(22,'../weights/zmstf.wgt',nwords,ierr)
      if(ierr.ne.0) then
        call zmfillw(nwords)
        call zmdumpw(22,'../weights/zmstf.wgt')
      endif
      write(lunout,'(/'' ZMSTF: words used ='',I10)') nwords 
      
      call zmwords(nztot,nzuse)
      write(6,'(/'' nztot, nzuse = '',2I10)') nztot,nzuse
      
C--   Scale variation
*      call setabr(2.D0,0.D0)
*      call zmdefq2(2.D0,1.D0) 

C--   Evolve
      call cpu_time(tim1)
      iter = 1
      do i = 1,iter
        call EvolFG(1,func,pdfin,iq0,epsi)
      enddo
      call cpu_time(tim2)
      
C--   Charge squared weighted for proton
      do i = -6,6
        pro(i) = pro(i)/9.D0
      enddo
      
      ichk = 1
      iset = 1
      
      write(lunout,'(/''     mu2       x         uv         '//
     +             'dv        del        uds         gl'')')

      do iq = 1,nqtab,4
        q = qqtab(iq)
        write(lunout,'('' '')')
        do ix = 1,nxtab,4
          x  = xxtab(ix)
          uv = sumfxq(iset,upv,1,x,q,ichk)
          dv = sumfxq(iset,dnv,1,x,q,ichk)
          de = sumfxq(iset,del,1,x,q,ichk)
          ud = sumfxq(iset,uds,1,x,q,ichk)
          gl = fvalxq(iset,  0,x,q,ichk)
          write(lunout,18) q,x,uv,dv,de,ud,gl
        enddo
      enddo
 
      write(lunout,'(/''     mu2       x         pr        '//
     +             'F2         FL        FLp        xF3'')')
     
      call zswitch(iset)

      do iq = 1,nqtab,4
        q  = qqtab(iq)
        qq = q
        write(lunout,'('' '')')
        do ix = 1,nxtab,4
          x  = xxtab(ix)
          xx = x
          pr = sumfxq(iset,pro,1,x,q,ichk)
C--       Ideally, zmstfun should not be called in a loop...
          call zmstfun(2,pro,xx,qq,F2, 1,ichk)
          call zmstfun(1,pro,xx,qq,FL, 1,ichk)
          call zmstfun(4,pro,xx,qq,FLp,1,ichk)
          call zmstfun(3,duv,xx,qq,xF3,1,ichk)
*          call zmslowf(2,pro,xx,qq,F2, 1,ichk)
*          call zmslowf(1,pro,xx,qq,FL, 1,ichk)
*          call zmslowf(4,pro,xx,qq,FLp,1,ichk)
*          call zmslowf(3,duv,xx,qq,xF3,1,ichk)
          write(lunout,18) q,x,pr,F2,FL,FLp,xF3
        enddo
      enddo

      call cpu_time(tim3)

C--   Make a list of interpolation points      
      nlist = 0
      do iq = 1,nqtab,4
        do ix = 1,nxtab,4
          nlist        = nlist+1
          xlist(nlist) = xxtab(ix)
          qlist(nlist) = qqtab(iq)
        enddo
      enddo
C--   ZmStfun with the list of points, instead of calling it in a loop
C--   This is to see how much faster that is (factor 10!)  
      call zmstfun(2,pro,xlist,qlist,flist,nlist,ichk)
      call zmstfun(1,pro,xlist,qlist,flist,nlist,ichk)
      call zmstfun(4,pro,xlist,qlist,flist,nlist,ichk)
      call zmstfun(3,pro,xlist,qlist,flist,nlist,ichk)
      
      call cpu_time(tim4)  

C--   Print alfas 
      write( lunout,'(/''     mu2      alfas'')')
      do iq = 1,nqtab,4
        q = qqtab(iq)
        a = asfunc(q,nf,ierr)
        write(lunout,19) q,a,ierr
      enddo

      write(lunout,'(/''time spent evol:'',F7.4)') (tim2-tim1)/iter
      write(lunout,'( ''time spent stfs:'',F7.4)') (tim3-tim2)
      write(lunout,'( ''time spent stfs:'',F7.4)') (tim4-tim3)
      
  18  format(1X,2(1PE9.1),1X,5(1PE11.4))
  19  format(1X,  1PE9.1 ,1X,  1PE12.5,1X,I5)

      end
      
C     ==============================
      subroutine mypdfs(x,qmu2,xpdf)
C     ==============================

      implicit double precision (a-h,o-z)
      
      dimension xpdf(-6:6)
      
      call fpdfxq(1,x,qmu2,xpdf,1)
*      xpdf(7) = fsnsxq(1,0,x,qmu2,1)
      
      return
      end

      
C     ====================
      subroutine makepgrid
C     ====================

      implicit double precision (a-h,o-z)

      common /tabout/ xxtab(99),qqtab(99),
     +                xmi,xma,qmi,qma,nxtab,nqtab

C--   coarse x - Q2 grid
 
      xxtab( 1) = 1.0D-5
      xxtab( 2) = 2.0D-5
      xxtab( 3) = 5.0D-5
      xxtab( 4) = 1.0D-4
      xxtab( 5) = 2.0D-4
      xxtab( 6) = 5.0D-4
      xxtab( 7) = 1.0D-3
      xxtab( 8) = 2.0D-3
      xxtab( 9) = 5.0D-3
      xxtab(10) = 1.0D-2
      xxtab(11) = 2.0D-2
      xxtab(12) = 5.0D-2
      xxtab(13) = 1.0D-1
      xxtab(14) = 1.5D-1
      xxtab(15) = 2.0D-1
      xxtab(16) = 3.0D-1
      xxtab(17) = 4.0D-1
      xxtab(18) = 5.5D-1
      xxtab(19) = 7.0D-1
      xxtab(20) = 9.0D-1
      nxtab     = 20
      xmi       = xxtab(1)
      xma       = xxtab(nxtab)

      qqtab( 1) = 2.0D0
      qqtab( 2) = 2.7D0
      qqtab( 3) = 3.6D0
      qqtab( 4) = 5.0D0
      qqtab( 5) = 7.0D0
      qqtab( 6) = 1.0D1
      qqtab( 7) = 1.4D1
      qqtab( 8) = 2.0D1
      qqtab( 9) = 3.0D1
      qqtab(10) = 5.0D1
      qqtab(11) = 7.0D1
      qqtab(12) = 1.0D2
      qqtab(13) = 2.0D2
      qqtab(14) = 5.0D2
      qqtab(15) = 1.0D3
      qqtab(16) = 3.0D3
      qqtab(17) = 1.0D4
      qqtab(18) = 4.0D4
      qqtab(19) = 2.0D5
      qqtab(20) = 1.0D6
      nqtab     = 20
      qmi       = qqtab(1)
      qma       = qqtab(nqtab)

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

C     ===================================
      double precision function dnor(b,c)
C     ===================================

      implicit double precision (a-h,o-z)

C--   Int_0^1 x**b (1-x)**c dx

      dnor = dmb_gamma(b+1.d0)*dmb_gamma(c+1.d0)/dmb_gamma(b+c+2.d0)

      return
      end
