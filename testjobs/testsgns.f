C     ------------------------------------------------------------------
      program testsgns
C     ------------------------------------------------------------------
C--   Compare EVOLSG andf EVSGNS
C     ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C--   ------------------------------------------------------------------
C--   QCDNUM declarations and parameter settings
C--   ------------------------------------------------------------------
C--   LO/NLO/NNLO,  unpol/pol/frag,  starting scale, debug printout
      data iord/3/, itype/1/, q0/2.D0/, idbug/0/

C--   FFNS,VFNS,MFNS and flavour thresholds
      data nfix/6/, q2c/3.D0/, q2b/25.D0/, q2t/1.D11/

C--   Input alphas
      data as0/0.364/, r20/2.D0/

C--   Renormalisation scale and factorisation scale
      data aar/1.D0/, bbr/0.D0/, aaq/1.D0/, bbq/0.D0/

C--   x-grid
      dimension xmin(5), iwt(5)
      data nxsubg/5/, nxin/100/, iosp/3/
      data xmin/1.D-5,0.2D0,0.4D0,0.6D0,0.75D0/
      data iwt / 1, 2, 4, 8, 16 /

C--   mu2-grid
      dimension qq(2),wt(2)
      data qq/2.D0,1.D4/, wt/1.D0,1.D0/, nqin/60/

C--   Input parton distributions and flavour composition of these
      dimension def(-6:6,12)
      data def  /
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,   !dval
     + 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,   !uval
     + 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,   !sval
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,   !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,   !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !sbar
     + 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,   !cval
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !cbar
     + 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,   !bval
     + 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !bbar
     +-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,   !tval
     + 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. /  !tbar

C--   Weight file
      character*26 fnamw(3)
C--                12345678901234567890123456
      data fnamw /'../weights/unpolarised.wgt',
     +            '../weights/polarised.wgt  ',
     +            '../weights/timelike.wgt   '/

C--   Evolution type definitions for EVSGNS: 1=SG, -1=V, 2=NS+, -2=NS-
      dimension isns(12)
C--               1   2   3   4   5   6   7   8   9  10  11  12
      data isns/  1,  2,  2,  2,  2,  2, -1, -2, -2, -2, -2, -2 /

C--   Start functions
      external func1,func2

C--   ------------------------------------------------------------------
C--   QCDNUM set-up
C--   ------------------------------------------------------------------

      lun    = 6
      call QCINIT(lun,' ')                               !initialisation
      call GXMAKE(xmin,iwt,nxsubg,nxin,nx,iosp)                  !x-grid
      call GQMAKE(qq,wt,2,nqin,nq)                             !mu2-grid
      call WTFILE(itype,fnamw(itype))                !make weight tables
      call SETORD(iord)                                   !LO, NLO, NNLO
      call SETALF(as0,r20)                                 !input alphas
      iqc  = IQFRMQ(q2c)                                !charm threshold
      iqb  = IQFRMQ(q2b)                               !bottom threshold
      if(nfix.ge.0) then
        call SETCBT(nfix,iqc,iqb,0)              !thresholds in the VFNS
      else
        call MIXFNS(abs(nfix),q2c,q2b,0.D0)      !thresholds in the MFNS
      endif
      call SETABR(aar,bbr)                        !renormalisation scale
      iq0  = IQFRMQ(q0)                                  !starting scale

      call GETINT('lunq',lunout)             !output logical unit number

C--   Print settings
      write(lunout,'('' xmin ='',1PE11.2,3X,''subgrids  ='',I4)')
     &  xmin(1), nxsubg
      write(lunout,'('' qmin ='',1PE11.2,3X,''qmax      ='',1PE11.2)')
     &  qq(1), qq(2)
      write(lunout,'('' nfix ='',I4,10X,''q2c, b, t ='',3(1PE11.2))')
     &  nfix, q2c, q2b, q2t
      write(lunout,'('' alfa ='',1PE11.2,3X,''r20       ='',1PE11.2)')
     &  as0, r20
      write(lunout,'('' aar  ='',1PE11.2,3X,''bbr       ='',1PE11.2)')
     &  aar, bbr
      write(lunout,'('' q20  ='',1PE11.2)') q0
      write(lunout,'(/)')

C--   ------------------------------------------------------------------
C--   Evolutions
C--   ------------------------------------------------------------------

C--   Evolfg
      iset1  = 1                                       !evolve pdf set 1
      jtype  = 10*iset1+itype
      call EVOLFG(jtype,func1,def,iq0,eps)

C--   Evsgns func2 gets the start values at iq0 from iset1
      call QSTORE('write',1,dble(iset1))    !pass iset1 via common block
      call QSTORE('write',2,dble(iq0))        !pass iq0 via common block

C--   Evsgns
      iset2 = 2                                        !evolve pdf set 2
      jtype = 10*iset2+itype
      call setint('edbg',idbug)
      call EVSGNS(jtype,func2,isns,12,iq0,eps)

C--   ------------------------------------------------------------------
C--   Compare
C--   ------------------------------------------------------------------

      write(lunout,'(/'' Compare EVOLFG and EVSGNS''/)')
      id1 = 0
      id2 = 12
      iqq = 0
      ixm = 1
      do id = id1,id2
        call compa(lunout,iset1,iset2,id,iqq,ixm)
      enddo
      write(lunout,'(/)')

      end

C--   ------------------------------------------------------------------
C--   Comparison routine
C--   ------------------------------------------------------------------

C     =============================================
      subroutine compa(lun,jset1,jset2,id,jq,ixmin)
C     =============================================

C--   Compare basis pdfs
C--
C--   jset1,2   (in) : pdf sets to compare
C--   id        (in) : basis pdf id [0,12]
C--   jq        (in) : 0 = all mu2; != 0 --> compare for this jq only
C--   ixmin     (in) : print range [1,ixmin] when jq != 0

      implicit double precision (a-h,o-z)

      call grpars(nx,xmi,xma,nq,qmi,qma,iosp)

      if(jq.eq.0) then
        iq1   = 1
        iq2   = nq
        jxmin = 0
      else
        iq1   = jq
        iq2   = jq
        jxmin = ixmin
      endif

      dif = 0.D0
      do iq = iq1,iq2
        do ix = 1,nx
          val1 = bvalij(jset1,id,ix,iq,1)
          val2 = bvalij(jset2,id,ix,iq,1)
          difi = val1-val2
          dif  = max(dif,abs(difi))
          if(ix.le.jxmin) write(lun,'(2I3,2F10.4,E13.3,2X,3I3)')
     +    id,ix,val1,val2,difi,jq,jset1,jset2
        enddo
      enddo
      if(jxmin.eq.0) then
        write(lun,'('' id = '',I2,'' dif = '',E15.5)') id,dif
      endif

      return
      end

C--   ------------------------------------------------------------------
C--   Pdf start functions
C--   ------------------------------------------------------------------

C     =======================================
      double precision function func1(ipdf,x)                !for EVOLFG
C     =======================================

C--   This function sets heavy quark input xt = xb = xc = xs

      implicit double precision (a-h,o-z)

                     func1 = 0.D0
      if(ipdf.eq. 0) func1 = xglu(x)
      if(ipdf.eq. 1) func1 = xdnv(x)
      if(ipdf.eq. 2) func1 = xupv(x)
      if(ipdf.eq. 3) func1 = 0.D0
      if(ipdf.eq. 4) func1 = xdbar(x)
      if(ipdf.eq. 5) func1 = xubar(x)
      if(ipdf.eq. 6) func1 = xsbar(x)
      if(ipdf.eq. 7) func1 = 0.D0
      if(ipdf.eq. 8) func1 = xsbar(x)
      if(ipdf.eq. 9) func1 = 0.D0
      if(ipdf.eq.10) func1 = xsbar(x)
      if(ipdf.eq.11) func1 = 0.D0
      if(ipdf.eq.12) func1 = xsbar(x)

      return
      end

C     =======================================
      double precision function func2(ipdf,x)                !for EVSGNS
C     =======================================

C--   This function takes the input pdfs from iset

      implicit double precision (a-h,o-z)

      save iset,iq0

      func2 = 0.D0

      if(ipdf.eq.-1) then
C--     Read pdf set and input scale
        call QSTORE('read',1,val)
        iset = int(val)
        call QSTORE('read',2,val)
        iq0  = int(val)
      else
C--     Get start pdf from iset
        func2 = BVALIJ(iset,ipdf,IXFRMX(x),iq0,1)
      endif

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







