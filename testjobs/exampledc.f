C     ------------------------------------------------------------------
      program example
C     ------------------------------------------------------------------
C--   25-03-15  Basic QCDNUM example job with datacards
C     ------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      data x/1.D-3/, q/1.D3/, qmz2/8315.25D0/             !output scales
      dimension pdf(-6:6)                                        !pdfout

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

      external mycards
      common /pass/ pdef(-6:6,12)

C--   Initialise
      lun    = 6
      lunout = abs(lun)
      call qcinit(lun,' ')

C--   Copy def to common block
      do j = 1,12
        do i = -6,6
          pdef(i,j) = def(i,j)
        enddo
      enddo

C--   Book EVOLFG key
      call qcbook( 'Add', 'EVOLFG' )

C--   Dummy read for printout
      call qcards( mycards, '../dcards/example.dcards', -1)
C--   Process the datacards
      call qcards( mycards, '../dcards/example.dcards',  0)

C--   Results
      call allfxq(1,x,q,pdf,0,1)
      csea = 2.D0*pdf(-4)
      asmz = asfunc(qmz2,nfout,ierr)
      write(lunout,'('' x, q, CharmSea ='',3E13.5)') x,q,csea
      write(lunout,'('' as(mz2)        ='', E13.5)') asmz

      end
      
C     ------------------------------------------------------------------

C     =======================================================
      subroutine mycards ( key, nk, pars, np, fmt, nf, ierr )
C     =======================================================

C--   ierr (out) : 1 = parameter list read error
C--                2 = cannot process datacard
C--                3 = unknown key

      implicit double precision (a-h,o-z)

      character*(*) key, pars, fmt
      external func
      common /pass/ pdef(-6:6,12)

      idum = nk + np + nf     !avoid compiler warning
      
      if(key .eq. 'EVOLFG') then
        read(unit=pars,fmt=fmt,err=100,end=100) itype,q0
        call evolfg(itype,func,pdef,iqfrmq(q0),eps)
      else
        ierr = 3
      endif
      
      return
       
 100  ierr = 1
      return
      end

C     ------------------------------------------------------------------

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







