C     ----------------------------------------------------------------
      program longlist
C     ----------------------------------------------------------------
C--   Interpolate arbitrary long list of pdfs using the fast engine.
C--   This illustrates how to circumvent the mpt0 limit on the max
C--   number of interpolations in one call.
C--   The result is compared with the QCDNUM routine FTABLE.
C     ----------------------------------------------------------------

      implicit double precision (a-h,o-z)
      
C--   ----------------------------------------------------------------
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

      external mycards
      common /pass/ pdef(-6:6,12)

      real tim1, tim2

C--   ----------------------------------------------------------------
C--   Define here which pdf to interpolate (dvalence in this example)
      dimension pdf(-6:6)
C--             tbar  bbar  cbar  sbar  ubar  dbar  glue
C--               -6    -5    -4    -3    -2    -1
      data pdf /0.D0, 0.D0, 0.D0, 0.D0, 0.D0,-1.D0, 0.D0,
C--             down    up   str   chm   bot   top 
C--                1     2     3     4     5     6      
     +          1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /
     
C--   ----------------------------------------------------------------     
C--   Define here the output grid to be filled with pdf values            
      parameter( nxu = 100 )
      dimension xgu(nxu)
      parameter( nqu = 100 )
      dimension qgu(nqu)
      dimension ffu(nxu,nqu)

C--   These are lists of x, mu2 and pdf values (to be filled by mypdf)
      dimension xlist(nxu*nqu), qlist(nxu*nqu), flist(nxu*nqu)
      
C--   ----------------------------------------------------------------      
C--   Initialise
      call qcinit(6,' ')

C--   Copy def to common block
      do j = 1,12
        do i = -6,6
          pdef(i,j) = def(i,j)
        enddo
      enddo

C--   Book EVOLFG key
      call qcbook( 'Add', 'EVOLFG' )

C--   Dummy read for printout
      call qcards( mycards, '../dcards/longlist.dcards', -1)
C--   Process the datacards
      call qcards( mycards, '../dcards/longlist.dcards',  0)

      call getint('lunq',lunout)
      call grpars(nx,xmin,xmax,nq,qmin,qmax,iosp)

C--   ----------------------------------------------------------------      
C--   Now generate the output grid (logarithmic)
      xmil = log(xfrmix(1))
      xmal = 0.D0
      bw   = (xmal-xmil)/nxu
      do i = 1,nxu
        xgu(i) = exp(xmil + (i-0.5)*bw)
      enddo   
      qmil = log(qfrmiq(1))
      qmal = log(qfrmiq(nq))
      bw   = (qmal-qmil)/nqu
      do i = 1,nqu
        qgu(i) = exp(qmil + (i-0.5)*bw)
      enddo
      
C--   Make the list of interpolation points
      nlist = 0
      do iq = 1,nqu
        do ix = 1,nxu
          nlist = nlist+1
          xlist(nlist) = xgu(ix)
          qlist(nlist) = qgu(iq)
        enddo
      enddo

C--   Yes/no check boundaries
      ichk = 1
      
C--   ----------------------------------------------------------------
C--   Interpolation
      call cpu_time(tim1)
      call mypdfs(1,pdf,xlist,qlist,flist,nlist,ichk)            !mypdfs
      call cpu_time(tim2)
      write(lunout,'('' mypdfs cpu time       = '',F13.5)') tim2-tim1
      call cpu_time(tim1)
      call ftable(1,pdf,1,xgu,nxu,qgu,nqu,ffu,ichk)              !qcdnum
      call cpu_time(tim2)
      write(lunout,'('' ftable cpu time       = '',F13.5)') tim2-tim1
      
C--   ----------------------------------------------------------------  
C--   Now check if OK 
      dmax  = -1.D0
      nlist = 0
      do iq = 1,nqu
        do ix = 1,nxu
          nlist = nlist+1
          dmax = max(dmax,abs(flist(nlist)-ffu(ix,iq)))
        enddo
      enddo
      write(lunout,'(/'' This should be small ...'',E13.5)') dmax
                       
      end
      
C     ----------------------------------------------------------------      
      
C     ========================================      
      subroutine mypdfs(iset,def,x,q,f,n,ichk)
C     ========================================

C--   Interpolation of linear combination of pdfs using fast engine.
C--   The number of interpolations can be larger than mpt0 since
C--   this routine autmatically buffers in chunks of mpt0 words.
C--
C--   iset       (in)   pdf set [1-9]
C--   def(-6:6)  (in)   coefficients, ..., sb, ub, db, g, d, u, s, ...
C--                     for gluon  set def(0) = non-zero
C--                     for quarks set def(0) = 0.D0   
C--   x          (in)   list of x-points
C--   q          (in)   list of mu2 points
C--   f          (out)  list of interpolated pdfs
C--   n          (in)   number of items in the list
C--   ichk       (in)   0/1  no/yes check grid boundary  
      
      implicit double precision (a-h,o-z)      
      
      dimension def(-6:6), coef(0:12,3:6)
      dimension x(*), q(*), f(*)


C--   Check if n is OK
      if(n.le.0) stop 'MYPDFS: n.le.0 --> STOP'
      
      call setUmsg('MYPDFS')

C--   Translate the coefficients def(-6:6) from flavour space to 
C--   singlet/non-singlet space. These coefs are n_f dependent.
       
      if(def(0).eq.0.D0) then                !quarks
        do nf = 3,6
          coef(0,nf) = 0.D0                  
          call efromqq(def, coef(1,nf), nf) 
        enddo
      else                                   !gluon
        do nf = 3,6
          coef(0,nf) = def(0)
          do i = 1,12
            coef(i,nf) = 0.D0 
          enddo  
        enddo
      endif  
      
C--   Fill output array f in batches of nmax (=mpt0) words
      call getint('mpt0',nmax)
*      nmax  = 10
      nlast = 0
      ntodo = min(n,nmax)
      do while( ntodo.gt.0 )
        i1    = nlast+1
        call interpolate(iset, coef, x(i1), q(i1), f(i1), ntodo, ichk)
        nlast = nlast+ntodo
        ntodo = min(n-nlast,nmax)
      enddo
        
      call clrUmsg
      
      return
      end

C     ==============================================
      subroutine interpolate(iset,coef,x,q,f,n,ichk)
C     ==============================================

C--   Interpolate linear combination of pdfs in internal memory
C--
C--   iset    (in) : pdf set 1=unpol, 2=pol, 3=timelike, 5-9=external
C--   coef    (in) : nf dependent coefficients of the linear combibnation
C--   x,q     (in) : list of interpolation points
C--   f      (out) : list of interpolated values
C--   n       (in) : number of items in x,q,f (always less than mpto)

      implicit double precision (a-h,o-z)

      dimension coef(0:12,3:6)
      dimension x(*), q(*), f(*)

      call fastini(x,q,n,ichk)
      call fastsum(iset,coef,-1)           !fill sparse buffer 1
      call fastfxq(1,f,n)

      return
      end
      
C     ----------------------------------------------------------------

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

