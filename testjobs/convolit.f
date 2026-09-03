C     ----------------------------------------------------------------
      program convolit
C     ----------------------------------------------------------------
C--   Test convolution engine against Gauss integration 
C--   VFNS w/backward NNLO evolution on multiple grids 
C     ----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      
      external func                                !input parton dists
      dimension def(-6:6,12)                       !flavor composition
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
     + 52*0.    /

C--   Weight file
      character*26 fnamw(3)
C--                12345678901234567890123456
      data fnamw /'../weights/unpolarised.wgt',
     +            '../weights/polarised.wgt  ',
     +            '../weights/timelike.wgt   '/

      external mycards
      common /evol/ q20, itype

      real time1, time2
      
      dimension itypes(6)
      data itypes /0,0,0,0,0,0/
      parameter ( nwds = 9956 )
      dimension wx(nwds)

      dimension idw(4)
      data idw   /0,0,0,0/
      
      dimension gxp(3),dbl(3),gxg(3)

      dimension coef(3:6)
      data coef /4*1.D0/

      external cvol,pvol,dfun,qpdvol,afun,achi
      common /pass/ x, qmu2, iset, ipdf
      
C--   lun = 6 stdout, -6 stdout w/out banner page
      lun    = 6
      lunout = abs(lun) 
      
      call qcinit(lun,' ')

C--   Book user datacard (not used here)
      call qcBook('Add','EVPARS')

C--   Dummy read for printout
      call qcards( mycards, '../dcards/convolit.dcards', -1)
C--   Process the datacards
      call qcards( mycards, '../dcards/convolit.dcards',  0)

      iq0 = iqfrmq(q20)
      call cpu_time(time1)                   !start timing
      call evolfg(itype,func,def,iq0,eps)    !evolve all pdf's
      call cpu_time(time2)                   !end timing

C--   Book and fill weight tables
      itypes(1) = 2                          !two type-1 tables (x)
      itypes(2) = 2                          !two type-2 tables (x,nf)
      call MakeTab(wx,nwds,itypes,0,0,isetx,nwrds)     !book tables

      id    = 1000*isetx+101
      call MakeWtX(wx,id)                    !weights for f * f

      id    = 1000*isetx+102
      call MakeWtA(wx,id,afun,achi)          !convolution kernel K

      idPGQ = idspfun('PGQ',1,1)             !Pgq in QCDNUM memory
      id    = 1000*isetx+201
      call CopyWgt(wx,idPGQ,id,0)            !copy Pgq to weight store

      ida   = 1000*isetx+201
      idb   = 1000*isetx+102
      idc   = 1000*isetx+202
      call WcrossW(wx,ida,idb,idc,0)         !weights for Pgq * K

      iset  = 1                              !unpolarised
      ipdf  = 0                              !gluon
      ix    = 1
      iq    = iq0
      x     = xfrmix(ix)                     !interpolation point
      qmu2  = qfrmiq(iq)                     !interpolation point
      
C--   Gauss integration 
      gxp(1) = dmb_gauss(pvol,x,1.D0,1.D-5)           !Pgq * gluon
      dbl(1) = dmb_gauss(qpdvol,x,1.D0,1.D-5)         !K * Pgq * gluon
      gxg(1) = dmb_gauss(cvol,x,1.D0,1.D-5)           !gluon * gluon

C--   Slow convolution engine
      jpdf   = iPdfTab(iset,ipdf)
      id     = 1000*isetx+201
      gxp(2) = FcrossK(wx,id,iset,jpdf,ix,iq)         !Pgq * gluon

      id     = 1000*isetx+202
      dbl(2) = FcrossK(wx,id,iset,jpdf,ix,iq)         !K * Pgq * gluon

      id     = 1000*isetx+101
      gxg(2) = FcrossF(wx,id,iset,jpdf,jpdf,ix,iq)    !gluon * gluon

C--   Fast convolution engine
      call FastIni(x,qmu2,1,1)               !pass x,q to fast engine
      jpdf   = iPdfTab(iset,ipdf)
      call FastEpm(iset,jpdf,1)              !1 = gluon fast buffer

      idw(1) = 1000*isetx+201                !Pgq weight table id
      call FastFxK(wx,idw,1,2)               !2 = Pgq * gluon
      call FastFxq(2,gxp(3),1)               !Interpolate

      idw(1) = 1000*isetx+201                !Pgq weight table id
      call FastFxK(wx,idw,1,-2)              !2 = Pgq * gluon
      idw(1) = 1000*isetx+102                !K weight table id
      call FastFxK(wx,idw,2,3)               !3 = K * Pgq * gluon
      call FastFxq(3,dbl(3),1)               !Interpolate

      id = 1000*isetx+101
      call FastFxF(wx,id,1,1,2)              !2 = gluon * gluon
      call FastFxq(2,gxg(3),1)               !Interpolate

C--   Compare
      write(lunout,'('' gauss, qcdnum, fast = '',3F14.5)') gxp
      write(lunout,'('' gauss, qcdnum, fast = '',3F14.5)') dbl
      write(lunout,'('' gauss, qcdnum, fast = '',3F14.5)') gxg
      
      write(lunout,'(/''time spent:'',F7.4)') time2-time1

      end
      
C     ----------------------------------------------------------------

C     =================================
      double precision function cvol(z)
C     =================================

C--   Gluon convolution

      implicit double precision (a-h,o-z)
      
      common /pass/ x, qmu2, iset, ipdf
      
      r    = x/z
      cvol = bvalxq(iset,ipdf,z,qmu2,1) * 
     +       bvalxq(iset,ipdf,r,qmu2,1) / (z)
      
      return
      end
      
C     ----------------------------------------------------------------

C     =================================
      double precision function pvol(z)
C     =================================

C--   Gluon convolution with PGQ

      implicit double precision (a-h,o-z)
      
      common /pass/ x, qmu2, iset, ipdf
      
      r    = x/z
      pvol = bvalxq(iset,ipdf,z,qmu2,1) * r * dqcP0gfa(r,3) / (z)
      
      return
      end
 
C     ----------------------------------------------------------------

C     ===================================
      double precision function qpdvol(z)
C     ===================================

C--   Gluon convolution with PGQ and dfun

      implicit double precision (a-h,o-z)
      
      common /pass/ x, qmu2, iset, ipdf
      
      r      = x/z
      qpdvol = bvalxq(iset,ipdf,z,qmu2,1) * pdconv(r) / (z)
      
      return
      end

           
C     ===================================      
      double precision function pdconv(x)
C     ===================================

      implicit double precision (a-h,o-z)
      
      external pdfun
      
      common /pass2/ y 
      
      y      =  x
      pdconv =  dmb_gauss(pdfun,x,1.D0,1.D-5)
      
      return
      end
     
C     ==================================       
      double precision function pdfun(z)
C     ==================================

C--   Dfun convolution with PGQ

      implicit double precision (a-h,o-z)
      
      common /pass2/ y
      
      r     = y/z
      pdfun = z * dfun(z) * r * dqcP0gfa(r,3) /z
      
      return
      end
                  
      
C     =================================      
      double precision function dfun(x)
C     =================================

C--   Mickey Mouse kernel

      implicit double precision (a-h,o-z)
      
      dfun = x*(1.D0-x)
      
      return
      end
      
C     ----------------------------------------------------------------

C     ====================================
      double precision function achi(qmu2)      
C     ====================================

      implicit double precision (a-h,o-z)
      
      ddum = qmu2
      achi = 1.D0
      
      return
      end                        
     
C     =========================================     
      double precision function afun(x,qmu2,nf)
C     =========================================

      implicit double precision (a-h,o-z)
      
      ddum = qmu2
      idum = nf
      afun = dfun(x)
      
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
      common /evol/ q20, itype

      idum = nk + np + nf     !avoid compiler warning

      if(key .eq. 'EVPARS') then
        read(unit=pars,fmt=fmt,err=100,end=100) itype, q20
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







