
c These are the functions that give the asymptotic dependence of the
c coefficient functions with the appropriate factors. xi = mq2/m2 (Q^2/m2) 
c If xi is small, the regular routines have convergence
c problems and we take the limit. (not anymore after code update 03/06/96).

c     ==========================================
      double precision function Gfun_L(dummy,xi)
c     ==========================================

c     Longitudinal: equation (19) in PLB347 (1995) 143 - 151
c     This function is called asymp_l in the original code.

      implicit none
      double precision xilast, store
      double precision dummy, ddum
      double precision xi, pi, term1
      double precision fii, fjj
*     double precision fii_lim, fjj_lim
      parameter (pi = 3.14159265359d0)

      save xilast, store

      data xilast, store /0.D0, 0.D0/

      ddum = dummy !avoid compiler warning
      
      if(xi.eq.xilast) then
        Gfun_L = store
        return
      endif

*     term1 = 1.d0/(1.d0 + 0.25d0*xi)

*     if (xi .le. 1.d-1) then
*        Gfun_L = 1.d0/6.d0/pi*(-4.d0/3.d0*term1 + 
*    #        (1.d0  - 1.d0/6.d0*term1)*fjj_lim(xi) -
*    #        2.d0* (-1.d0/3.d0 + xi/15.d0 - xi**2/70.d0) +
*    #        0.25d0*term1*fii_lim(xi) -
*    #        3.d0* (1.d0/3.d0 - xi/10.d0 + 11.d0*xi**2/420.d0))
*     else
*        Gfun_L = 1.d0/6.d0/pi*(4.d0/xi - 4.d0/3.d0*term1 
*    #        + (1.d0 - 2.d0/xi - 1.d0/6.d0*term1)*fjj(xi)
*    #        - (3.d0/xi + 0.25d0*term1)*fii(xi))
*     endif

      term1   = 1.d0/(1.d0 + 0.25d0*xi)

      Gfun_L  = 1.d0/6.d0/pi*(4.d0/xi - 4.d0/3.d0*term1 
     #     + (1.d0 - 2.d0/xi - 1.d0/6.d0*term1)*fjj(xi)
     #     - (3.d0/xi + 0.25d0*term1)*fii(xi))

      xilast = xi
      store  = Gfun_L

      return
      end

c     ==========================================
      double precision function Gfun_T(dummy,xi)
c     ==========================================

c     Transverse: equation (20) in PLB347 (1995) 143 - 151
c     This function is called asymp_t in the original code.

      implicit none
      double precision xilast, store
      double precision dummy, ddum
      double precision xi, pi, term1
      double precision fii, fjj
*     double precision fii_lim, fjj_lim
      parameter (pi = 3.14159265359d0)

      save xilast, store

      data xilast, store /0.D0, 0.D0/
      
      ddum = dummy !avoid compiler warning      

      if(xi.eq.xilast) then
        Gfun_T = store
        return
      endif

*     term1 = 1.d0/(1.d0 + 0.25d0*xi)

*     if (xi .le. 1.d-1) then
*        Gfun_T = 1.d0/6.d0/pi*(4.d0/3.d0*term1 + (7.d0/6.d0 +
*    #        1.d0/6.d0*term1)*fjj_lim(xi) + 1/3.d0*
*    #        (-1.d0/3.d0 + xi/15.d0 - xi**2/70.d0) +
*    #        (1.d0 + 0.25d0*term1)*fii_lim(xi) + 2.d0*
*    #        (1.d0/3.d0 - xi/10.d0 + 11.d0*xi**2/420.d0))
*     else
*        Gfun_T = 1.d0/6.d0/pi*(-2.d0/3.d0/xi + 4.d0/3.d0*term1
*    #        + (7.d0/6.d0 + 1.d0/3.d0/xi + 1.d0/6.d0*term1)*fjj(xi)
*    #        + (1.d0 + 2.d0/xi + 0.25d0*term1)*fii(xi))
*     endif

      term1   = 1.d0/(1.d0 + 0.25d0*xi)

      Gfun_t  = 1.d0/6.d0/pi*(-2.d0/3.d0/xi + 4.d0/3.d0*term1
     #     + (7.d0/6.d0 + 1.d0/3.d0/xi + 1.d0/6.d0*term1)*fjj(xi)
     #     + (1.d0 + 2.d0/xi + 0.25d0*term1)*fii(xi))

      xilast = xi
      store  = Gfun_T

      return
      end

c     ==========================================
      double precision function Gbar_L(dummy,xi)
c     ==========================================

c     Longitudinal mass factorization: (21) in PLB347 (1995) 143 - 151
c     This function is called asympbar_l in the original code.

      implicit none
      double precision xilast, store
      double precision dummy, ddum
      double precision xi, pi, term1
      double precision fjj
*     double precision fjj_lim
      parameter (pi = 3.14159265359d0)

      save xilast, store

      data xilast, store /0.D0, 0.D0/

      ddum = dummy !avoid compiler warning
      
      if(xi.eq.xilast) then
        Gbar_L = store
        return
      endif

*     term1 = 1.d0/(1.d0 + 0.25d0*xi)

*     if (xi .le. 1.d-1) then
*        Gbar_L = 1.d0/6.d0/pi*(0.5d0*term1 + 
*    #        0.25d0*term1*fjj_lim(xi) +
*    #        3.d0* (-1.d0/3.d0 + xi/15.d0 - xi**2/70.d0))
*     else
*        Gbar_L = 1.d0/6.d0/pi*(-6.d0/xi + 0.5d0*term1
*    #        + (3.d0/xi + 0.25d0*term1)*fjj(xi))
*     endif

      term1   = 1.d0/(1.d0 + 0.25d0*xi)

      Gbar_L  = 1.d0/6.d0/pi*(-6.d0/xi + 0.5d0*term1
     #     + (3.d0/xi + 0.25d0*term1)*fjj(xi))

      xilast = xi
      store  = Gbar_L

      return
      end

c     ==========================================
      double precision function Gbar_T(dummy,xi)
c     ==========================================

c     transverse mass factorization: (22) in PLB347 (1995) 143 - 151
c     This function is called asympbar_t in the original code.

      implicit none
      double precision xilast, store
      double precision dummy, ddum
      double precision xi, pi, term1
      double precision fjj
*     double precision fjj_lim
      parameter (pi = 3.14159265359d0)

      save xilast, store

      data xilast, store /0.D0, 0.D0/

      ddum = dummy !avoid compiler warning      

      if(xi.eq.xilast) then
        Gbar_T = store
        return
      endif

*     term1 = 1.d0/(1.d0 + 0.25d0*xi)

*     if (xi .le. 1.d-1) then
*        Gbar_T = 1.d0/6.d0/pi*(-.5d0*term1 -
*    #        (1.d0 + 0.25d0*term1)*fjj_lim(xi) -
*    #        2.d0* (-1.d0/3.d0 + xi/15.d0 - xi**2/70.d0))
*     else
*        Gbar_T = 1.d0/6.d0/pi*(4.d0/xi - 0.5d0*term1
*    #        - (1.d0 + 2.d0/xi + 0.25d0*term1)*fjj(xi))
*     endif

      term1   = 1.d0/(1.d0 + 0.25d0*xi)

      Gbar_T  = 1.d0/6.d0/pi*(4.d0/xi - 0.5d0*term1
     #     - (1.d0 + 2.d0/xi + 0.25d0*term1)*fjj(xi))

      xilast = xi
      store  = Gbar_T

      return
      end

c     =================================
      double precision function fii(xi)
c     =================================

c     Equation (24) in PLB347 (1995) 143 - 151

      implicit none
      double precision pi, term1, term2, xi, di_log
      parameter (pi = 3.14159265359d0)

      term1 = dsqrt(xi)
      term2 = dsqrt(4.d0 + xi)
      fii = 4.d0/term1/term2*(-pi*pi/6.d0 
     #      - 0.5d0*(dlog((term2 + term1)/(term2 - term1)))**2
     #      + (dlog(0.5d0*(1.d0 - term1/term2)))**2 
     #      + 2.d0*di_log(0.5d0*(1.d0 - term1/term2)))

      return
      end

c     =================================
      double precision function fjj(xi)
c     =================================

c     Equation (23) in PLB347 (1995) 143 - 151

      implicit none
      double precision pi, xi, term1, term2
      parameter (pi = 3.14159265359d0)

      term1 = dsqrt(xi)
      term2 = dsqrt(4.d0 + xi)
      fjj = 4.d0/term1/term2*dlog((term2 + term1)/(term2 - term1))

      return
      end

c     =====================================
      double precision function fii_lim(xi)
c     =====================================

c     this gives fii(xi) in the limit that xi -> 0 up to xi**2

      implicit none
      double precision xi

      fii_lim = xi/3.d0 - xi**2/10.d0

      return
      end

c     =====================================
      double precision function fjj_lim(xi)
c     =====================================

c     this gives fjj(xi) in the limit that xi -> 0 up to xi**2

      implicit none
      double precision xi

      fjj_lim = 2.d0 - xi/3.d0 + xi**2/15.d0

      return
      end

c     ===================================
      double precision function di_log(x)
c     ===================================

c     Equation (25) in PLB347 (1995) 143 - 151

      implicit double precision  (a-z)
      dimension b(8)
      integer ncall
      data ncall/0/,pi6/1.644934066848226d+00/,een,vier/1.d+00,.25d+00/
 
      u     = 1
      ncall = 0
      if(ncall.eq.0)go to 2
1     if(x.lt.0)go to 3
      if(x.gt.0.5)go to 4
      z=-dlog(1.-x)
7     z2=z*z
      di_log=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
      if(x.gt.een)di_log=-di_log-.5*u*u+2.*pi6
      return
2     b(1)=een
      b(2)=een/36.
      b(3)=-een/3600.
      b(4)=een/211680.
      b(5)=-een/(30.*362880.d+00)
      b(6)=5./(66.*39916800.d+00)
      b(7)=-691./(2730.*39916800.d+00*156.)
      b(8)=een/(39916800.d+00*28080.)
      ncall=1
      go to 1
3     if(x.gt.-een)go to 5
      y=een/(een-x)
      z=-dlog(een-y)
      z2=z*z
      u=dlog(y)
      di_log=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier-u*(z+.5*u)-pi6
      return
4     if(x.ge.een)go to 10
      y=een-x
      z=-dlog(x)
6     u=dlog(y)
      z2=z*z
      di_log=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een-u)+z2*vier+pi6
      if(x.gt.een)di_log=-di_log-.5*z*z+pi6*2.
      return
5     y=een/(een-x)
      z=-dlog(y)
      z2=z*z
      di_log=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
      return
10    if(x.eq.een)go to 20
      xx=1./x
      if(x.gt.2.)go to 11
      z=dlog(x)
      y=1.-xx
      go to 6
11    u=dlog(x)
      z=-dlog(1.-xx)
      go to 7
20    di_log=pi6

      return
      end
