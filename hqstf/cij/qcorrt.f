
c     ========================================
      double precision function h1_HTq(eta,xi)
c     ========================================

c     Eq (26) in PLB347 (1995) 143 - 151 for the transverse piece.
c     MSbar scheme.
c     This routine is called subcqht in the original code.
c     Called schqt in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1_HTq.inc'

      h1_HTq = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     ===========================================
      double precision function h1bar_HTq(eta,xi)
c     ===========================================

c     Eq (27) in PLB347 (1995) 143 - 151 for the transverse piece.
c     MSbar scheme.
c     This routine is called subcqhtbar in the original code.
c     Called sqtbar in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1bar_HTq.inc'

      h1bar_HTq = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     =========================================
      double precision function h1f_LTq(eta,xi)
c     =========================================

c     Eq (28) in PLB347 (1995) 143 - 151 for the transverse piece.
c     This also takes into account the additional mass factorizations
c     necessary from a low Q^2 photon coupling to the light quark.
c     MSbar scheme.
c     This routine is called subd1tqf in the original code.
c     Gives h1_LTq for Q2 < 1.5 GeV2 (use h1_LTq for Q2 > 1.5 GeV2)

      implicit none
      integer neta, nxi
      parameter (neta = 45, nxi = 15)
      double precision calcpts(neta, nxi)
      double precision aeta(neta), axi(nxi)
      double precision eta, xi, huge, small
      double precision t, u, y1, y2, y3, y4
      parameter (small = 1.d-8, huge = 1.d10)
      integer ieta, ixi

      include 'h1f_LTq.inc'
c
c  here we have to choose the array elements that will go into the 
c  interpolation.
      call locate(aeta, neta, eta, ieta)
      call locate(axi, nxi, xi, ixi)
      if (ieta .le. 1) ieta = 1
      if (ieta .gt. (neta - 1)) ieta = neta - 1
      if (ixi .le. 1) ixi = 1
      if (ixi .gt. (nxi - 1)) ixi = nxi - 1
      y1 = calcpts(ieta,ixi)
      y2 = calcpts(ieta+1,ixi)
      y3 = calcpts(ieta+1,ixi+1)
      y4 = calcpts(ieta,ixi+1)
c interpolating between the points
      t = (eta - aeta(ieta))/(aeta(ieta + 1) - aeta(ieta))
      u = (xi - axi(ixi))/(axi(ixi + 1) - axi(ixi))
      h1f_LTq = (1.d0 - t)*(1.d0 - u)*y1 + t*(1.d0 - u)*y2 +
     #           t*u*y3 + (1.d0 - t)*u*y4

      return
      end

c     =========================================
      double precision function h1_LTq(eta,xi)
c     =========================================

c     Eq (28) in PLB347 (1995) 143 - 151 for the transverse piece.
c     MSbar scheme.
c     This routine is called subd1tq in the original code.
c     Gives h1_LTq for Q2 > 1.5 GeV2 (use h1f_LTq for Q2 < 1.5 GeV2)
c     Called sclqt in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1_LTq.inc'

      h1_LTq = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     ===========================================
      double precision function h1bar_LTq(eta,xi)
c     ===========================================

c     Eq (29) in PLB347 (1995) 143 - 151 only necessary for the 
c     transverse piece.
c     MSbar scheme.
c     This routine is called subd1bar in the original code.
c     Gives h1bar_LTq for Q2 < 1.5 GeV2 ( = 0 for  Q2 > 1.5 GeV2)
      
      implicit none
      integer neta, nxi
      parameter (neta = 45, nxi = 15)
      double precision calcpts(neta, nxi)
      double precision aeta(neta), axi(nxi)
      double precision t, u, y1, y2, y3, y4
      double precision eta, xi, huge, small
      parameter (small = 1.d-8, huge = 1.d10)
      integer ieta, ixi

      include 'h1bar_LTq.inc'
c
c  here we have to choose the array elements that will go into the 
c  interpolation.
      call locate(aeta, neta, eta, ieta)
      call locate(axi, nxi, xi, ixi)
      if (ieta .le. 1) ieta = 1
      if (ieta .gt. (neta - 1)) ieta = neta - 1
      if (ixi .le. 1) ixi = 1
      if (ixi .gt. (nxi - 1)) ixi = nxi - 1
      y1 = calcpts(ieta,ixi)
      y2 = calcpts(ieta+1,ixi)
      y3 = calcpts(ieta+1,ixi+1)
      y4 = calcpts(ieta,ixi+1)
c interpolating between the points
      t = (eta - aeta(ieta))/(aeta(ieta + 1) - aeta(ieta))
      u = (xi - axi(ixi))/(axi(ixi + 1) - axi(ixi))
      h1bar_LTq = (1.d0 - t)*(1.d0 - u)*y1 + t*(1.d0 - u)*y2 +
     #             t*u*y3 + (1.d0 - t)*u*y4

      return
      end
