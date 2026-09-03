
c     ===========================================
      double precision function h1_HLq(eta,xi)
c     ===========================================

c     Eq (26) in PLB347 (1995) 143 - 151 for the transverse piece.
c     MSbar scheme.
c     This routine is called subcqhl in the original code.
c     Called schql in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1_HLq.inc'

      h1_HLq = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     ===========================================
      double precision function h1bar_HLq(eta,xi)
c     ===========================================

c     Eq (27) in PLB347 (1995) 143 - 151 for the transverse piece.
c     MSbar scheme.
c     This routine is called subcqhlbar in the original code.
c     Called sqlbar in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1bar_HLq.inc'

      h1bar_HLq = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     =========================================
      double precision function h1f_LLq(eta,xi)
c     =========================================

c     Eq (28) in PLB347 (1995) 143 - 151 for the longitudinal piece
c     This also takes into account the additional mass factorizations
c     necessary from a low Q^2 photon coupling to the light quark.
c     MSbar scheme.
c     This routine is called subd1lqf in the original code.
c     Gives h1_LLq for Q2 < 1.5 GeV2 (use h1_LLq for Q2 > 1.5 GeV2).

      implicit none
      integer neta, nxi
      parameter (neta = 45, nxi = 15)
      double precision calcpts(neta, nxi)
      double precision aeta(neta), axi(nxi)
      double precision eta, xi, huge, small
      double precision t, u, y1, y2, y3, y4
      parameter (small = 1.d-8, huge = 1.d10)
      integer ieta, ixi

      include 'h1f_LLq.inc'
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
      h1f_LLq = (1.d0 - t)*(1.d0 - u)*y1 + t*(1.d0 - u)*y2 +
     #           t*u*y3 + (1.d0 - t)*u*y4

      return
      end

c     ========================================
      double precision function h1_LLq(eta,xi)
c     ========================================

c     Eq (28) in PLB347 (1995) 143 - 151 for the longitudinal piece.
c     MSbar scheme.
c     This routine is called subd1lq in the original code.
c     Gives h1_LLq for Q2 > 1.5 GeV2 (use h1f_LLq for Q2 < 1.5 GeV2).
c     Called sclql in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1_LLq.inc'

      h1_LLq = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end
