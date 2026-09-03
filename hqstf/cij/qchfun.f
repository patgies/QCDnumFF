
c     =============================================================
      double precision function qchfun(eta,xi,calcpts,dlaeta,dlaxi)
c     =============================================================

c     Interpolation routine used in most, but not all, heavy
c     quark coefficient functions. The interpolation tables
c     calcpts, dlaeta and dlaxi are passed as arguments.

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi, dleta, dlxi
      double precision pxi, peta, f(-1:1), delxi, deleta
      integer ieta, ixi

      dleta = dlog10(eta)
      dlxi = dlog10(xi)
      if (dlxi .le. dlaxi(1)) dlxi = dlaxi(1)
      if (dlxi .ge. dlaxi(nxi)) dlxi = dlaxi(nxi)
      if (dleta .ge. dlaeta(neta)) dleta = dlaeta(neta)
      if (dleta .le. dlaeta(1)) dleta = dlaeta(1)
      call locate(dlaeta,neta, dleta, ieta)
      call locate(dlaxi, nxi, dlxi, ixi)
c     interpolating between the appropriate points
      delxi = 1d0/6d0
      deleta = 1d0/6d0
c  lagrange 3-pt.
      if (ixi .le. 2) ixi = 2
      if (ixi .ge. 48) ixi = 48
      if (ieta .le. 2) ieta = 2
      if (ieta .ge. 72) ieta = 72
      pxi = (dlxi - dlaxi(ixi))/delxi
      f(-1) = pxi*(pxi-1d0)/2d0*calcpts(ieta-1,ixi-1) +
     #     (1d0 - pxi**2)*calcpts(ieta-1,ixi) +
     #     pxi*(pxi+1d0)/2d0*calcpts(ieta-1,ixi+1)
      f(0) = pxi*(pxi-1d0)/2d0*calcpts(ieta,ixi-1) +
     #     (1d0 - pxi**2)*calcpts(ieta,ixi) +
     #     pxi*(pxi+1d0)/2d0*calcpts(ieta,ixi+1)
      f(1) = pxi*(pxi-1d0)/2d0*calcpts(ieta+1,ixi-1) +
     #     (1d0 - pxi**2)*calcpts(ieta+1,ixi) +
     #     pxi*(pxi+1d0)/2d0*calcpts(ieta+1,ixi+1)
      peta = (dleta - dlaeta(ieta))/deleta
      qchfun = peta*(peta-1d0)/2d0*f(-1) +
     #     (1d0 - peta**2)*f(0) +
CMB  #     + peta*(peta+1d0)/2d0*f(1)
     #     peta*(peta+1d0)/2d0*f(1)
      return
      end
