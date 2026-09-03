
c     ========================================
      double precision function h1_ALg(eta,xi)
c     ========================================

c     Eq (9) in PLB347 (1995) 143 - 151 for the longitudinal piece.
c     MSbar scheme.
c     This routine is called subclca in the original code.
c     Called sclca in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1_ALg.inc'

      h1_ALg = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     ========================================
      double precision function h1_FLg(eta,xi)
c     ========================================

c     Eq (10) in PLB347 (1995) 143 - 151 for the longitudinal piece.
c     MSbar scheme.
c     This routine is called subclcf in the original code.
c     Called sclcf in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1_FLg.inc'

      h1_FLg = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     ==========================================
      double precision function h1bar_Lg(eta,xi)
c     ==========================================

c     Eq (12) in PLB347 (1995) 143 - 151 for the longitudinal piece.
c     MSbar scheme.
c     This routine is called subclbar in the original code.
c     Called sclbar in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1bar_Lg.inc'

      h1bar_Lg = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end
