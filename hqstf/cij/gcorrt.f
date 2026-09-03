
c     ========================================
      double precision function h1_ATg(eta,xi)
c     ========================================

c     Eq (9) in PLB347 (1995) 143 - 151 for the transverse piece.
c     MSbar scheme.
c     This routine is called subctca in the original code.
c     Called sctca in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1_ATg.inc'

      h1_ATg = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     ========================================
      double precision function h1_FTg(eta,xi)
c     ========================================

c     Eq (10) in PLB347 (1995) 143 - 151 for the transverse piece.
c     MSbar scheme.
c     This routine is called subctcf in the original code.
c     Called sctcf in updated code (03/06/96).

      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1_FTg.inc'

      h1_FTg = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end

c     ==========================================
      double precision function h1bar_Tg(eta,xi)
c     ==========================================

c     Eq (12) in PLB347 (1995) 143 - 151 for the transverse piece.
c     MSbar scheme.
c     This routine is called subctbar in the original code.
c     Called sctbar in updated code (03/06/96).
 
      implicit none
      integer neta, nxi
      parameter (neta = 73, nxi = 49)
      double precision calcpts(neta, nxi)
      double precision dlaeta(neta), dlaxi(nxi)
      double precision eta, xi
      double precision qchfun

      include 'h1bar_Tg.inc'

      h1bar_Tg = qchfun(eta,xi,calcpts,dlaeta,dlaxi)

      return
      end
