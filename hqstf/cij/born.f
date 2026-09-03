c This gives the Born coefficients 
c For QCD take tf = 1d0/2d0, for QED take  tf = 1d0.
c eta = (s - 4d0*m2)/4d0/m2, s is the gamma* gluon (gamma) CM Energy
c xi = Q^2/m2

c     =======================================
      double precision function C0_Lg(eta,xi)
c     =======================================

c     Longitudinal coefficient function: PL B347(1995)143 eq. (7).
c     This function is called born_l in the original code.

      implicit none
      double precision eta, xi, pi, tf
*     common/group/ca, cf, tf
      parameter(tf = 0.5d0)
      parameter(pi = 3.14159265359d0)

      C0_Lg  = 0.5d0*pi*tf*xi*(1.d0 + eta + 0.25d0*xi)**(-3.d0)*
     #         (2.d0*dsqrt(eta*(1.d0 + eta)) -
     #         dlog((dsqrt(1.d0 + eta) + dsqrt(eta))/
     #              (dsqrt(1.d0 + eta) - dsqrt(eta))))

      return
      end

c     =======================================
      double precision function C0_Tg(eta,xi)
c     =======================================

c     Transverse coefficient function: PL B347(1995)143 eq. (8).
c     This function is called born_t in the original code.

      implicit none
      double precision eta, xi, pi, tf
*     common/group/ca, cf, tf
      parameter(tf = 0.5d0)
      parameter(pi = 3.14159265359d0)

      C0_Tg  = 0.5d0*pi*tf*(1.d0 + eta + 0.25d0*xi)**(-3)*
     #         (-2.d0*((1.d0 + eta - 0.25d0*xi)**2 + eta + 1.d0)*
     #         dsqrt(eta/(1.d0 + eta)) + (2.d0*(1.d0 + eta)**2 +
     #         0.125d0*xi**2 + 2.d0*eta + 1.d0)*
     #         dlog((dsqrt(1.d0 + eta) + dsqrt(eta))/
     #              (dsqrt(1.d0 + eta) - dsqrt(eta))))

      return
      end
