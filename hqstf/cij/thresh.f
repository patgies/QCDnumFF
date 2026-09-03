
c These are the functions that give the threshold dependence of the 
c coefficient functions with the appropriate factors.
c eta = (W^2 - 4d0*m2)/4d0/m2  where W is the CM energy of the 
c gamma* parton system. xi = mq2/m2 (Q^2/m2)

c     =========================================
      double precision function Efun_LF(eta,xi)
c     =========================================

c     Longitudinal CF group structure: eq (13) in PLB347 (1995) 143 - 151
c     This function is called threshf_l in the original code.

      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)

      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      Efun_LF = 1.d0/6.d0/pi*xi*term1**3*beta*beta*pi*pi/2.d0

      return
      end

c     =========================================
      double precision function Efun_TF(eta,xi)
c     =========================================

c     Transverse CF group structure: eq (14) in PLB347 (1995) 143 - 151
c     This function is called threshf_t in the original code.

      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)

      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      Efun_TF = 0.25d0/pi*term1*pi*pi/2.d0

      return
      end

c     =========================================
      double precision function Efun_LA(eta,xi)
c     =========================================

c     Longitudinal CA group structure: eq (15) in PLB347 (1995) 143 - 151
c     This function is called thresha_l in the original code.

      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)

      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      Efun_LA = 1.d0/6.d0/pi*xi*term1**3*beta**2*
     #     (beta*(dlog(8.d0*beta*beta))**2
     #     - 5.d0*beta*dlog(8.d0*beta*beta) - 0.25d0*pi*pi)

      return
      end

c     =========================================
      double precision function Efun_TA(eta,xi)
c     =========================================

c     Transverse CA group structure: eq (16) in PLB347 (1995) 143 - 151
c     This function is called thresha_t in the original code.
      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)

      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      Efun_TA = 0.25d0/pi*term1*(beta*(dlog(8.d0*beta*beta))**2
     #     - 5.d0*beta*dlog(8.d0*beta*beta) - 0.25d0*pi*pi)

      return
      end

c     =========================================
      double precision function Ebar_LA(eta,xi)
c     =========================================

c     Longitudinal CA group structure for the mass factorization piece: 
c     equation (17) in PLB347 (1995) 143 - 151
c     This function is called threshbar_l in the original code.

      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)

      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      Ebar_LA = 1.d0/6.d0/pi*xi*term1**3*beta**3*
     #     (-dlog(4.d0*beta*beta))

      return
      end

c     =========================================
      double precision function Ebar_TA(eta,xi)
c     =========================================

c     Transverse CA group structure for the mass factorization piece: 
c     equation (18) in PLB347 (1995) 143 - 151
c     This function is called threshbar_t in the original code.

      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)

      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      Ebar_TA = 0.25d0/pi*term1*beta*(-dlog(4.d0*beta*beta))

      return
      end
