
C     ======================================
      DOUBLE PRECISION FUNCTION DEEJ3Q(X,NF)
C     ======================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C--   C3Q = C2Q - D3Q

      IDUM    = NF            !avoid compiler warning
      DEEJ3Q  = (4.D0/3.D0)*(1.+X)
 
      RETURN
      END

C     ======================================
      DOUBLE PRECISION FUNCTION CEEJLQ(X,NF)
C     ======================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      IDUM    = NF        !avoid compiler warning
      CEEJLQ  = (8.D0/3.D0)*X
 
      RETURN
      END

C     ======================================
      DOUBLE PRECISION FUNCTION CEEJLG(X,NF)
C     ======================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      CEEJLG = NF*4.*X*(1.-X)
 
      RETURN
      END

C     ======================================
      DOUBLE PRECISION FUNCTION CEEJ2Q(X,NF)
C     ======================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      IDUM   = NF          !avoind compiler warning
      C1MX   = 1.-X
      CEEJ2Q = 3. + (5.D0/3.D0)*X + ((4.D0/3.D0)*LOG(C1MX/X)-1.) * 
     &        (1.+X**2) / C1MX
 
      RETURN
      END

C     ======================================
      DOUBLE PRECISION FUNCTION CEEJ2G(X,NF)
C     ======================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      C1MX    = 1. - X
      CEEJ2G  = -.5 + 4.*X*C1MX + .5 * (X**2+C1MX**2) * LOG(C1MX/X)
      CEEJ2G  = 2.*NF*CEEJ2G
 
      RETURN
      END
