
C-MB  ------------------------------------------------------------------
C-MB  Historical note: this is the only piece of original QCDNUM code
C-MB  left, written by Ouraou and Virchaud in 1988 ...
C-MB  ------------------------------------------------------------------

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGS(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'

      IDUM     = NF        !avoid compiler warning      
      dqcP0GGS = 1.D0 / ( 1.D0 - X )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGR(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'

      IDUM     = NF        !avoid compiler warning       
      dqcP0GGR = 6.D0 * X 
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGA(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'

      IDUM     = NF        !avoid compiler warning       
      ONEMX    = 1.D0 - X
      dqcP0GGA = 6.D0 * ( ONEMX/X + X*ONEMX )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GFA(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'
 
      IDUM     = NF        !avoid compiler warning 
      dqcP0GFA = 4.D0 * ( 1.D0 + (1.D0-X)*(1.D0-X) ) / ( 3.D0*X )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FGA(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'
 
      IDUM     = NF        !avoid compiler warning 
      dqcP0FGA = 0.5D0 * ( X*X + (1.D0-X)*(1.D0-X) )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FFS(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'
 
      IDUM     = NF        !avoid compiler warning 
      dqcP0FFS = 1.D0 / (1.D0-X)
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FFR(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'

      IDUM     = NF        !avoid compiler warning 
      dqcP0FFR = 4.D0 * ( 1.D0 + X*X ) / 3.D0
 
      RETURN
      END
