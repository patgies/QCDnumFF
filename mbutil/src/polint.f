C--
C--   This is the file polint.f with fast interpolation routines
C--
C--   subroutine smb_polwgt(x,xi,n,w)
C--   double precision function dmb_polin1(w,fi,n)
C--   double precision function dmb_polin2(wx,nx,wy,ny,fij,m)
C--
C--   subroutine smb_spline(n,x,y,b,c,d)
C--   double precision function dmb_seval(n,u,x,y,b,c,d)
C--
C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  MBUTIL interpolation routines from polint.f             */
CXXHDR    /************************************************************/
CXXHDR
CXXHDR    // Already defined in file wspace.f
CXXHDR    // inline int iaFtoC(int ia) { return ia-1; };
CXXHDR    // inline int iaCtoF(int ia) { return ia+1; };
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  MBUTIL interpolation routines from polint.f               */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  MBUTIL interpolation routines from polint.f               */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C==   ==================================================================
C==   Local polynomial interpolation ===================================
C==   ==================================================================

C     ===============================
      subroutine smb_polwgt(x,xi,n,w)
C     ===============================

C--   Setup weights for 1, 2, 3-point interpolation

C--   Author: M. Botje h24@nikhef.nl

C--   Bug fixed at 27-04-16

      implicit double precision (a-h,o-z)

      dimension xi(*), w(*)

      if(n.eq.3)     then
        w(1) = (xi(2)-x)/(xi(2)-xi(1))
        w(2) = 1.D0 - w(1)
        w(3) = (xi(3)-x)/(xi(3)-xi(2))
        w(4) = 1.D0 - w(3)
        w(5) = (xi(3)-x)/(xi(3)-xi(1))
        w(6) = 1.D0 - w(5)
      elseif(n.eq.2) then
        w(1) = (xi(2)-x)/(xi(2)-xi(1))
        w(2) = 1.D0 - w(1)
      elseif(n.eq.1) then
        w(1) = 1.D0
      else
        stop 'SMB_POLWGT: invalid interpolation order'
      endif

      return
      end

C     ============================================
      double precision function dmb_polin1(w,fi,n)
C     ============================================

C--   1-dim interpolation of order 1,2,3

C--   Author: M. Botje h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dimension w(*), fi(*)

      if(n.eq.3)     then
        g          = w(1)*fi(1) + w(2)*fi(2)
        h          = w(3)*fi(2) + w(4)*fi(3)
        dmb_polin1 = w(5)*g + w(6)*h
      elseif(n.eq.2) then
        dmb_polin1 = w(1)*fi(1) + w(2)*fi(2)
      elseif(n.eq.1) then
        dmb_polin1 = fi(1)
      else
        stop 'SMB_POLIN1: invalid interpolation order'
      endif

      return
      end

C     =======================================================
      double precision function dmb_polin2(wx,nx,wy,ny,fij,m)
C     =======================================================

C--   2-dim interpolation on an nx * ny mesh

C--   Author: M. Botje h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dimension wx(*), wy(*), fij(m,*), fi(3)

      if(nx.lt.1 .or. nx.gt.3)
     + stop 'SMB_POLIN2: invalid interpolation order in x'
      if(ny.lt.1 .or. ny.gt.3)
     + stop 'SMB_POLIN2: invalid interpolation order in y'

       do j = 1,ny
          fi(j) = dmb_polin1(wx,fij(1,j),nx)
       enddo
       dmb_polin2 = dmb_polin1(wy,fi,ny)

       return
       end

C==   ==================================================================
C==   Cubic spline interpolation =======================================
C==   ==================================================================

C-----------------------------------------------------------------------
CXXHDR    void smb_spline(int n, double *x, double *y, double *b, double *c, double *d);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_spline FC_FUNC(smb_spline,SMB_SPLINE)
CXXHFW    void fsmb_spline(int*, double*, double*, double*, double*, double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    void smb_spline(int n, double *x, double *y, double *b, double *c, double *d)
CXXWRP    {
CXXWRP      fsmb_spline(&n, x, y, b, c, d);
CXXWRP    }
C-----------------------------------------------------------------------

C     =======================================
      subroutine smb_spline(n, x, y, b, c, d)
C     =======================================

      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) stop 'SMB_SPLINE: need at least two node points'
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dmb_seval(int n, double u, double *x, double *y, double *b, double *c, double *d);
C-----------------------------------------------------------------------
CXXHFW  #define fdmb_seval FC_FUNC(dmb_seval,DMB_SEVAL)
CXXHFW    double fdmb_seval(int*, double*, double*, double*, double*, double*, double*);
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    double dmb_seval(int n, double u, double *x, double *y, double *b, double *c, double *d)
CXXWRP    {
CXXWRP      return fdmb_seval(&n, &u, x, y, b, c, d);
CXXWRP    }
C-----------------------------------------------------------------------

C     ========================================================
      double precision function dmb_seval(n, u, x, y, b, c, d)
C     ========================================================

      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
c
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
      integer i, j, k
      double precision dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      dmb_seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

      return
      end

