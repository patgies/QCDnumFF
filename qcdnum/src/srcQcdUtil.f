
C--   This is the file srcQcdUtil.f containing the qcdnum utility routines
C--   A utility routine has the characteristic that it does not 
C--   communicate via common blocks so that it can be used anywhere
C--
C--   subroutine sqcABmult(a,b,c,n)
C--   subroutine sqcNSmulti(a,nbnd,b,ci,i,ndim)
C--   subroutine sqcNSmult(a,nbnd,b,c,ndim)
C--   subroutine sqcNSeqsi(a,nbnd,x,i1,i2,b,n)
C--   subroutine sqcNSeqs(a,nbnd,x,b,n)
C--   subroutine sqcNSiter(a,nbnd,x,b,n,iter)
C--   subroutine sqcSGmulti(a,b,c,d,nbnd,e,f,gi,hi,i,ndim)
C--   subroutine sqcSGmult(a,b,c,d,nbnd,e,f,g,h,ndim)
C--   subroutine sqcSGeqsi(a,b,c,d,f,g,i1,i2,r,s,n)
C--   subroutine sqcSGeqs(a,b,c,d,f,g,r,s,n)
C--   subroutine sqcSGiter(a,b,c,d,f,g,r,s,n,iter)
C--   subroutine sqcLBeqs(a,m,nbnd,x,b,n)
C--   subroutine sqcUBeqs(a,m,nbnd,x,b,n)
C--
C--   subroutine sqcQHalf(iosp,acoef,sval,n)
C--   subroutine sqcLHalf(iosp,acoef,sval,n)
C--   subroutine sqcDHalf(iosp,acoef,epsi,n)
C--
C--   subroutine sqcOrtInv(amat,ainv,na,nf)

C--   subroutine sqcBrackit(func,x1,x2,ierr)
C--   double precision function dqcBiSect(func,x1,x2,epsi,ierr)

C--   subroutine sqcRange(v,n,vmi,vma,epsi,imi,ima,ierr)

C--   subroutine sqcGetLinCoef(ix1,iy1,ix2,iy2,a,b)
C--   integer function iqcIyFromIx(ix,a,b)

C--   integer function iqcLunFree

C     =============================      
      subroutine sqcABmult(a,b,c,n)
C     =============================

C--   Multiplies two lower triangular Toeplitz matrices A*B = C
C--   The matrix C is then also lower triangular Toeplitz
C--
C--   a    (in)   lower triangular matrix stored in vector (a1,...,an)
C--   b    (in)   lower triangular matrix stored in vector (b1,...,bn)
C--   c    (out)  lower triangular matrix stored in vector (c1,...,cn)
C--   n    (in)   dimension of a, b and c

      implicit double precision (a-h,o-z)
      
      dimension a(*),b(*),c(*)
      
      do i = 1,n
        ci = 0.D0
        do k = 1,i
          ci = ci + a(i-k+1)*b(k)
        enddo
        c(i) = ci
      enddo
      
      return
      end        

C     =========================================
      subroutine sqcNSmulti(a,nbnd,b,ci,i,ndim)
C     =========================================

C--   Returns the element ci of the matrix multiplication A*b = c
C--   A is a lower diagonal Toeplitz matrix
C--
C--       |a1             |   |b1|     |c1|
C--       |a2 a1          |   |b2|     |c2|
C--       |a3 a2 a1       | * |b3|  =  |c3|
C--       |   a3 a2 a1    |   |b4|     |c4|
C--       |      a3 a2 a1 |   |b5|     |c5|
C--
C--   a     (in)  band matrix stored in vector (a1,a2,a3,...)
C--   nbnd  (in)  number of elements in A (3 in the example above)
C--   b     (in)  vector to be multiplied with a
C--   ci    (out) element (i) of the result vector
C--   i     (in)  index of ci
C--   ndim  (in)  dimension of b and c (5 in the example above)

      implicit double precision (a-h,o-z)

      dimension a(*),b(*)    

      if(i.lt.1 .or. i.gt.ndim) stop 
     +          'sqcNSmulti: i out of range ---> STOP'
      ci  = 0.D0
      j1  = max(i+1-nbnd,1)
      ip1 = i+1
      do j = j1,i
        ci = ci + a(ip1-j)*b(j)
      enddo

      return
      end 

C     =====================================
      subroutine sqcNSmult(a,nbnd,b,c,ndim)
C     =====================================

C--   Returns the vector c of the matrix multiplication A*b = c
C--   A is a lower diagonal Toeplitz matrix
C--
C--       |a1             |   |b1|     |c1|
C--       |a2 a1          |   |b2|     |c2|
C--       |a3 a2 a1       | * |b3|  =  |c3|
C--       |   a3 a2 a1    |   |b4|     |c4|
C--       |      a3 a2 a1 |   |b5|     |c5|
C--
C--   a     (in)  band matrix stored in vector (a1,a2,a3,...)
C--   nbnd  (in)  number of elements in A (3 in the example above)
C--   b     (in)  vector to be multiplied with a
C--   c     (out) result vector
C--   ndim  (in)  dimension of b and c (5 in the example above)

      implicit double precision (a-h,o-z)

      dimension a(*),b(*),c(*)    

      do i = 1,ndim
        ci  = 0.D0
        j1  = max(i+1-nbnd,1)
        ip1 = i+1
        do j = j1,i
          ci = ci + a(ip1-j)*b(j)
        enddo
        c(i) = ci
      enddo

      return
      end

C     ========================================
      subroutine sqcNSeqsi(a,nbnd,x,i1,i2,b,n)
C     ======================================== 

C--   Not used at present

C--   Fill x(i1)--x(i2) of the solution x of the matrix equation Ax = b
C--   The elements x(1)--x(i1-1) must have been filled before 
C--   A is a lower diagonal Toeplitz matrix
C--
C--                                       | a1  0  0  0  0 |
C--                                       | a2 a1  0  0  0 |
C--   e.g. for 5 x 5 and bandwidth 3, A = | a3 a2 a1  0  0 |
C--                                       |  0 a3 a2 a1  0 |
C--                                       |  0  0 a3 a2 a1 |
C--
C--   a(n)   (in)    input vector with matrixelements
C--   nbnd   (in)    bandwidth (3 in the example above)
C--   x(n)   (inout) output solution vector x(i1)--x(i2)
C--   i1,i2  (in)    range of x elelemts to be filled
C--   b(n)   (in)    input right-hand side vector
C--   n      (in)    n dimension of a, x and b (5 in the example above)

      implicit double precision (a-h,o-z)

      dimension a(*), x(*), b(*)

*mb   n is a dummy variable
      idum = n

      if(i1.eq.1) then
        x(1) = b(1)/a(1)
      endif

      do i = min(i1,2),i2
        sum = 0.D0
        ip1 = i+1
        j1  = max(i+1-nbnd,1)
        do j = j1,i-1
          sum = sum + x(j)*a(ip1-j)
        enddo
        x(i) = (b(i)-sum)/a(1)
      enddo

      return
      end

C     =================================
      subroutine sqcNSeqs(a,nbnd,x,b,n)
C     ================================= 

C--   Solves the n x n matrix equation Ax = b where 
C--   A is a lower diagonal Toeplitz matrix
C--
C--                                       | a1             |
C--                                       | a2 a1          |
C--   e.g. for 5 x 5 and bandwidth 3, A = | a3 a2 a1       |
C--                                       |    a3 a2 a1    |
C--                                       |       a3 a2 a1 |
C--
C--   a(n)   (in)  input vector with matrixelements
C--   nbnd   (in)  bandwidth (3 in the example above)
C--   x(n)   (out) output solution vector
C--   b(n)   (in)  input right-hand side vector
C--   n      (in)  n dimension of a, x and b (5 in the example above)

      implicit double precision (a-h,o-z)

      dimension a(*), x(*), b(*)

      if(a(1).eq.0.D0) stop 'sqcNSeqs: singular matrix ---> STOP'

      x(1) = b(1)/a(1)

      do i = 2,n
        sum = 0.D0
        ip1 = i+1
        j1  = max(i+1-nbnd,1)
        do j = j1,i-1
          sum = sum + x(j)*a(ip1-j)
        enddo
        x(i) = (b(i)-sum)/a(1)
      enddo

      return
      end

C     =======================================
      subroutine sqcNSiter(a,nbnd,x,b,n,iter)
C     ======================================= 

C--   Not used at present
C--   iter = 0: 1-pass solution
C--   iter # 0: 2-pass solution

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension a(*), x(*), b(*)
      dimension db(mxx0),dx(mxx0)

      call sqcNSeqs(a,nbnd,x,b,n)
      if(iter.eq.0) return
      call sqcNSmult(a,nbnd,x,db,n)
      do i = 1,n
        db(i) = db(i)-b(i)
      enddo
      call sqcNSeqs(a,nbnd,dx,db,n)
      do i = 1,n
        x(i) = x(i)-dx(i)
      enddo

      return
      end

C--   ====================================================
      subroutine sqcSGmulti(a,b,c,d,nbnd,e,f,gi,hi,i,ndim)
C--   ====================================================

C--   Return the elements g(i) and h(i) of the matrix multiplication
C--
C--             | a b | * |e| = |g| 
C--             | c d |   |f|   |h|
C--
C--   a,b,c,d are lower diagonal ndim*ndim Toeplitz matrices 
C--
C--   |a1             |
C--   |a2 a1          |
C--   |a3 a2 a1       |  with, in this example, nbnd = 3 and ndim = 5 
C--   |   a3 a2 a1    |
C--   |      a3 a2 a1 |
C--
C--   a--d  (in)  band matrix stored in vector (a1,a2,a3,...)
C--   nbnd  (in)  number of elements in A (3 in the example above)
C--   e,f   (in)  vector to be multiplied with a
C--   gi,hi (out) elements g(i) and h(i) result vector
C--   i     (in)  index of g(i) and h(i)
C--   ndim  (in)  dimension (5 in the example above)

      implicit double precision (a-h,o-z)

      dimension a(*),b(*),c(*),d(*),e(*),f(*)    

      if(i.lt.1 .or. i.gt.ndim) stop 
     +              'sqcSGmulti: i out of range ---> STOP'
      gi = 0.D0
      hi = 0.D0
      j1 = max(i+1-nbnd,1)
      do j = j1,i
        gi = gi + a(i-j+1)*e(j) + b(i-j+1)*f(j)
        hi = hi + c(i-j+1)*e(j) + d(i-j+1)*f(j)
      enddo

      return
      end

C--   ===============================================
      subroutine sqcSGmult(a,b,c,d,nbnd,e,f,g,h,ndim)
C--   ===============================================

C--   Multiply  | a b | * |e| = |g| 
C--             | c d |   |f|   |h|
C--
C--   a,b,c,d are lower diagonal Toeplitz
C--
C--   |a1             |
C--   |a2 a1          |
C--   |a3 a2 a1       |  with, in this example, nbnd = 3 and ndim = 5 
C--   |   a3 a2 a1    |
C--   |      a3 a2 a1 |
C--
C--   a--d  (in)  band matrix stored in vector (a1,a2,a3,...)
C--   nbnd  (in)  number of elements in A (3 in the example above)
C--   e,f   (in)  vector to be multiplied with a
C--   g,h   (out) result vector
C--   ndim  (in)  dimension (5 in the example above)

      implicit double precision (a-h,o-z)

      dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)    

      do i = 1,ndim
        gi = 0.D0
        hi = 0.D0
        j1 = max(i+1-nbnd,1)
        do j = j1,i
          gi = gi + a(i-j+1)*e(j) + b(i-j+1)*f(j)
          hi = hi + c(i-j+1)*e(j) + d(i-j+1)*f(j)
        enddo
        g(i) = gi
        h(i) = hi
      enddo

      return
      end

C     =============================================
      subroutine sqcSGeqsi(a,b,c,d,f,g,i1,i2,r,s,n)
C     =============================================

C--   Not used at present

C--   Fill the elements i1--i2 of the solution vectors f and g of 
C--   the coupled singlet-gluon equations
C--
C--            |ff fg| |f|    |a b| |f|   |r|
C--            |gf gg| |g| =  |c d| |g| = |s|
C--
C--   f(1)--f(i1-1) and g(1)--g(i1-1) must have been filled before
C--
C--   ff, fg, gf and gg are n x n lower triangle Toeplitz matrices
C--
C--                                          | a1  0  0  0  0 |
C--                                          | a2 a1  0  0  0 |
C--   A(i,j<=i) = a(i-j+1)  e.g. for 5 x 5:  | a3 a2 a1  0  0 |
C--                                          | a4 a3 a2 a1  0 |
C--                                          | a5 a4 a3 a2 a1 |
C--
C--   a(n),b(n),c(n),d(n) (in)   input vectors with matrixelements
C--   f(n),g(n)           (out)  output solution vectors
C--   i1,i2               (in)   range of elements filled in f and g 
C--   r(n),s(n)           (in)   input right-hand side vectors
C--   n                   (in)   dimension of the matrices

      implicit double precision (a-h,o-z)

      dimension a(*), b(*), c(*), d(*)
      dimension f(*), g(*), r(*), s(*)

*mb   n is a dummy variable
      idum = n

      a1     = a(1)
      b1     = b(1)
      c1     = c(1)
      d1     = d(1)
      det    = (a1*d1-b1*c1)
      if(det.eq.0.D0) stop 'sqcSGeqs: singular matrix ---> STOP'
      detinv = 1.D0/det

      if(i1.eq.1) then
        f(1) = detinv*(r(1)*d1-s(1)*b1)
        g(1) = detinv*(s(1)*a1-r(1)*c1)
      endif

      do i = min(i1,2),i2
        rsum = r(i)
        ssum = s(i)
        ip1  = i+1
        do j = 1,i-1
          rsum = rsum - f(j)*a(ip1-j) - g(j)*b(ip1-j)
          ssum = ssum - f(j)*c(ip1-j) - g(j)*d(ip1-j)
        enddo
        f(i) = detinv*(rsum*d1-ssum*b1)
        g(i) = detinv*(ssum*a1-rsum*c1)
      enddo

      return
      end

C     ===============================================
      subroutine sqcQSGeqs(aa,bb,cc,dd,ff,gg,rr,ss,n)
C     ===============================================

C--   Solves coupled singlet-gluon equations
C--
C--            |ff fg| |f|    |a b| |f|   |r|
C--            |gf gg| |g| =  |c d| |g| = |s|
C--
C--   where ff, fg, gf and gg are n x n lower triangle Toeplitz matrices
C--
C--                                          | a1  0  0  0  0 |
C--                                          | a2 a1  0  0  0 |
C--   A(i,j<=i) = a(i-j+1)  e.g. for 5 x 5:  | a3 a2 a1  0  0 |
C--                                          | a4 a3 a2 a1  0 |
C--                                          | a5 a4 a3 a2 a1 |
C--
C--   a(n),b(n),c(n),d(n) (in)   input vectors with matrixelements
C--   f(n),g(n)           (out)  output solution vectors
C--   r(n),s(n)           (in)   input right-hand side vectors
C--   n                   (in)   dimension of the matrices

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension aa(*), bb(*), cc(*), dd(*)
      dimension ff(*), gg(*), rr(*), ss(*)

      double precision a(mxx0), b(mxx0), c(mxx0), d(mxx0)
      double precision f(mxx0), g(mxx0), r(mxx0), s(mxx0)
      double precision a1, b1, c1, d1, det, detinv, rsum, ssum
      double precision delr, dels

      do i = 1,n
        a(i) = aa(i)
        b(i) = bb(i)
        c(i) = cc(i)
        d(i) = dd(i)
      enddo
      r(1)   = rr(1)
      s(1)   = ss(1)
      do i = 2,n
        r(i) = rr(i)-rr(i-1)
        s(i) = ss(i)-ss(i-1)
      enddo

      a1     = a(1)
      b1     = b(1)
      c1     = c(1)
      d1     = d(1)
      det    = (a1*d1-b1*c1)
      if(det.eq.0.D0) stop 'sqcQSGeqs: singular matrix ---> STOP'
      detinv = 1.D0/det

      f(1) = detinv*(r(1)*d1-s(1)*b1)
      g(1) = detinv*(s(1)*a1-r(1)*c1)
      delr = r(1) - a1*f(1) - b1*g(1)
      dels = s(1) - c1*f(1) + d1*g(1)

      do i = 2,n
        rsum = r(i)
        ssum = s(i)
        ip1  = i+1
        do j = 1,i-1
          rsum = rsum - f(j)*a(ip1-j) - g(j)*b(ip1-j)
          ssum = ssum - f(j)*c(ip1-j) - g(j)*d(ip1-j)
        enddo
        write(6,'(I3,2E15.5)') i, rsum, ssum
        f(i) = detinv*(rsum*d1-ssum*b1)
        g(i) = detinv*(ssum*a1-rsum*c1)
      enddo

      ff(1) = dble(f(1))
      gg(1) = dble(g(1))
      do i = 2,n
        ff(i) = dble(f(i))+ff(i-1)
        gg(i) = dble(g(i))+gg(i-1)
      enddo

      return
      end

C     ======================================
      subroutine sqcSGeqs(a,b,c,d,f,g,r,s,n)
C     ======================================

C--   Solves coupled singlet-gluon equations
C--
C--            |ff fg| |f|    |a b| |f|   |r|
C--            |gf gg| |g| =  |c d| |g| = |s|
C--
C--   where ff, fg, gf and gg are n x n lower triangle Toeplitz matrices
C--
C--                                          | a1  0  0  0  0 |
C--                                          | a2 a1  0  0  0 |
C--   A(i,j<=i) = a(i-j+1)  e.g. for 5 x 5:  | a3 a2 a1  0  0 |
C--                                          | a4 a3 a2 a1  0 |
C--                                          | a5 a4 a3 a2 a1 |
C--
C--   a(n),b(n),c(n),d(n) (in)   input vectors with matrixelements
C--   f(n),g(n)           (out)  output solution vectors
C--   r(n),s(n)           (in)   input right-hand side vectors
C--   n                   (in)   dimension of the matrices

      implicit double precision (a-h,o-z)

      dimension a(*), b(*), c(*), d(*)
      dimension f(*), g(*), r(*), s(*)

      a1     = a(1)
      b1     = b(1)
      c1     = c(1)
      d1     = d(1)
      det    = (a1*d1-b1*c1)
      if(det.eq.0.D0) stop 'sqcSGeqs: singular matrix ---> STOP'
      detinv = 1.D0/det

      f(1) = detinv*(r(1)*d1-s(1)*b1)
      g(1) = detinv*(s(1)*a1-r(1)*c1)

      do i = 2,n
        rsum = r(i)
        ssum = s(i)
        ip1  = i+1
        do j = 1,i-1
          rsum = rsum - f(j)*a(ip1-j) - g(j)*b(ip1-j)
          ssum = ssum - f(j)*c(ip1-j) - g(j)*d(ip1-j)
        enddo
        f(i) = detinv*(rsum*d1-ssum*b1)
        g(i) = detinv*(ssum*a1-rsum*c1)
      enddo

      return
      end

C     ============================================
      subroutine sqcSGiter(a,b,c,d,f,g,r,s,n,iter)
C     ============================================ 

C--   Not used at present
C--   iter = 0: 1-pass solution
C--   iter # 0: 2-pass solution

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension a(*), b(*), c(*), d(*)
      dimension f(*), g(*), r(*), s(*)
      dimension df(mxx0),dg(mxx0),dr(mxx0),ds(mxx0)

      call sqcSGeqs(a,b,c,d,f,g,r,s,n)
      if(iter.eq.0) return
      call sqcSGmult(a,b,c,d,n,f,g,dr,ds,n)
      do i = 1,n
        dr(i) = dr(i)-r(i)
        ds(i) = ds(i)-s(i)
      enddo
      call sqcSGeqs(a,b,c,d,df,dg,dr,ds,n)
      do i = 1,n
        f(i) = f(i)-df(i)
        g(i) = g(i)-dg(i)
      enddo

      return
      end

C     ===================================
      subroutine sqcLBeqs(a,m,nbnd,x,b,n)
C     =================================== 

C--   Solves the n x n matrix equation Ax = b where A is a lower diagonal
C--   band matrix. See s/r sqcNSeqs when A has equal elements along the
C--   diagonals (Toeplitz matrix).
C--
C--                                  | a11                 |
C--                                  | a21 a22             |
C--   For 5 x 5 and bandwidth 3, A = | a31 a32 a33         |
C--                                  |     a42 a43 a44     |
C--                                  |         a53 a54 a55 |
C--
C--   a    (in)  input 2-dim matrix with n*n submatrix filled as example above
C--   m    (in)  dimension a(m,m) declared in the calling routine (m >= n)
C--   nbnd (in)  bandwidth (3 in the example above, nbnd <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of a, x and b (5 in the example above)

      implicit double precision (a-h,o-z)

      dimension a(m,m), x(*), b(*)

      x(1) = b(1)/a(1,1)

      do i = 2,n
        sum = 0.D0
        j1  = max(i+1-nbnd,1)
        do j = j1,i-1
          sum = sum + x(j)*a(i,j)
        enddo
        x(i) = (b(i)-sum)/a(i,i)
      enddo

      return
      end

C     ===================================
      subroutine sqcUBeqs(a,m,nbnd,x,b,n)
C     =================================== 

C--   Solves the n x n matrix equation Ax = b where A is a upper diagonal
C--   band matrix.
C--
C--                                  | a11 a12 a13         |
C--                                  |     a22 a23 a24     |
C--   For 5 x 5 and bandwidth 3, A = |         a33 a34 a35 |
C--                                  |             a44 a45 |
C--                                  |                 a55 |
C--
C--   a    (in)  input 2-dim matrix with n*n submatrix filled as example above
C--   m    (in)  dimension a(m,m) declared in the calling routine (m >= n)
C--   nbnd (in)  bandwidth (3 in the example above, nbnd <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of a, x and b (5 in the example above)

      implicit double precision (a-h,o-z)

      dimension a(m,m), x(*), b(*)

      x(n) = b(n)/a(n,n)

      do i = n-1,1,-1
        sum = 0.D0
        j2  = min(i-1+nbnd,n)
        do j = i+1,j2
          sum = sum + x(j)*a(i,j)
        enddo
        x(i) = (b(i)-sum)/a(i,i)
      enddo

      return
      end
      
C     ======================================
      subroutine sqcQHalf(iosp,acoef,sval,n)
C     ======================================

C--   Calculate quad spline interpolation at midpoints
C--
C--   iosp      (in)  :  spline order (should be 3 = quad)
C--   acoef(n)  (in)  :  vector of quad (!) spline coefficients
C--   sval(n)   (out) :  vector of interpolated results
C--   n         (in)  :  number of points in acoef and sval

      implicit double precision (a-h,o-z)
      
      save dmat
      dimension dmat(3), acoef(*), sval(*)
      
      logical first
      save first
      data first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,dmat)
      
      if(iosp.ne.3) stop 'sqcQHalf : not quad interpolation'
      
      if(first) then
        dmat(1) = dqcBsplyy(2,1,0.5D0)
        dmat(2) = dqcBsplyy(2,1,1.5D0)
        dmat(3) = dmat(1)
        first   = .false.
      endif
      
      call sqcNSmult(dmat,3,acoef,sval,n)  !sval = dmat*acoef
      
      return
      end
      
C     ======================================
      subroutine sqcLHalf(iosp,acoef,sval,n)
C     ======================================

C--   Calculate lin interpolation of a quad spline at midpoints
C--
C--   iosp      (in)  :  spline order (should be 3 = quad)
C--   acoef(n)  (in)  :  vector of quad (!) spline coefficients
C--   sval(n)   (out) :  vector of interpolated results
C--   n         (in)  :  number of points in acoef and sval

      implicit double precision (a-h,o-z)
      
      save emat
      dimension emat(3), acoef(*), sval(*)
      
      logical first
      save first
      data first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,emat)
      
      if(iosp.ne.3) stop 'sqcQHalf : not quad interpolation'
      
      if(first) then
        emat(1) = dqcBsplyy(2,1,1.0D0)/2.D0
        emat(2) = dqcBsplyy(2,1,1.0D0)
        emat(3) = emat(1)
        first   = .false.
      endif
      
      call sqcNSmult(emat,3,acoef,sval,n)  !sval = emat*acoef
      
      return
      end
      
C     ======================================
      subroutine sqcDHalf(iosp,acoef,epsi,n)
C     ======================================

C--   Calculate quad minus lin interpolation of a quad spline at midpoints
C--
C--   iosp      (in)  :  spline order (should be 3 = quad)
C--   acoef(n)  (in)  :  vector of quad (!) spline coefficients
C--   epsi(n)   (out) :  vector of deviations
C--   n         (in)  :  number of points in acoef and epsi

      implicit double precision (a-h,o-z)
      
      save dmat, emat, dmine
      dimension dmat(3), emat(3), dmine(3), acoef(*), epsi(*)
      
      logical first
      save first
      data first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,dmat,emat,dmine)
      
      if(iosp.ne.3) stop 'sqcQHalf : not quad interpolation'
      
      if(first) then
        dmat(1) = dqcBsplyy(2,1,0.5D0)
        dmat(2) = dqcBsplyy(2,1,1.5D0)
        dmat(3) = dmat(1)
        emat(1) = dqcBsplyy(2,1,1.0D0)/2.D0
        emat(2) = dqcBsplyy(2,1,1.0D0)
        emat(3) = emat(1)
        do i = 1,3
        dmine(i) = dmat(i) - emat(i)
        enddo
        first   = .false.
      endif
      
      call sqcNSmult(dmine,3,acoef,epsi,n)  !epsi = dmine*acoef
      
      return
      end


C     =====================================
      subroutine sqcOrtInv(amat,ainv,na,nf)
C     =====================================

C--   Invert an row-wise orthogonal (not necessarily orthonormal) matrix A. 
C--   Let S be the diagonal scaling matrix s_i = sqrt{1/(Sum_j^n A_ij^2}, 
C--   then the matrix B = SA is an ortonormal matrix (rotation) for which
C--   the inverse is the transpose. Scaling the transpose again we find for
C--   the inverse of A: Ainv = A^{transpose} S^2.
C--
C--   amat   input  matrix dimensioned amat(na,na) in the calling routine
C--   ainv   output matrix dimensioned ainv(na,na) in the calling routine
C--   na     first and second dimension of amat and ainv
C--   nf     dimension of the nf*nf submatrix of amat to be inverted  

      implicit double precision (a-h,o-z)

      dimension amat(na,na),ainv(na,na),scale(na)

C--   Protect
      if(nf.le.0 .or. na.lt.nf) 
     &stop 'sqcOrtInv: wrong input dimensions --> STOP'
C--   Initialize
      do i = 1,na
        scale(i) = 0.D0
        do j = 1,na
          ainv(i,j) = 0.D0
        enddo
      enddo
C--   Scaling matrix squared
      do i = 1,nf
        sc = 0.D0
        do j = 1,nf
          sc  = sc + amat(i,j)*amat(i,j)
        enddo
        if(sc.le.0.D0) stop 'sqcOrtInv: singular matrix --> STOP'
        scale(i) = 1.D0/sc
      enddo
C--   Ainv = Atranspose*scale^2
      do i = 1,nf
        do j = 1,nf
          ainv(i,j) = amat(j,i)*scale(j)
        enddo
      enddo

      return
      end

C     ======================================
      subroutine sqcBrackit(func,x1,x2,ierr)
C     ======================================

C--   Find the range [x1,x2] in which func(x) = 0 [from Numerical Recipes]
C--
C--   func(x) :  function to be declared external in the calling routine
C--   x1      : (inout) initial guess on entry, final lower limit on exit
C--   x2      : (inout) initial guess on entry, final upper limit on exit
C--   ierr    : (out)   0 = OK; 1 = cannot find range where func(x) = 0 

      implicit double precision (a-h,o-z)

      external func

      parameter (factor = 1.6, ntry = 50)

      if(x1.eq.x2) stop 'sqcBrackit: x1 = x2 not allowed ---> STOP'

      f1   = func(x1)
      f2   = func(x2)
      ierr = 0
      do j = 1,ntry
        if(f1*f2.lt.0.D0) return
        if(abs(f1).lt.abs(f2)) then
          x1 = x1 + factor*(x1-x2)
          f1 = func(x1)
        else
          x2 = x2 + factor*(x2-x1)
          f2 = func(x2)
        endif
      enddo
      ierr = 1

      return
      end

C     =========================================================
      double precision function dqcBiSect(func,x1,x2,epsi,ierr)
C     =========================================================

C--   Find root of func(x) within range [x1,x2] [from Numerical Recipes]
C--
C--   func(x) :  function to be declared external in the calling routine
C--   x1      : (in)  lower limit of search range
C--   x2      : (in)  upper limit of search range
C--   epsi    : (in)  iterate till root has an accuracy +- epsi
C--   ierr    : (out) 0 = OK; 1 = accuracy not reached

      implicit double precision (a-h,o-z)

      external func

      parameter (jmax = 40)

      fmid = func(x2)
      f    = func(x1)
      if(f*fmid.ge.0.D0) stop 
     +          'dqcBiSect: [x1,x2] does not bracket root ---> STOP'
      if(f.lt.0.D0) then
        dqcbisect = x1
        dx        = x2-x1
      else
        dqcbisect = x2
        dx        = x1-x2
      endif
      ierr = 0
      do j = 1,jmax
        dx   = 0.5D0*dx
        xmid = dqcbisect + dx
        fmid = func(xmid)
        if(fmid.le.0.D0) dqcbisect = xmid
        if(abs(dx).lt.epsi .or. fmid.eq.0.D0) return
      enddo
C--   Land here if root not found after jmax bisections
      dqcbisect = 1.D11
      ierr      = 1
      return
      end
      
C     ==================================================      
      subroutine sqcRange(v,n,vmi,vma,epsi,imi,ima,ierr)
C     ==================================================

C--   Find elements of an array inbetween cuts vmi and vma
C--
C--   v      (in)   array in ascending order
C--   n      (in)   dimension of v
C--   vmi    (in)   lower cut
C--   vma    (in)   upper cut
C--   epsi   (in)   x == y when |x-y| < epsi
C--   imi    (out)  index of first element .ge. vmi (0 if none)
C--   ima    (out)  index of last  element .le. vma (0 if none)
C--   ierr   (out)  0 all OK
C--                 1 no elements .ge. vmi and .le. vma (imi = ima = 0)
C--                 2 array not in ascending order      (imi = ima = 0)

      implicit double precision (a-h,o-z)
      
      logical lmb_ge, lmb_le
      
      dimension v(*)
      
      if(n.le.0) stop 'sqcRange: n .le. 0  ---> STOP'
      
      if(lmb_ge(vmi,vma,epsi))
     +           stop 'sqcRange: vmi .ge. vma ---> STOP'
  
      imi  = 0
      ima  = 0
      ierr = 0
      if(lmb_ge(v(1),vmi,epsi))     imi = 1
      if(lmb_le(v(n),vma,epsi))     ima = n
      do i = 2,n
        if(lmb_ge(v(i-1),v(i),epsi)) then
          imi  = 0
          ima  = 0
          ierr = 2
          return
        endif  
        if(imi.eq.0 .and. lmb_ge(v(i),vmi,epsi))     imi = i
        if(ima.eq.0 .and. lmb_le(v(n+1-i),vma,epsi)) ima = n+1-i
      enddo 
      
      if(imi.eq.0 .or. ima.eq.0) then
        imi  = 0
        ima  = 0
        ierr = 1
      endif  
      
      return
      end   

C     =============================================
      subroutine sqcGetLinCoef(ix1,iy1,ix2,iy2,a,b)
C     =============================================

C--   Calculate coefficients of the transformation iy = a*ix + b

      implicit double precision (a-h,o-z)

      xdif = ix1-ix2
      ydif = iy1-iy2
      xxyy = ix1*iy2-ix2*iy1
      a    = ydif/xdif
      b    = 2.D0*xxyy/xdif

      return
      end

C     ====================================
      integer function iqcIyFromIx(ix,a,b)
C     ====================================

C--   Calculate y = a*ix + b and round to nearest integer

      implicit double precision (a-h,o-z)

      aa = a  !avoid compiler warning
      bb = b  !avoid compiler warning
      jx = ix !avoid compiler warning
*mb      iqcIyFromIx = anint(a*ix+b)
      iqcIyFromIx = 0

      return
      end

C     =================================
      integer function iqcLunFree(lmin)
C     =================================

C--   Find first free lun in the range [max(lmin,10)-99]
C--   Returns zero if no free lun found

      logical lopened

      iqcLunFree = 0
      lopened    = .true.
      lun        = max(lmin,10)-1
      do while(lopened .eqv. .true. .and. lun .le. 99)
        lun = lun + 1
        inquire(unit=lun,opened=lopened)
      enddo
      if(lun.ne.100) iqcLunFree = lun

      return
      end


