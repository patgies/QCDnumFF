
C--  This is the file srcTransform.f with pdf transformation routines

C--   subroutine sqcPDFtoEPM(tmatqf,pdf,epm,nf)

c--   subroutine sqcGetMatQF(tmatfq,tmatqf,ierr)
C--   subroutine sqcGetMatEQ(tmateq,nf)
C--   subroutine sqcGetMatQE(tmatqe,nf)

C--   subroutine sqcQQBtoEPM(qqb,epm,nf)
C--   subroutine sqcEPMtoQQB(qqb,epm,nf)
C--   subroutine sqcQQBtoQPM(qqb,qpm)
C--   subroutine sqcQPMtoQQB(qpm,qqb)
C--   subroutine sqcQPMtoEPM(qpm,epm,nf)
C--   subroutine sqcEPMtoQPM(epm,qpm,nf)
C--   subroutine sqcQtoE6(qval,eval,nf)
C--   subroutine sqcEtoQ6(eval,qval,nf)

C======================================================================
C==   Store transformation matrices for later use
C=======================================================================

C     =====================
      subroutine sqcPdfMatn
C     =====================

C--          1   2   3   4   5   6   7   8   9  10  11  12  13
C--   qqb = tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   epm = e0  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2- e3- e4- e5- e6-

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpdfs7.inc'

C--   Store qqb(i) = Sum_j tmatqe(i,j) epm(j)
      do nf = 3,6
        call sqcGetMatQE(tmatqen7(1,1,nf),nf)
      enddo

      return
      end

C=======================================================================
C==   Transform user input pdfs
C=======================================================================

C     =========================================
      subroutine sqcPDFtoEPM(tmatqf,pdf,epm,nf)
C     =========================================

C--   Transform 13 pdf(x,q2) --> epm(x,q2)
C--
C--   pdf =  g  f1  f2  f3  f4  f5  f6  f7  f8  f9 f10 f11 f12
C--          1   2   3   4   5   6   7   8   9  10  11  12  13
C--   qqb = tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   epm = e0  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2- e3- e4- e5- e6-
C--
C--   tmatqf(13,13)  (in) : matrix qqb(i) = Sum_j tmatqf(i,j) * pdf(j)
C--   pdf(13)        (in) : array with pdf_i values at some (x,q2)
C--   emm(13)       (out) : array with epm_i values at some (x,q2)
C--   nf             (in) : number of active flavours

      implicit double precision (a-h,o-z)

      dimension tmatqf(13,13), pdf(*), epm(*), qqb(13)

C--   Get qqb
      do i = 1,13
        qqb(i) = 0.D0
        do j = 1,13
          qqb(i) = qqb(i) + tmatqf(i,j)*pdf(j)
        enddo
      enddo
C--   Get epm
      call sqcQQBtoEPM(qqb,epm,nf)

      return
      end

C=======================================================================
C==   Transformation matrices  =========================================
C=======================================================================

C     ==========================================
      subroutine sqcGetMatQF(tmatfq,tmatqf,ierr)
C     ==========================================

C--   Invert tmatfq (in) --> tmatqf (out)
C--
C--   pdf =  g  f1  f2  f3  f4  f5  f6  f7  f8  f9 f10 f11 f12
C--          1   2   3   4   5   6   7   8   9  10  11  12  13
C--   qqb = tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--
C--   Note: pdf(i) = Sum_j tmatfq(i,j) * qqb(j)
C--         qqb(i) = Sum_j tmatqf(i,j) * pdf(j)
C--
C--   Out:  ierr !=0 : input pdfs not independent

      implicit double precision (a-h,o-z)

      dimension tmatfq(13,13), tmatqf(13,13), iwork(13)

C--   smb_dinv is destructive so first copy
      do i = 1,13
        do j = 1,13
          tmatqf(i,j) = tmatfq(i,j)
        enddo
      enddo

      call smb_dminv(13,tmatqf,13,iwork,ierr)

      return
      end

C     =================================
      subroutine sqcGetMatEQ(tmateq,nf)
C     =================================

C--   qqb = tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--          1   2   3   4   5   6   7   8   9  10  11  12  13
C--   epm = e0  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2- e3- e4- e5- e6-
C--
C--   epm(i) = Sum_j tmateq(i,j) * qqb(j)
C--
C--   NB: tmateq columns are images of vectors (1,0,0,...), (0,1,0,...)

      implicit double precision (a-h,o-z)

      dimension tmateq(13,13), qqb(13)

      call smb_VFill(qqb,13,0.D0)

      do i = 1,13
        qqb(i) = 1.D0
        call sqcQQBtoEPM(qqb,tmateq(1,i),nf)
        qqb(i) = 0.D0
      enddo

      return
      end

C     =================================
      subroutine sqcGetMatQE(tmatqe,nf)
C     =================================

C--   qqb = tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--          1   2   3   4   5   6   7   8   9  10  11  12  13
C--   epm = e0  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2- e3- e4- e5- e6-
C--
C--   qqb(i) = Sum_j tmatqe(i,j) * epm(j)
C--
C--   NB: tmatqe columns are images of vectors (1,0,0,...), (0,1,0,...)

      implicit double precision (a-h,o-z)

      dimension tmatqe(13,13), epm(13)

      call smb_VFill(epm,13,0.D0)

      do i = 1,13
        epm(i) = 1.D0
        call sqcEPMtoQQB(epm,tmatqe(1,i),nf)
        epm(i) = 0.D0
      enddo

      return
      end

C=======================================================================
C==   Transformation routines  =========================================
C=======================================================================

C     ==================================
      subroutine sqcQQBtoEPM(qqb,epm,nf)
C     ==================================

C--   Transform 13 qqb(x,q2) --> epm(x,q2)
C--
C--   qqb = tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--          1   2   3   4   5   6   7   8   9   10  11  12  13
C--   epm = e0  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2-  e3- e4- e5- e6-
C--
C--   qqb(13)   (in) : array with qqb_i values at some (x,q2)
C--   epm(13)  (out) : array with epm_i values at some (x,q2)
C--   nf        (in) : number of active flavours

      implicit double precision (a-h,o-z)

      dimension qqb(*), epm(*), qpm(13)

      call sqcQQBtoQPM(qqb,qpm)
      call sqcQPMtoEPM(qpm,epm,nf)

      return
      end

C     ==================================
      subroutine sqcEPMtoQQB(epm,qqb,nf)
C     ==================================

C--   Transform 13 epm(x,q2) --> qqb(x,q2)
C--
C--   epm = e0  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2-  e3- e4- e5- e6-
C--          1   2   3   4   5   6   7   8   9   10  11  12  13
C--   qqb = tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--
C--   epm(13)   (in) : array with epm_i values at some (x,q2)
C--   qqb(13)  (out) : array with qqb_i values at some (x,q2)
C--   nf        (in) : number of active flavours

      implicit double precision (a-h,o-z)

      dimension epm(*), qqb(*), qpm(13)

      call sqcEPMtoQPM(epm,qpm,nf)
      call sqcQPMtoQQB(qpm,qqb)

      return
      end

C     ===============================
      subroutine sqcQQBtoQPM(qqb,qpm)
C     ===============================

C--   Transform 13 qqb(x,q2) --> qpm(x,q2)
C--
C--   qqb = tb  bb  cb  sb  ub  db  g   d   u   s   c   b   t
C--          1  2   3   4   5   6   7   8   9   10  11  12  13
C--   qpm =  g  d+  u+  s+  c+  b+  t+  d-  u-  s-  c-  b-  t-
C--
C--   qqb(13)   (in) : array with qqb_i values at some (x,q2)
C--   qpm(13)  (out) : array with qpm_i values at some (x,q2)

      implicit double precision (a-h,o-z)

      dimension qqb(*), qpm(*)

      qpm(1) = qqb(7)

      do i = 1,6
        qpm(i+1) = qqb(7+i) + qqb(7-i)
        qpm(i+7) = qqb(7+i) - qqb(7-i)
      enddo

      return
      end

C     ===============================
      subroutine sqcQPMtoQQB(qpm,qqb)
C     ===============================

C--   Transform 13 qpm(x,q2) --> qqb(x,q2)
C--
C--   qpm =  g  d+  u+  s+  c+  b+  t+  d-  u-  s-  c-  b-  t-
C--          1  2   3   4   5   6   7   8   9   10  11  12  13
C--   qqb = tb  bb  cb  sb  ub  db  g   d   u   s   c   b   t
C--
C--   qpm(13)   (in) : array with qpm_i values at some (x,q2)
C--   qqb(13)  (out) : array with qqb_i values at some (x,q2)

      implicit double precision (a-h,o-z)

      dimension qpm(*), qqb(*)

      qqb(7) = qpm(1)

      do i = 1,6
        qqb(7+i) = 0.5D0 * ( qpm(i+1) + qpm(i+7) )
        qqb(7-i) = 0.5D0 * ( qpm(i+1) - qpm(i+7) )
      enddo

      return
      end

C     ==================================
      subroutine sqcQPMtoEPM(qpm,epm,nf)
C     ==================================

C--   Transform 13 qpm(x,q2) --> epm(x,q2)
C--
C--   qpm =  g  d+  u+  s+  c+  b+  t+  d-  u-  s-  c-  b-  t-
C--          1  2   3   4   5   6   7   8   9   10  11  12  13
C--   epm = e0 e1+ e2+ e3+ e4+ e5+ e6+ e1- e2-  e3- e4- e5- e6-
C--
C--   qpm(13)   (in) : array with qpm_i values at some (x,q2)
C--   epm(13)  (out) : array with epm_i values at some (x,q2)
C--   nf        (in) : number of active flavours


      implicit double precision (a-h,o-z)

      dimension qpm(*), epm(*)

      epm(1) = qpm(1)

      call sqcQtoE6(qpm(2),epm(2),nf)
      call sqcQtoE6(qpm(8),epm(8),nf)

      return
      end

C     ==================================
      subroutine sqcEPMtoQPM(epm,qpm,nf)
C     ==================================

C--   Transform 13 epm(x,q2) --> qpm(x,q2)
C--
C--   epm = e0 e1+ e2+ e3+ e4+ e5+ e6+ e1- e2-  e3- e4- e5- e6-
C--          1  2   3   4   5   6   7   8   9   10  11  12  13
C--   qpm =  g  d+  u+  s+  c+  b+  t+  d-  u-  s-  c-  b-  t-
C--
C--   epm(13)   (in) : array with epm_i values at some (x,q2)
C--   qpm(13)  (out) : array with qpm_i values at some (x,q2)
C--   nf        (in) : number of active flavours

      implicit double precision (a-h,o-z)

      dimension epm(*), qpm(*)

      qpm(1) = epm(1)

      call sqcEtoQ6(epm(2),qpm(2),nf)
      call sqcEtoQ6(epm(8),qpm(8),nf)

      return
      end

C     =================================
      subroutine sqcQtoE6(qval,eval,nf)
C     =================================

C--   Transform 6 qpm(x,q2) --> epm(x,q2)
C--
C--   qpm --> epm  for i = 1,..,nf
C--   qpm  =  epm  for i = nf+1,..,6
C--
C--   qval(6)  (in): qval(i) = qi+(x,q2) or qi-(x,q2)   i = 1,...,6
C--   eval(6) (out): eval(i) = ei+(x,q2) or ei-(x,q2)   i = 1,...,6
C--   nf       (in): number of active flavours

      implicit double precision (a-h,o-z)

      dimension qval(*), eval(*)

      eval(1) = qval(1)
      do i = 2,nf
        eval(1) = eval(1) + qval(i)
        eval(i) = eval(1) - i*qval(i)
      enddo
      do i = nf+1,6
        eval(i) = qval(i)
      enddo

      return
      end

C     =================================
      subroutine sqcEtoQ6(eval,qval,nf)
C     =================================

C--   Transform 6 epm(x,q2) --> qpm(x,q2)
C--
C--   epm --> qpm  for i = 1,..,nf
C--   epm  =  qpm  for i = nf+1,..,6
C--
C--   eval(6)  (in): eval(i) = ei+(x,q2) or ei-(x,q2)   i = 1,...,6
C--   qval(6) (out): qval(i) = qi+(x,q2) or qi-(x,q2)   i = 1,...,6
C--   nf       (in): number of active flavours

      implicit double precision (a-h,o-z)

      dimension eval(*), qval(*)

      qval(1) = eval(1)
      do i = nf,2,-1
        qval(i) = (qval(1) - eval(i))/i
        qval(1) =  qval(1) - qval(i)
      enddo
      do i = nf+1,6
        qval(i) = eval(i)
      enddo

      return
      end


