
C--   This is the file vectors.f with MBUTIL vector operations
C--
c--   subroutine smb_Ifill(ia,n,ival)
C--   subroutine smb_Vfill(a,n,val)
C--   subroutine smb_Vmult(a,n,val)
C--   subroutine smb_Vcopy(a,b,n)
C--   subroutine smb_Icopy(ia,ib,n)
C--   subroutine smb_Vitod(ia,b,n)
C--   subroutine smb_Vdtoi(a,ib,n)
C--   subroutine smb_VaddV(a,b,c,n)
C--   subroutine smb_VminV(a,b,c,n)
C--   double precision function dmb_VdotV(a,b,n)
C--   double precision function dmb_Vnorm(m,a,n)
C--   logical function lmb-Vcomp(a,b,n,epsi)
C--   double precision function dmb_Vpsum(a,w,n)
C--   subroutine smbAddPairs(w,n)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  MBUTIL vector routines from vectors.f                   */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  MBUTIL vector routines from vectors.f                     */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  MBUTIL vector routines from vectors.f                     */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    void smb_ifill(int *ia, int n, int ival);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_ifill FC_FUNC(smb_ifill,SMB_IFILL)
CXXHFW    void fsmb_ifill(int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_ifill(int *ia, int n, int ival)
CXXWRP  {
CXXWRP    fsmb_ifill(ia, &n, &ival);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===============================
      subroutine smb_Ifill(ia,n,ival)
C     ===============================

C--   Set n elements of IA to ival
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension ia(*)

      if(n.le.0) stop 'SMB_IFILL(ia,n,ival) input n is zero or negative'

      do i = 1,n
        ia(i) = ival
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_vfill(double *a, int n, double val);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_vfill FC_FUNC(smb_vfill,SMB_VFILL)
CXXHFW    void fsmb_vfill(double*,int*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_vfill(double *a, int n, double val)
CXXWRP  {
CXXWRP    fsmb_vfill(a, &n, &val);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      subroutine smb_Vfill(a,n,val)
C     =============================

C--   Set a_i to val
C--
C--   a(n)        (in) : vector
C--   val         (in) : value to be set
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension a(*)

      if(n.le.0) stop 'SMB_VFILL(a,n,val) input n is zero or negative'

      do i = 1,n
        a(i) = val
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_vmult(double *a, int n, double val);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_vmult FC_FUNC(smb_vmult,SMB_VMULT)
CXXHFW    void fsmb_vmult(double*,int*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_vmult(double *a, int n, double val)
CXXWRP  {
CXXWRP    fsmb_vmult(a, &n, &val);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      subroutine smb_Vmult(a,n,val)
C     =============================

C--   Multiply a_i by val
C--
C--   a(n)        (in) : vector
C--   val         (in) : multiplication factor
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension a(*)

      if(n.le.0) stop 'SMB_VMULT(a,n,val) input n is zero or negative'

      do i = 1,n
        a(i) = val*a(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_vcopy(double *a, double* b, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_vcopy FC_FUNC(smb_vcopy,SMB_VCOPY)
CXXHFW    void fsmb_vcopy(double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_vcopy(double *a, double* b, int n)
CXXWRP  {
CXXWRP    fsmb_vcopy(a, b, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     ===========================
      subroutine smb_Vcopy(a,b,n)
C     ===========================

C--   Set b = a
C--
C--   a(n)        (in) : vector
C--   b(n)       (out) : copy of a
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension a(*), b(*)

      if(n.le.0) stop 'SMB_VCOPY(a,b,n) input n is zero or negative'

      do i = 1,n
        b(i) = a(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_icopy(int *ia, int* ib, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_icopy FC_FUNC(smb_icopy,SMB_ICOPY)
CXXHFW    void fsmb_icopy(int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_icopy(int *ia, int* ib, int n)
CXXWRP  {
CXXWRP    fsmb_icopy(ia, ib, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      subroutine smb_Icopy(ia,ib,n)
C     =============================

C--   Set ib = ia
C--
C--   ia(n)        (in) : vector
C--   ib(n)       (out) : copy of ia
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-02-21

      implicit double precision (a-h,o-z)

      dimension ia(*), ib(*)

      if(n.le.0) stop 'SMB_ICOPY(ia,ib,n) input n is zero or negative'

      do i = 1,n
        ib(i) = ia(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_vitod(int *ia, double *b, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_vitod FC_FUNC(smb_vitod,SMB_VITOD)
CXXHFW    void fsmb_vitod(int*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_vitod(int *ia, double *b, int n)
CXXWRP  {
CXXWRP    fsmb_vitod(ia, b, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     ============================
      subroutine smb_Vitod(ia,b,n)
C     ============================

C--   Set b = dble(ia)
C--
C--   ia(n)       (in) : vector
C--   b(n)       (out) : copy of ia
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension ia(*), b(*)

      if(n.le.0) stop 'SMB_VITOD(ia,b,n) input n is zero or negative'

      do i = 1,n
        b(i) = dble(ia(i))
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_vdtoi(double *a, int *ib, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_vdtoi FC_FUNC(smb_vdtoi,SMB_VDTOI)
CXXHFW    void fsmb_vdtoi(double*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_vdtoi(double *a, int *ib, int n)
CXXWRP  {
CXXWRP    fsmb_vdtoi(a, ib, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     ============================
      subroutine smb_Vdtoi(a,ib,n)
C     ============================

C--   Set ib = int(a)
C--
C--   a(n)        (in) : vector
C--   ib(n)      (out) : copy of a
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension a(*), ib(*)

      if(n.le.0) stop 'SMB_VDTOI(a,ib,n) input n is zero or negative'

      do i = 1,n
        ib(i) = int(a(i))
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_vaddv(double *a, double *b, double *c, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_vaddv FC_FUNC(smb_vaddv,SMB_VADDV)
CXXHFW    void fsmb_vaddv(double*,double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_vaddv(double *a, double *b, double *c, int n)
CXXWRP  {
CXXWRP    fsmb_vaddv(a, b, c, &n);
CXXWRP  }
C-----------------------------------------------------------------------


C     =============================
      subroutine smb_VaddV(a,b,c,n)
C     =============================

C--   Vector sum  c = a + b
C--
C--   a(n), b(n)  (in) : vectors
C--   c(n)       (out) : vector (allowed to overwrite a or b)
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension a(*), b(*), c(*)

      if(n.le.0) stop 'SMB_VADDV(a,b,c,n) input n is zero or negative'

      do i = 1,n
        c(i) = a(i) + b(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void smb_vminv(double *a, double *b, double *c, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fsmb_vminv FC_FUNC(smb_vminv,SMB_VMINV)
CXXHFW    void fsmb_vminv(double*,double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void smb_vminv(double *a, double *b, double *c, int n)
CXXWRP  {
CXXWRP    fsmb_vminv(a, b, c, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     =============================
      subroutine smb_VminV(a,b,c,n)
C     =============================

C--   Vector difference  c = a - b
C--
C--   a(n), b(n)  (in) : vectors
C--   c(n)       (out) : vector (allowed to overwrite a or b)
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension a(*), b(*), c(*)

      if(n.le.0) stop 'SMB_VMINV(a,b,c,n) input n is zero or negative'

      do i = 1,n
        c(i) = a(i) - b(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dmb_vdotv(double *a, double *b, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fdmb_vdotv FC_FUNC(dmb_vdotv,DMB_VDOTV)
CXXHFW    double fdmb_vdotv(double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double dmb_vdotv(double *a, double *b, int n)
CXXWRP  {
CXXWRP    return fdmb_vdotv(a, b, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================================
      double precision function dmb_VdotV(a,b,n)
C     ==========================================

C--   Vector dot product a.b
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension a(*), b(*)

      if(n.le.0) stop 'DMB_VDOTV(a,b,n) input n is zero or negative'

      sum = 0.D0
      do i = 1,n
        sum = sum + a(i)*b(i)
      enddo
      dmb_VdotV = sum

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dmb_vnorm(int m, double *a, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fdmb_vnorm FC_FUNC(dmb_vnorm,DMB_VNORM)
CXXHFW    double fdmb_vnorm(int*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double dmb_vnorm(int m, double *a, int n)
CXXWRP  {
CXXWRP    return fdmb_vnorm(&m, a, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================================
      double precision function dmb_Vnorm(m,a,n)
C     ==========================================

C--   Various norms of a
C--
C--   Author: Michiel Botje h24@nikhef.nl   02-11-18

      implicit double precision (a-h,o-z)

      dimension a(*)

      if(n.le.0) stop 'DMB_VNORM(m,a,n) input n is zero or negative'

      if(m.eq.2) then
C--     Euclidian norm
        sum = 0.D0
        do i = 1,n
          sum = sum + a(i)*a(i)
        enddo
        dmb_Vnorm = sqrt(sum)
      elseif(m.eq.0) then
C--     Max-norm
        vmax = 0.D0
        do i = 1,n
          vmax = max(vmax,abs(a(i)))
        enddo
        dmb_Vnorm = vmax
      elseif(m.eq.1) then
C--     City-block norm
        sum = 0.D0
        do i = 1,n
          sum = sum + abs(a(i))
        enddo
        dmb_Vnorm = sum
      elseif(m.gt.2) then
C--     m-norm
        sum = 0.D0
        do i = 1,n
          sum = sum + abs(a(i))**m
        enddo
        dmb_Vnorm = sum**(1.D0/m)
      else
        stop 'DMB_VNORM(m,a,n) input m is negative'
      endif

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int lmb_vcomp(double *a, double *b, int n, double epsi);
C-----------------------------------------------------------------------
CXXHFW  #define flmb_vcomp FC_FUNC(lmb_vcomp,LMB_VCOMP)
CXXHFW    int flmb_vcomp(double*,double*,int*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int lmb_vcomp(double *a, double *b, int n, double epsi)
CXXWRP  {
CXXWRP    return flmb_vcomp(a, b, &n, &epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     ======================================
      logical function lmb_Vcomp(a,b,n,epsi)
C     ======================================

C--   True if a == b within tolerance epsi
C--
C--   Author: Michiel Botje h24@nikhef.nl   27-02-19

      implicit double precision (a-h,o-z)

      logical lmb_ne

      dimension a(*), b(*)

      if(n.le.0) stop 'LMB_VCOMP(a,b,n,eps) input n is zero or negative'

      lmb_Vcomp = .false.
      do i = 1,n
        if( lmb_ne( a(i), b(i), epsi ) ) return
      enddo
      lmb_Vcomp = .true.

      return
      end

C-----------------------------------------------------------------------
CXXHDR    double dmb_vpsum(double *a, double *w, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fdmb_vpsum FC_FUNC(dmb_vpsum,DMB_VPSUM)
CXXHFW    double fdmb_vpsum(double*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double dmb_vpsum(double *a, double *w, int n)
CXXWRP  {
CXXWRP    return fdmb_vpsum(a, w, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================================
      double precision function dmb_Vpsum(a,w,n)
C     ==========================================

C--   Pairwise addition of the first n elements of a. The user must
C--   provide a work array w, dimensioned to at least n+1 in
C--   the calling routine.

C--   Author: Michiel Botje h24@nikhef.nl   24-12-21

      implicit double precision (a-h,o-z)

      dimension a(*), w(*)

      if(n.le.0) stop 'DMB_VPSUM(a,w,n) input n is zero or negative'

      call smb_Vcopy(a,w,n)
      m = n
      do while(m.gt.1)
        call smbAddPairs(w,m)
      enddo
      dmb_Vpsum = w(1)

      return
      end

C     ===========================
      subroutine smbAddPairs(w,n)
C     ===========================

C--   Compute pairwise sums w(i)+w(i+1) and store these sums in the
C     first j = n/2 elements of w. On exit, n is set to j.
C--   In the calling routine w should be dimensioned to at least n+1
C--   to correctly handle an odd number of terms.
C--   The original content of w is destroyed in the process.

C--   Author: Michiel Botje h24@nikhef.nl   24-12-21

      implicit double precision (a-h,o-z)

      dimension w(*)

      w(n+1) = 0.D0
      j      = 0
      do i = 1,n,2
        j    = j+1
        w(j) = w(i)+w(i+1)
      enddo
      n = j

      return
      end


