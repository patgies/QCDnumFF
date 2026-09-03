C--   This is the file compa.f with floating-point comparisons

C--   logical function lmb_eq(a,b,epsi)
C--   logical function lmb_ne(a,b,epsi)
C--   logical function lmb_ge(a,b,epsi)
C--   logical function lmb_le(a,b,epsi)
C--   logical function lmb_gt(a,b,epsi)
C--   logical function lmb_lt(a,b,epsi)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  MBUTIL comparison routines from compa.f                 */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  MBUTIL comparison routines from compa.f                   */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  MBUTIL comparison routines from compa.f                   */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    int lmb_eq(double a, double b, double epsi);
C-----------------------------------------------------------------------
CXXHFW  #define flmb_eq FC_FUNC(lmb_eq,LMB_EQ)
CXXHFW    int flmb_eq(double*,double*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int lmb_eq(double a, double b, double epsi)
CXXWRP  {
CXXWRP    return flmb_eq(&a, &b, &epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      logical function lmb_eq(a,b,epsi)
C     =================================

C--   Author: Michiel Botje h24@nikhef.nl   12-12-18

      implicit double precision(a-h,o-z)

C--   Absolute or relative, thats the question
      if(epsi.ge.0) then
        eps = epsi
      else
        eps = abs(epsi)*max(abs(a),abs(b))
      endif
C--   Compare
      lmb_eq = abs(a-b).le.eps

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int lmb_ne(double a, double b, double epsi);
C-----------------------------------------------------------------------
CXXHFW  #define flmb_ne FC_FUNC(lmb_ne,LMB_NE)
CXXHFW    int flmb_ne(double*,double*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int lmb_ne(double a, double b, double epsi)
CXXWRP  {
CXXWRP    return flmb_ne(&a, &b, &epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      logical function lmb_ne(a,b,epsi)
C     =================================

C--   Author: Michiel Botje h24@nikhef.nl   12-12-18

      implicit double precision(a-h,o-z)

C--   Absolute or relative, thats the question
      if(epsi.ge.0) then
        eps = epsi
      else
        eps = abs(epsi)*max(abs(a),abs(b))
      endif
C--   Compare
      lmb_ne = abs(a-b).gt.eps

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int lmb_ge(double a, double b, double epsi);
C-----------------------------------------------------------------------
CXXHFW  #define flmb_ge FC_FUNC(lmb_ge,LMB_GE)
CXXHFW    int flmb_ge(double*,double*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int lmb_ge(double a, double b, double epsi)
CXXWRP  {
CXXWRP    return flmb_ge(&a, &b, &epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      logical function lmb_ge(a,b,epsi)
C     =================================

C--   Author: Michiel Botje h24@nikhef.nl   12-12-18

      implicit double precision(a-h,o-z)

C--   Absolute or relative, thats the question
      if(epsi.ge.0) then
        eps = epsi
      else
        eps = abs(epsi)*max(abs(a),abs(b))
      endif
C--   Compare
      lmb_ge = (abs(a-b).le.eps) .or. ((a-b).gt.0.D0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int lmb_le(double a, double b, double epsi);
C-----------------------------------------------------------------------
CXXHFW  #define flmb_le FC_FUNC(lmb_le,LMB_LE)
CXXHFW    int flmb_le(double*,double*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int lmb_le(double a, double b, double epsi)
CXXWRP  {
CXXWRP    return flmb_le(&a, &b, &epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      logical function lmb_le(a,b,epsi)
C     =================================

C--   Author: Michiel Botje h24@nikhef.nl   12-12-18

      implicit double precision(a-h,o-z)

C--   Absolute or relative, thats the question
      if(epsi.ge.0) then
        eps = epsi
      else
        eps = abs(epsi)*max(abs(a),abs(b))
      endif
C--   Compare
      lmb_le = (abs(a-b).le.eps) .or. ((a-b).lt.0.D0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int lmb_gt(double a, double b, double epsi);
C-----------------------------------------------------------------------
CXXHFW  #define flmb_gt FC_FUNC(lmb_gt,LMB_GT)
CXXHFW    int flmb_gt(double*,double*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int lmb_gt(double a, double b, double epsi)
CXXWRP  {
CXXWRP    return flmb_gt(&a, &b, &epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      logical function lmb_gt(a,b,epsi)
C     =================================

C--   Author: Michiel Botje h24@nikhef.nl   12-12-18

      implicit double precision(a-h,o-z)

C--   Absolute or relative, thats the question
      if(epsi.ge.0) then
        eps = epsi
      else
        eps = abs(epsi)*max(abs(a),abs(b))
      endif
C--   Compare
      lmb_gt = (abs(a-b).gt.eps) .and. ((a-b).gt.0.D0)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    int lmb_lt(double a, double b, double epsi);
C-----------------------------------------------------------------------
CXXHFW  #define flmb_lt FC_FUNC(lmb_lt,LMB_LT)
CXXHFW    int flmb_lt(double*,double*,double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int lmb_lt(double a, double b, double epsi)
CXXWRP  {
CXXWRP    return flmb_lt(&a, &b, &epsi);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      logical function lmb_lt(a,b,epsi)
C     =================================

C--   Author: Michiel Botje h24@nikhef.nl   12-12-18

      implicit double precision(a-h,o-z)

C--   Absolute or relative, thats the question
      if(epsi.ge.0) then
        eps = epsi
      else
        eps = abs(epsi)*max(abs(a),abs(b))
      endif
C--   Compare
      lmb_lt = (abs(a-b).gt.eps) .and. ((a-b).lt.0.D0)

      return
      end
