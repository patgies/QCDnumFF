
C--   This is the file usrGlobalId with conversions to global addressing

C--   integer function iPdftab(iset,id)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM address function from usrGlobalId.f              */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM address function from usrGlobalId.f                */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM address function from usrGlobalId.f                */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C=======================================================================
C==== User callable internal pdf addres function
C=======================================================================

C-----------------------------------------------------------------------
CXXHDR    int ipdftab(int iset, int id);
C-----------------------------------------------------------------------
CXXHFW  #define fipdftab FC_FUNC(ipdftab,IPDFTAB)
CXXHFW    int fipdftab(int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  int ipdftab(int iset, int id)
CXXWRP  {
CXXWRP    return fipdftab(&iset, &id);
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================
      integer function iPdftab(iset,id)
C     =================================

C--   Returns minus the global pdf identifier in the internal store
C--
C--   iset  (in): pdf set identifier [1,mset0]
C--   id    (in): pdf identifier     [0,mpdf0]
C--
C--   Returns -1990XX if iset not in range   (ierr = 1)
C--           -2XX099 if id   not in range   (ierr = 2)
C--           -3XX0YY if iset does not exist (ierr = 3)
C--           -4XX0YY if id   does not exist (ierr = 4)
C--           -5XX0YY if id   is not filled  (ierr = 5)

      implicit double precision (a-h,o-z)      
      
      include 'qcdnum.inc'
      include 'qstor7.inc'

      logical lqcIsFilled

      idg(ierr,iset,id) = 100000*ierr+1000*iset+id      !inline function

      ierr = 0

C--   Bring input iset in range (2-digit positive number)
      if(iset.ge.1. .and. iset.le.mset0) then
        jset = iset
      else
        jset = 99
        ierr = 1
      endif
C--   Bring input id in range (2-digit positive number)
      if(id.ge.0 .and. id.le.mpdf0) then
        jd   = id
      else
        jd   = 99
        ierr = 2
      endif

C--   Error: iset and/or id were out of range
      if(ierr.ne.0) then
        iPdftab = -idg(ierr,jset,jd)
        return
      endif

C--   Find out if jset exists
      if(jset.ge.1 .and. jset.le.mset0) then
        kset = isetf7(jset)                  !kset might be 0
      else
        kset = 0
      endif

      if(kset.eq.0) then
C--     Error: jset does not exist
        ierr    =  3
        iPdftab = -idg(ierr,jset,jd)
        return
      elseif(jd.ge.ifrst7(jset) .and. jd.le.ilast7(jset)) then
C--     OK: jd does exist
        kd = jd - ifrst7(jset) + 501
      else
C--     Error: jd does not exist
        ierr    =  4
        iPdftab = -idg(ierr,jset,jd)
        return
      endif

C--   Get global identifier
      kdg = idg(ierr,kset,kd)

      if(.not.lqcIsFilled(stor7,kdg)) then
C--     Error: pdf table is empty
        ierr = 5
        iPdftab = -idg(ierr,jset,jd)
      else
C--     OK: pdf table is filled
        iPdftab = -kdg
      endif

      return
      end

