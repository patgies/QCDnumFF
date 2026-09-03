
C--   This is the file obsolete.f with obsolete routines

C--*   integer function Nflavs(iq,ithrs)                     [17-01-15]
C--   subroutine SetCut(xmi,qmi,qma,dum)                    [17-01-15]
C--   subroutine GetCut(xmi,qmi,qma,dum)                    [17-01-15]
C--   subroutine PdfExt(subr,jset,n,offset,epsi)            [17-01-14]
C--   double precision function fsumxq(jset,def,xx,qq,jchk) [17-01-13]
C--   subroutine PdfTab(jset,def,xx,nx,qq,nq,pdf,jchk)      [17-01-13]
C--   subroutine PdfLst(jset,def,x,q,f,n,jchk)              [17-01-13]
C--   subroutine PushPar                                    [17-01-13]
C--   subroutine PullPar                                    [17-01-13]
C--   subroutine fpdfxq(jset,xx,qq,pdf,jchk)                [17-01-13]
C--   subroutine fpdfij(jset,ix,iq,pdf,jchk)                [17-01-13]
C--   double precision function fsnsxq(jset,id,xx,qq,jchk)  [17-01-13]
C--   double precision function fsnsij(jset,id,ix,iq,jchk)  [17-01-13]
C--   subroutine SavePar                                    [obsolete]
C--   subroutine UnsaveP                                    [obsolete]
C--   integer function Nflavor(iq)                          [obsolete]
C--   subroutine EvFCopy(w,id,def,jset)                     [obsolete]
C--   logical function chkpdf(itype)                        [obsolete]
C--   double precision function GetAlfN(jq,n,ierr)          [obsolete]
C--   subroutine PdfInp(subr,jset,offset,epsi,nwlast)       [obsolete]
C--   logical function lpassc(x,qmu2,ifail,jchk)            [obsolete]
C--   subroutine TabDump(w,lun,file,key)                    [obsolete]
C--   subroutine TabRead(w,nw,lun,file,key,nwords,ierr)     [obsolete]
C--   subroutine SetWpar(w,par,n)                           [obsolete]
C--   subroutine GetWpar(w,par,n)                           [obsolete]
C--   subroutine BookTab(w,nw,jtypes,nwds)                  [obsolete]
C--   function alfunc(iord,as0,r20,r2,ierr)                 [obsolete]
C--   double precision function AsEvol
C--                  (iord,as0,r20,r2,iqcdnum,nfout,ierr)   [obsolete]
C--   subroutine Evolve(func,def,iq0,epsi)                  [obsolete]
C--   subroutine EvolFF(func,def,iq0,epsi)                  [obsolete}
C--   subroutine NSevol(ityp,func,idf,iq0,epsi)             [obsolete]
C--   subroutine EvolNS(ityp,func,idf,iq0,epsi)             [obsolete]
C--   subroutine SGEvol(funf,idf,fung,idg,iq0,epsi)         [obsolete]
C--   subroutine EvolSG(funf,idf,fung,idg,iq0,epsi)         [obsolete]
C--   subroutine getids(idmin,idmax,nwords)                 [obsolete]
C--   subroutine dumpwt(lun,file)                           [obsolete]
C--   subroutine FastFac(id,funxq)                          [obsolete]
C--   subroutine FtimesW(w,fun,jd1,id2,iadd)                [obsolete]
C--   double precision function FcrossC(w,id,idf,ix,iq)     [obsolete]
C--   subroutine SetPdf(itype)                              [obsolete]
C--   subroutine GetPdf(itype)                              [obsolete]
C--   subroutine FastPdf(id,coef)                           [obsolete]
C--   subroutine FastFcC(w,id,id1,id2)                      [obsolete]
C--   subroutine FastAdd(id1,id2)                           [obsolete]
C--   subroutine setabq(aq,bq)                              [obsolete]
C--   subroutine getabq(aq,bq)                              [obsolete]
C--   subroutine setthr(nfix,q2c,q2b,q2t)                   [obsolete]
C--   subroutine getthr(nfix,q2c,q2b,q2t,iqcdnum)           [obsolete]
C--   subroutine pdfsum(id,wt,n,idout)                      [obsolete]
C--   double precision function pdfval(idf,xx,qq,mode)      [obsolete]
C--   double precision function pgluon(x,qmu2,jchk)         [obsolete]
C--   double precision function onepdf(x,qmu2,def,jchk)     [obsolete]
C--   subroutine allpdf(x,qmu2,pdf,imode)                   [obsolete]
C--   subroutine pdfsxq(xx,qq,pdf,jchk)                     [obsolete]
C--   subroutine pdfsij(ix,iq,pdf,jchk)                     [obsolete]
C--   double precision function sgnsxq(id,xx,qq,jchk)       [obsolete]
C--   double precision function sgnsij(id,ix,iq,jchk)       [obsolete]
C--   double precision function pfunxq(id,xx,qq,jchk)       [obsolete]
C--   double precision function pfunij(id,ix,iq,jchk)       [obsolete]
C--   double precision function psumxq(def,xx,qq,jchk)      [obsolete]
C--   double precision function psumij(def,ix,iq,jchk)      [obsolete]
C--   subroutine stfval(istf,ityp,id,x,q,f,n,mode)          [obsolete]
C--   subroutine StrFun(istf,user,x,q,f,n,mode)             [obsolete]
C--   subroutine StampIt('string')                          [obsolete]
C--   integer function idSplij(string)                      [obsolete]

*C     =================================
*      integer function Nflavs(iq,ithrs)
*C     =================================
*
*C--   [obsolete]
*
*      implicit double precision (a-h,o-z)
*
*      character*80 subnam
*      data subnam  /'NFLAVS ( IQ, ITHRESH )'/
*
*      idum   = iq
*      idum   = ithrs
*      nflavs = 0
*
*      call sqcErrMsg(subnam,
*     +              'NFLAVS obsolete, please use NFRMIQ instead')
*
*      return
*      end

C     ==================================
      subroutine SetCut(xmi,qmi,qma,dum)
C     ==================================

C--   [obsolete]

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'SETCUT ( XMI, Q2MI, Q2MA, DUMMY )'/

      ddum = xmi
      ddum = qmi
      ddum = qma
      ddum = dum

      call sqcErrMsg(subnam,
     +              'SETCUT obsolete, please use SETLIM instead')

      return
      end

C     ==================================
      subroutine GetCut(xmi,qmi,qma,dum)
C     ==================================

C--   [obsolete]

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'GETCUT ( XMI, Q2MI, Q2MA, DUMMY )'/

      ddum = xmi
      ddum = qmi
      ddum = qma
      ddum = dum

      call sqcErrMsg(subnam,
     +              'GETCUT obsolete, please use GETLIM instead')

      return
      end

C     ==========================================
      subroutine PdfExt(subr,jset,n,offset,epsi)
C     ==========================================

C--   [obsolete]

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'PDFEXT ( SUBR, ISET, N, DELTA, EPSI )'/

      ddum = subr
      idum = jset
      idum = n
      ddum = offset
      ddum = epsi

      call sqcErrMsg(subnam,
     +              'PDFEXT obsolete, please use EXTPDF instead')

      return
      end

C     =====================================================
      double precision function fsumxq(jset,def,xx,qq,jchk)
C     =====================================================

C--   [obsolete]

      implicit double precision (a-h,o-z)

      dimension def(-6:6)

      character*80 subnam
      data subnam /'FSUMXQ ( ISET, C, X, QMU2, ICHK )'/

      idum = jset
      ddum = def(1)
      ddum = xx
      ddum = qq
      idum = jchk

      call sqcErrMsg(subnam,
     +              'FSUMXQ obsolete, please use SUMFXQ instead')

      fsumxq = 0.D0

      return
      end

C     =====================================================
      double precision function fsumij(jset,def,xx,qq,jchk)
C     =====================================================

C--   [obsolete]

      implicit double precision (a-h,o-z)

      dimension def(-6:6)

      character*80 subnam
      data subnam /'FSUMIJ ( ISET, C, IX, IQ, ICHK )'/

      idum = jset
      ddum = def(1)
      ddum = xx
      ddum = qq
      idum = jchk

      call sqcErrMsg(subnam,
     +              'FSUMIJ obsolete, please use SUMFIJ instead')

      fsumij = 0.D0

      return
      end

C     ================================================
      subroutine PdfTab(jset,def,xx,nx,qq,nq,pdf,jchk)
C     ================================================

C--   [obsolete]

      implicit double precision (a-h,o-z)

      dimension def(-6:6)
 
      dimension xx(nx), qq(nq), pdf(nx,nq)
      
      character*80 subnam
      data subnam /'PDFTAB ( ISET, DEF, XX, NX, QQ, NQ, PDF, ICHK )'/

      idum = jset
      ddum = def(0)
      ddum = xx(nx)
      ddum = qq(nq)
      ddum = pdf(1,1)
      idum = jchk

      call sqcErrMsg(subnam,
     +              'PDFTAB obsolete, please use FTABLE instead')

      return
      end

C     ========================================      
      subroutine PdfLst(jset,def,x,q,f,n,jchk)
C     ========================================

C--   [obsolete]

      implicit double precision (a-h,o-z)
      
      dimension def(-6:6)
      dimension x(*), q(*), f(*)
      
      character*80 subnam
      data subnam /'PDFLST ( ISET, DEF, X, Q, F, N, ICHK )'/

      idum = jset
      ddum = def(0)
      ddum = x(1)
      ddum = q(1)
      ddum = f(1)
      idum = n
      idum = jchk

      call sqcErrMsg(subnam,
     +              'PDFLST obsolete, please use FFLIST instead')

      return
      end

C     ==================
      subroutine PushPar
C     ==================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam  /'PUSHPAR'/

      call sqcErrMsg(subnam,
     +              'PUSHPAR obsolete, please use PUSHCP instead')

      return
      end

C     ==================
      subroutine PullPar
C     ==================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam  /'PULLPAR'/

      call sqcErrMsg(subnam,
     +              'PULLPAR obsolete, please use PULLCP instead')

      return
      end


C     ======================================
      subroutine fpdfxq(jset,xx,qq,pdf,jchk)
C     ======================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'FPDFXQ ( ISET, X, QMU2, PDFS, ICHK )'/

      dimension pdf(-6:6)

      idum = jset
      ddum = xx
      ddum = qq
      ddum = pdf(0)
      idum = jchk

      call sqcErrMsg(subnam,
     +              'FPDFXQ obsolete, please use ALLFXQ instead')

      return
      end

C     ======================================
      subroutine fpdfij(jset,ix,iq,pdf,jchk)
C     ======================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'FPDFIJ ( ISET, IX, IQ, PDFS, ICHK )'/

      dimension pdf(-6:6)

      idum = jset
      ddum = ix
      ddum = iq
      ddum = pdf(0)
      idum = jchk

      call sqcErrMsg(subnam,
     +              'FPDFIJ obsolete, please use ALLFIJ instead')

      return
      end

C     ====================================================
      double precision function fsnsxq(jset,id,xx,qq,jchk)
C     ====================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'FSNSXQ ( ISET, ID, X, QMU2, ICHK )'/

      idum = jset
      idum = id
      ddum = xx
      ddum = qq
      idum = jchk

      fsnsxq = 0.D0

      call sqcErrMsg(subnam,
     +              'FSNSXQ obsolete, please use BVALXQ instead')

      return
      end

C     ====================================================
      double precision function fsnsij(jset,id,ix,iq,jchk)
C     ====================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'FSNSIJ ( ISET, ID, IX, IQ, ICHK )'/

      idum = jset
      idum = id
      idum = ix
      idum = iq
      idum = jchk

      fsnsij = 0.D0

      call sqcErrMsg(subnam,
     +              'FSNSIJ obsolete, please use BVALIJ instead')

      return
      end

C     ==================
      subroutine SavePar
C     ==================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam  /'SAVEPAR'/

      call sqcErrMsg(subnam,
     +              'SAVEPAR obsolete, please use PUSHPAR instead')

      return
      end

C     ==================
      subroutine UnsaveP
C     ==================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam  /'UNSAVEP'/

      call sqcErrMsg(subnam,
     +              'UNSAVEP obsolete, please use PULLPAR instead')

      return
      end

C     ============================
      integer function Nflavor(iq)
C     ============================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam  /'NFLAVOR ( IQ )'/

      idum    = iq
      nflavor = 0

      call sqcErrMsg(subnam,
     +              'NFLAVOR obsolete, please use NFLAVS instead')

      return
      end

C     =================================
      subroutine EvFCopy(w,id,def,jset)
C     =================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'EVFCOPY ( W, IDF, DEF, ISET )'/

      dimension w(*)
      dimension id(*), def(-6:6,*)

      ddum = w(1)
      idum = id(1)
      ddum = def(0,1)
      idum = jset

      call sqcErrMsg(subnam,
     +              'EVFCOPY obsolete, please use EVPCOPY instead')

      return
      end

C     ==============================
      logical function chkpdf(itype)
C     ==============================      
      
C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'CHKPDF ( ITYPE )'/

      jtype  =  itype
      chkpdf = .false.

      call sqcErrMsg(subnam,
     +              'CHKPDF obsolete, please use NPTABS instead')

      return
      end

C     ============================================
      double precision function GetAlfN(jq,n,ierr)
C     ============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'GETALFN ( IQ, N, IERR )'/

      idum    = jq
      idum    = n
      idum    = ierr
      getalfn = 0.D0

      call sqcErrMsg(subnam,
     +              'GETALFN obsolete, please use ALTABN instead')

      return
      end

C     ===============================================
      subroutine PdfInp(subr,jset,offset,epsi,nwlast)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'PDFINP ( SUBR, ISET, DELTA, EPSI, NWDS )'/

      ddum = subr
      idum = jset
      ddum = offset
      ddum = epsi
      idum = nwlast

      call sqcErrMsg(subnam,
     +              'PDFINP obsolete, please use PDFEXT instead')

      return
      end

C     ==========================================
      logical function lpassc(x,qmu2,ifail,jchk)
C     ==========================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      character*80 subnam
      data subnam /'LPASSC ( X, QMU2, IFAIL, ICHK )'/

      lpassc = .false.
      dum    = x
      dum    = qmu2
      ifail  = 0
      idum   = jchk

      call sqcErrMsg(subnam,'LPASSC obsolete, has been removed')

      return
      end

C     =================================================
      subroutine TabRead(w,nw,lun,file,key,nwords,ierr)
C     =================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      dimension w(*)
      character*(*) file, key 

      character*80 subnam
      data subnam /'TABREAD ( W, NW, LUN, FILE, KEY, NWORDS, IERR )'/

      dum = w(1)
      nnn = nw
      lll = lun
      iii = len(file)
      jjj = len(key)
      nww = nwords
      iee = ierr

      call sqcErrMsg(subnam,
     +              'TABREAD obsolete, please use READTAB instead')


      return
      end

C     ==================================
      subroutine TabDump(w,lun,file,key)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      dimension w(*)
      character*(*) file, key 

      character*80 subnam
      data subnam /'TABDUMP ( W, LUN, FILE, KEY )'/

      dum = w(1)
      lll = lun
      iii = len(file)
      jjj = len(key)

      call sqcErrMsg(subnam,
     +              'TABDUMP obsolete, please use DUMPTAB instead')


      return
      end

C     ===========================
      subroutine SetWpar(w,par,n)
C     ===========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SETWPAR ( W, PAR, N )'/

      dimension w(*), par(*)

      dum1 = w(1)
      dum2 = par(1)
      nw   = n

      call sqcErrMsg(subnam,
     +              'SETWPAR obsolete, please use SETPARW instead')

      return
      end

C     ===========================
      subroutine GetWpar(w,par,n)
C     ===========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'GETWPAR ( W, PAR, N )'/

      dimension w(*), par(*)

      dum1 = w(1)
      dum2 = par(1)
      nw   = n

      call sqcErrMsg(subnam,
     +              'GETWPAR obsolete, please use GETPARW instead')

      return
      end

C     ====================================
      subroutine BookTab(w,nw,jtypes,nwds)
C     ====================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'BOOKTAB ( W, NW, JTYPES, NWDS )'/

      dimension w(*), jtypes(*)

      dum  = w(1)
      mw   = nw
      jt   = jtypes(1)
      mm   = nwds

      call sqcErrMsg(subnam,
     +              'BOOKTAB obsolete, please use MAKETAB instead')

      return
      end


C     ======================================================
      double precision function alfunc(iord,as0,r20,r2,ierr)
C     ======================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'ALFUNC ( IORD, AS0, R20, R2, IERR )'/

      alfunc = 0.D0
      idum   = iord
      dum    = as0
      dum    = r20
      dum    = r2
      ierr   = 999

      call sqcErrMsg(subnam,
     +              'Alfunc obsolete, pls use EvolAs instead')

      return
      end
      
C     ==========================================================
      double precision function AsEvol
     +                      (iord,as0,r20,r2,iqcdnum,nfout,ierr)
C     ==========================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam 
     +  /'ASEVOL ( IORD, AS0, R20, R2, INTERN, NFOUT, IERR )'/

      asevol = 0.D0
      idum   = iord
      dum    = as0
      dum    = r20
      dum    = r2
      idum   = iqcdnum
      nfout  = 0
      ierr   = 999
           
      call sqcErrMsg(subnam,
     +              'AsEvol obsolete, pls use EvolAs instead')

      return
      end      
      
C     ====================================
      subroutine Evolve(func,def,iq0,epsi)
C     ====================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'EVOLVE ( FUNC, DEF, IQ0, EPSI )'/

      dum  = func
      dum  = def
      idum = iq0
      dum  = epsi

      call sqcErrMsg(subnam,
     +              'Evolve obsolete, pls use EvolFG instead')

      return
      end

C     ====================================
      subroutine EvolFF(func,def,iq0,epsi)
C     ====================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'EVOLFF ( FUNC, DEF, IQ0, EPSI )'/

      dum  = func
      dum  = def
      idum = iq0
      dum  = epsi

      call sqcErrMsg(subnam,
     +              'EVOLFF obsolete, pls use EVOLFG instead')

      return
      end   
      
C     =========================================
      subroutine NSevol(ityp,func,idf,iq0,epsi)
C     =========================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'NSEVOL ( ITYP, FUNC, IDF, IQ0, EPSI )'/

      call sqcErrMsg(subnam,'NSEVOL obsolete, use EVOLFG instead')

      idum = ityp
      dum  = func
      idum = idf
      idum = iq0
      epsi = 1.D11

      return
      end

C     =========================================
      subroutine EvolNS(ityp,func,idf,iq0,epsi)
C     =========================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'EVOLNS ( ITYP, FUNC, IDF, IQ0, EPSI )'/

      call sqcErrMsg(subnam,'EVOLNS obsolete, use EVOLFG instead')

      idum = ityp
      dum  = func
      idum = idf
      idum = iq0
      epsi = 1.D11

      return
      end
            
C     =============================================
      subroutine SGEvol(funf,idf,fung,idg,iq0,epsi)
C     =============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SGEVOL ( SFUN, IDS, GFUN, IDG, IQ0, EPSI )'/

      call sqcErrMsg(subnam,'SGEVOL obsolete, use EVOLFG instead')

      dum  = funf
      idum = idf
      dum  = fung
      idum = idg
      idum = iq0
      epsi = 1.D11

      return
      end

C     =============================================
      subroutine EvolSG(funf,idf,fung,idg,iq0,epsi)
C     =============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'EVOLSG ( SFUN, IDS, GFUN, IDG, IQ0, EPSI )'/

      call sqcErrMsg(subnam,'EVOLSG obsolete, use EVOLFG instead')

      dum  = funf
      idum = idf
      dum  = fung
      idum = idg
      idum = iq0
      epsi = 1.D11

      return
      end                  

C     =====================================
      subroutine getids(idmin,idmax,nwords)
C     =====================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*80 subnam
      data subnam /'GETIDS ( IDMIN, IDMAX, NWORDS )'/

      call sqcErrMsg(subnam,
     +   'GETIDS obsolete, please use NWUSED( nwtot, nwuse, nwtab )')

      idum  = idmin               !avoid compiler warning
      idum  = idmax               !avoid compiler warning
      idum  = nwords              !avoid compiler warning

      return
      end

C     ===========================
      subroutine dumpwt(lun,file)
C     ===========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*(*) file

      character*80 subnam
      data subnam /'DUMPWT ( LUN, FILE )'/

      call sqcErrMsg(subnam,
     +   'DUMPWT obsolete, please use DMPWGT(itype,lun,file)')

      idum  = lun               !avoid compiler warning
      leng  = imb_lenoc(file)   !avoid compiler warning

      return
      end

C     ============================
      subroutine FastFac(id,funxq)
C     ============================

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FASTFAC ( ID, FUNXQ )'/

      call sqcErrMsg(subnam,
     + 'FASTFAC obsolete, please use FASTKIN')

      idum  = id              !avoid compiler warning
      ddum  = funxq           !avoid compiler warning

      return
      end

C     ======================================
      subroutine FtimesW(w,fun,jd1,id2,iadd)
C     ======================================

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FTIMESW ( W, FUN, ID1, ID2, IADD )'/

      call sqcErrMsg(subnam,
     + 'FTIMESW obsolete, please use WTIMESF')

      ddum  = w                !avoid compiler warning
      ddum  = fun              !avoid compiler warning
      idum  = jd1              !avoid compiler warning
      idum  = id2              !avoid compiler warning
      idum  = iadd             !avoid compiler warning

      return
      end
      
C     =================================================
      double precision function FcrossC(w,id,idf,ix,iq)
C     =================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*80 subnam
      data subnam /'FCROSSC ( W, IDW, IDF, IX, IQ )'/

      
      dum     = w
      idum    = id
      idum    = idf
      idum    = ix
      idum    = iq
      FcrossC = 0.D0

      call sqcErrMsg(subnam,'FCROSSC obsolete please use FCROSSK')

      return
      end      

C     ========================
      subroutine SetPdf(itype)
C     ========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SETPDF ( ITYPE )'/

      call sqcErrMsg(subnam,'SETPDF obsolete')

      idum  = itype            !avoid compiler warning

      return
      end
      
C     ========================
      subroutine GetPdf(itype)
C     ========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'GETPDF ( ITYPE )'/

      call sqcErrMsg(subnam,'GETPDF obsolete')

      idum  = itype            !avoid compiler warning

      return
      end      

C     ===========================
      subroutine FastPdf(id,coef)
C     ===========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FASTPDF ( ID, COEF )'/

      call sqcErrMsg(subnam,
     +   'FASTPDF obsolete, please use FASTSUM( iset, coef, id )')

      idum  = id               !avoid compiler warning
      ddum  = coef             !avoid compiler warning

      return
      end

C     ==================================
      subroutine FastFcC(w,idwt,id1,id2)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FASTFCC ( W, ID, ID1, ID2 )'/

      call sqcErrMsg(subnam,
     +   'FASTFCC obsolete, please use FASTFXK')

      ddum  = w          !avoid compiler warning
      idum  = idwt       !avoid compiler warning
      idum  = id1       !avoid compiler warning
      idum  = id2       !avoid compiler warning

      return
      end

C     ===========================
      subroutine FastAdd(id1,id2)
C     ===========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FASTADD ( ID1, ID2 )'/

      call sqcErrMsg(subnam,
     +   'FASTADD obsolete, please use FASTCPY ( id1, id2, iadd )')

      idum  = id1       !avoid compiler warning
      idum  = id2       !avoid compiler warning

      return
      end

C     ========================
      subroutine setabq(aq,bq)
C     ========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SETABQ ( AQ, BQ )'/

      call sqcErrMsg(subnam,'SETABQ obsolete')

      dum = aq
      dum = bq
      
      return
      end
      
C     ========================
      subroutine getabq(aq,bq)
C     ========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'GETABQ ( AQ, BQ )'/

      call sqcErrMsg(subnam,'GETABQ obsolete')

      dum = aq
      dum = bq
      
      return
      end

C     ===================================
      subroutine setthr(nfix,q2c,q2b,q2t)
C     ===================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SETTHR ( NFIX, Q2C, Q2B, Q2T )'/

      call sqcErrMsg(subnam,
     +              'SETTHR obsolete, pls use SETCBT instead')

      idum = nfix
      dum  = q2c
      dum  = q2t
      dum  = q2b
      
      return
      end
      
C     ==========================================
      subroutine getthr(nfix,q2c,q2b,q2t,intern)
C     ==========================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'GETTHR ( NFIX, Q2C, Q2B, Q2T, INTERN )'/

      call sqcErrMsg(subnam,
     +              'GETTHR obsolete, pls use GETCBT instead')

      idum = nfix
      dum  = q2c
      dum  = q2b
      dum  = q2t
      idum = intern

      return
      end

C     ================================
      subroutine pdfsum(id,wt,n,idout)
C     ================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PDFSUM ( IDIN, WGT, N, IDOUT )'/

      call sqcErrMsg(subnam,'PDFSUM obsolete')

      idum  = id
      dum   = wt
      idum  = n
      idout = 999 

      return
      end
      
C     ================================================
      double precision function pdfval(idf,xx,qq,mode)
C     ================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PDFVAL ( ID, X, Q2, MODE )'/

      pdfval = 0.D0
      idum   = idf
      dum    = xx
      dum    = qq
      idum   = mode

      call sqcErrMsg(subnam,'PDFVAL obsolete')

      return
      end
            
C     ============================================
      double precision function pgluon(xx,qq,jchk)
C     ============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PGLUON ( X, QMU2, ICHK )'/

      pgluon = 0.D0
      dum    = xx
      dum    = qq
      idfum  = jchk
   
      call sqcErrMsg(subnam,'PGLUON obsolete')

      return
      end

C     ================================================
      double precision function onepdf(xx,qq,def,jchk)
C     ================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'ONEPDF ( X, QMU2, DEF, ICHK )'/

      onepdf = 0.D0
      dum    = xx
      dum    = qq
      dum    = def
      jchk   = 999

      call sqcErrMsg(subnam,'ONEPDF obsolete')

      return
      end

C     ==================================
      subroutine allpdf(xx,qq,pdf,imode)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'ALLPDF ( X, QMU2, VAL, IMODE )'/

      dum   = xx
      dum   = qq
      dum   =  pdf
      idum  = imode

      call sqcErrMsg(subnam,'ALLPDF obsolete please call FPDFXQ')

      return
      end
      
C     ==================================
      subroutine pdfsxq(xx,qq,pdf,jchk)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PDFSXQ ( X, QMU2, PDF, ICHK )'/

      dum   = xx
      dum   = qq
      dum   = pdf
      idum  = jchk

      call sqcErrMsg(subnam,'PDFSXQ obsolete please call FPDFXQ')

      return
      end
      
C     ==================================
      subroutine pdfsij(ix,iq,pdf,jchk)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PDFSIJ ( IX, IQ, PDF, ICHK )'/

      idum  = ix
      idum  = iq
      dum   = pdf
      idum  = jchk

      call sqcErrMsg(subnam,'PDFSIJ obsolete please call FPDFIJ')

      return
      end                   

C     ===============================================
      double precision function sgnsxq(id,xx,qq,jchk)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SGNSXQ ( ID, X, QMU2, ICHK )'/

      sgnsxq = 0.D0
      idum   = id
      dum    = xx
      dum    = qq
      idum   = jchk

      call sqcErrMsg(subnam,'SGNSXQ obsolete please call FSNSXQ')

      return
      end
       
C     ===============================================
      double precision function sgnsij(id,ix,iq,jchk)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SGNSIJ ( ID, IX, IQ, ICHK )'/

      sgnsij = 0.D0
      idum   = id
      idum   = ix
      idum   = iq
      idum   = jchk

      call sqcErrMsg(subnam,'SGNSIJ obsolete please call FSNSIJ')

      return
      end                    

C     ===============================================
      double precision function pfunxq(id,xx,qq,jchk)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PFUNXQ ( ID, X, QMU2, ICHK )'/

      pfunxq = 0.D0
      idum   = id
      dum    = xx
      dum    = qq
      idum   = jchk

      call sqcErrMsg(subnam,'PFUNXQ obsolete please call FVALXQ')

      return
      end
      
C     ===============================================
      double precision function pfunij(id,ix,iq,jchk)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PFUNIJ ( ID, IX, IQ, ICHK )'/

      pfunij = 0.D0
      idum   = id
      idum   = ix
      idum   = iq
      idum   = jchk

      call sqcErrMsg(subnam,'PFUNIJ obsolete please call FVALIJ')

      return
      end         

C     ================================================
      double precision function psumxq(def,xx,qq,jchk)
C     ================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PSUMXQ ( DEF, X, QMU2, ICHK )'/

      psumxq = 0.D0
      dum    = def
      dum    = xx
      dum    = qq
      idum   = jchk

      call sqcErrMsg(subnam,'PSUMXQ obsolete please call FSUMXQ')

      return
      end
                        
C     ================================================
      double precision function psumij(def,ix,iq,jchk)
C     ================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PSUMIJ ( DEF, IX, IQ, ICHK )'/

      psumij = 0.D0
      dum    = def
      idum   = ix
      idum   = iq
      idum   = jchk

      call sqcErrMsg(subnam,'PSUMIJ obsolete please call FSUMIJ')

      return
      end        

C     ============================================
      subroutine stfval(istf,ityp,id,x,q,f,n,mode)
C     ============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'STFVAL ( ISTF, ITYP, ID, X, Q2, F, n, MODE )'/

      idum = istf
      idum = ityp
      idum = id
      dum  = x
      dum  = q
      dum  = f
      idum = n
      idum = mode

      call sqcErrMsg(subnam,'STFVal obsolete, pls use ZMSTF package')

      return
      end

C     =========================================    
      subroutine StrFun(istf,user,x,q,f,n,mode)
C     =========================================

C--   Obsolete

      implicit double precision (a-h,o-z)
      
      include 'qluns1.inc'

      character*80 subnam
      data subnam /'STRFUN ( ISTF, DEF, X, Q2, F, N, MODE )'/

      idum = istf
      dum  = user
      dum  = x
      dum  = q
      dum  = f
      idum = n
      idum = mode
      
      call sqcErrMsg(subnam,'STRFUN obsolete, pls use ZMSTF package')
      
      return
      end

C     ==========================
      subroutine StampIt(string)
C     ==========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*(*) string

      character*80 subnam
      data subnam /'STAMPIT ( STRING )'/
      len = imb_lenoc(string)   !avoid compiler warning
      call sqcErrMsg(subnam,'STAMPIT obsolete')

      return
      end

C     ================================
      integer function idSplij(string)
C     ================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*80 subnam
      data subnam /'IDSPLIJ ( STRING )'/
      
      character*(*) string

      len     = imb_lenoc(string)   !avoid compiler warning
      idSplij = 0

      call sqcErrMsg(subnam,'IDSPLIJ obsolete please use IDSPFUN')

      return
      end
      


      
      
