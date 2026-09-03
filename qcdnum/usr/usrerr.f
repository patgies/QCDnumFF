
C--   This is the file usrerr.f containing the error bits and messages

C--   subroutine sqcBitini
C--   subroutine sqcMakefl(subnam,ichk,iset,idel,iall)
C--   subroutine sqcChkIni(subnam)
C--   subroutine sqcChkflg(jset,ichk,subnam)
C--   subroutine sqcChekit(jset,ichk,jbit)
C--   subroutine sqcSetflg(iset,idel,iall)
C--   subroutine sqcSetbit(ibit,iword,n)
C--   subroutine sqcDelbit(ibit,iword,n)
C--   integer function iqcGetbit(ibit,iword,n)
C--
C--   subroutine sqcErrMsg(subnam,message)
C--   subroutine sqcErrMsg2(subnam,message1,message2)
C--   subroutine sqcPdfMsg(subnam,idnam,idin)
C--   subroutine sqcCutMsgD(subnam,vnam,cnam,dval,cut,noextra)
C--   subroutine sqcCutMsgI(subnam,vnam,cnam,ival,cut,noextra)
C--   subroutine sqcFstMsg(subnam)
C--   subroutine sqcSetMsg(subnam,setnam,iset)
C--   subroutine sqcParMsg(subnam,setnam,iset)
C--   subroutine sqcNtbMsg(subnam,setnam,iset)
C--   subroutine sqcMemMsg(subnam,nwords,ierr)
C--
C--   integer function iqcSjekId(subnam,idnam,w,idin,imi,ima,ifl,lint)
C--   subroutine sqcChkTyp12(subnam,cid1,cid2,id1,id2,mat12)
C--   subroutine sqcChkIoy12(subnam,cid1,cid2,w1,id1,w2,id2)
C--
C--   subroutine setUmsgCPP(name,ls)
C--   subroutine setUmsg(name)
C--   subroutine clrUmsg

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM error message routines from usrerr.f             */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM error message routines from usrerr.f               */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM error message routines from usrerr.f               */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C==   ==================================================================
C==   Manage status bits and error messages ============================
C==   ==================================================================
      
C     ====================
      subroutine sqcBitini
C     ====================

C--   Assign bits and define error messages
C--   Called by qcinit

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qibit4.inc'

C--   !!!!Warning!!!! increase dimension of errmsg3 if # bits > 10
      character*45    errmsg3,commsg3
      common /qemsg3/ errmsg3(10),commsg3(10)

C--                      '123456789012345678901234567890123456789012345'
      ibInit4          = 1
      errmsg3(ibInit4) = 'QCDNUM not initialized                       '
      commsg3(ibInit4) = 'Please call QCINIT                           '

      ibXgri4          = 2
      errmsg3(ibXgri4) = 'No x-grid available                          '
      commsg3(ibXgri4) = 'Please call GXMAKE                           '

      ibQgri4          = 3
      errmsg3(ibQgri4) = 'No Q2-grid available                         '
      commsg3(ibQgri4) = 'Please call GQMAKE                           '

      ibThrs4          = 4
      errmsg3(ibThrs4) = 'No thresholds defined                        '
      commsg3(ibThrs4) = 'Please call SETCBT or MIXFNS                 '

      ibPdfs4          = 5
      errmsg3(ibPdfs4) = 'Pdf set not filled                           '
      commsg3(ibPdfs4) = 'Please call EVOLFG, PDFCPY, EXTPDF or EVPCOPY'
      
      ibFbuf4          = 6
      errmsg3(ibFbuf4) = 'No fast buffers avaialble                    '
      commsg3(ibFbuf4) = 'Please call FASTINI                          '
      
C--   NB: if new bits are defined like ibXXX4 then do not forget
C--       to declare them in qibit4.inc !!!!       

      return
      end

C     ===========================================
      subroutine sqcMakefl(subnam,ichk,iset,idel)
C     ===========================================

C--   Set bitpatterns for user subroutines (large switchyard)
C--   Called at the first invocation of each user routine
C--
C--   Input:  subnam subroutine name
C--
C--   Output: ichk(mbp0)  pattern of check-bits (check memory bit) 
C--           iset(mbp0)  pattern of set-bits   (set   memory bit)
C--           idel(mbp0)  pattern of del-bits   (unset memory bit)
C--
C--   Common: The bit assignments are defined in sqcBitini and stored
C--           in the common block /qibit4/    

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'
      
      character*(*)   subnam

      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

C--   Check if qcdnum is initialized
      call sqcChkIni(subnam)

C--   Initialize
      do i = 1,mbp0
        ichk(i) = 0
        iset(i) = 0
        idel(i) = 0
      enddo

      if    (subnam(1:6) .eq. 'SETVAL') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'SETINT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'QSTORE') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'SETORD') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'SETALF') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'SETCBT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibThrs4,iset,mbp0)   !validate thresholds

      elseif(subnam(1:6) .eq. 'MIXFNS') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibThrs4,iset,mbp0)   !validate thresholds

      elseif(subnam(1:6) .eq. 'GETCBT') then
C--   ----------------------------------
        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

      elseif(subnam(1:7) .eq. 'NFRMIQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

      elseif(subnam(1:7) .eq. 'NFLAVS') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

      elseif(subnam(1:6) .eq. 'SETABR') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'CPYPAR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'GETPAR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'USEPAR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'KEYPAR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'KEYGRP') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'PUSHCP') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'PULLCP') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized


      elseif(subnam(1:6) .eq. 'GXMAKE') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

        call sqcSetbit(ibXgri4,iset,mbp0)   !validate x-grid

        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs
        call sqcSetbit(ibFbuf4,idel,mbp0)   !invalidate fast buffers

      elseif(subnam(1:6) .eq. 'IXFRMX') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:6) .eq. 'XFRMIX') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:6) .eq. 'XXATIX') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:6) .eq. 'GQMAKE') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

        call sqcSetbit(ibQgri4,iset,mbp0)   !validate q2-grid 
        
        call sqcSetbit(ibThrs4,idel,mbp0)   !invalidate thresholds
        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs
        call sqcSetbit(ibFbuf4,idel,mbp0)   !invalidate fast buffers

      elseif(subnam(1:6) .eq. 'IQFRMQ') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'QFRMIQ') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'QQATIQ') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'GRPARS') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'GXCOPY') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:6) .eq. 'GQCOPY') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'FILLWT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'FILLWC') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'DMPWGT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'READWT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'NWUSED') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
      elseif(subnam(1:6) .eq. 'SETLIM') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
      elseif(subnam(1:6) .eq. 'GETLIM') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
      elseif(subnam(1:6) .eq. 'ASFUNC') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

      elseif(subnam(1:6) .eq. 'EVOLAS') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

      elseif(subnam(1:7) .eq. 'ALTABN') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

      elseif(subnam(1:6) .eq. 'EVOLFG') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

        call sqcSetbit(ibPdfs4,iset,mbp0)   !validate pdfs

      elseif(subnam(1:6) .eq. 'EVSGNS') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

        call sqcSetbit(ibPdfs4,iset,mbp0)   !validate pdfs

      elseif(subnam(1:6) .eq. 'PDFCPY') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

        call sqcSetbit(ibPdfs4,iset,mbp0)   !validate pdfs

        
      elseif(subnam(1:6) .eq. 'EXTPDF') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

        call sqcSetbit(ibPdfs4,iset,mbp0)   !validate pdfs

      elseif(subnam(1:6) .eq. 'USRPDF') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

        call sqcSetbit(ibPdfs4,iset,mbp0)   !validate pdfs
        
      elseif(subnam(1:6) .eq. 'ALLFXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'ALLFIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'BVALXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'BVALIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FVALXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FVALIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'SUMFXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'SUMFIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FSPLNE') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist        
        
      elseif(subnam(1:6) .eq. 'SPLCHK') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist
        
      elseif(subnam(1:6) .eq. 'FFLIST') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist        
        
      elseif(subnam(1:6) .eq. 'FTABLE') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'QCARDS') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'QCBOOK') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'MAKETAB') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
      elseif(subnam(1:7) .eq. 'SETPARW') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'GETPARW') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'DUMPTAB') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'READTAB') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'MAKEWTA') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'MAKEWTB') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'MAKEWRS') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'MAKEWTD') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'MAKEWTX') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:7) .eq. 'SCALEWT') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'COPYWGT') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'WCROSSW') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'WTIMESF') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'EVFILLA') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

      elseif(subnam(1:7) .eq. 'EVGETAA') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'EVDGLAP') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist

      elseif(subnam(1:7) .eq. 'CPYPARW') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'USEPARW') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'KEYPARW') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'KEYGRPW') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'IDSCOPE') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'EVPDFIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'EVPLIST') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'EVTABLE') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'EVPCOPY') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibPdfs4,iset,mbp0)   !validate pdfs
                 
      elseif(subnam(1:7) .eq. 'FCROSSK') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'FCROSSF') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'EFROMQQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'QQFROME') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'STFUNXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'FASTINI') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibThrs4,ichk,mbp0)   !check thresholds exist
        
        call sqcSetbit(ibFbuf4,iset,mbp0)   !validate fast buffers
        
      elseif(subnam(1:7) .eq. 'FASTCLR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        
      elseif(subnam(1:7) .eq. 'FASTINP') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist

      elseif(subnam(1:7) .eq. 'FASTEPM') then
C--   -----------------------------------
        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist

      elseif(subnam(1:7) .eq. 'FASTSNS') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:7) .eq. 'FASTSUM') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist
        
      elseif(subnam(1:7) .eq. 'FASTFXK') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        
      elseif(subnam(1:7) .eq. 'FASTFXF') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist

      elseif(subnam(1:7) .eq. 'FASTKIN') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        
      elseif(subnam(1:7) .eq. 'FASTCPY') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist

      elseif(subnam(1:7) .eq. 'FASTFXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist

      else
C--   ----

        goto 510

      endif
C--   -----

      return

C--   Error messages

 510  continue
      write(lunerr1,'(/'' sqcMAKEFL: unknown subroutine '',A10,'//
     +        '  '' ---> STOP'')') subnam(1:7)
      stop

      end

C     ============================
      subroutine sqcChkIni(subnam)
C     ============================

C--   Check if QCDNUM is initialized

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      
      character*(*)   subnam

      if(iniflg4.ne.123456) then
C--     Error message
        leng = imb_lenoc(subnam)
        write(lundef1,'(/1X,70(''-''))')
        write(lundef1,'('' Error in '',A,'' ---> STOP'')')
     +      subnam(1:leng)
        write(lundef1,'( 1X,70(''-''))')
        write(lundef1,'(
     +    '' QCDNUM not initialized (no call to QCINIT)'')')
        stop
      endif

      return
      end

C     ======================================
      subroutine sqcChkflg(jset,ichk,subnam)
C     ======================================

C--   Check status bits; abort with error message if not OK.
C--   Called by each user subroutine.
C--
C--   Input:  ichk(mbp0)            pattern of status bits to be checked
C--
C--   Common: istat4(mbp0,mset0)    QCDNUM status passed through /qcsflg/
C--           subnam(1:7)           name of calling routine

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qibit4.inc'
      include 'qsflg4.inc'
      include 'qstor7.inc'

      character*(*)   subnam
      character*45    errmsg3,commsg3
      common /qemsg3/ errmsg3(10),commsg3(10)
      
      character*46 etxt1(5)
C                  1234567890123456789012345678901234567890123456
      data etxt1 /'ISET  X no weights available ',
     +            'ISET  X no weights available ',
     +            'ISET  X no weights available ',
     +            'ISET  X temporarily disabled ',
     +            'ISET  X cannot evolve external pdfs'/

      character*37 etxt2(5)
C                  1234567890123456789012345678901234567
      data etxt2 /'ISET  X does not exist ',
     +            'ISET  X does not exist ',
     +            'ISET  X does not exist ',
     +            'ISET  X does not exist ',
     +            'ISET  X does not exist '/

      dimension ichk(mbp0)

C--   Loop over status words and check status bits
      do i = 1,mbp0
        ierr = imb_test1(ichk(i),istat4(i,jset))
        if(ierr.ne.0) then
          jword = i
          goto  500
        endif
      enddo

      return

C--   Error messages

 500  continue

C--   First status bit which failed
      jbit = 0
      do i = 1,32
        if(imb_gbitn(ichk(jword),i)      .eq.1 .and.
     +     imb_gbitn(istat4(jword,jset),i).eq.0   ) then
           jbit = i
           goto 510
        endif
      enddo

 510  continue
C--   Error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,'('' Error in '',A,'' ---> STOP'')') 
     +     subnam(1:leng)
      write(lunerr1,'( 1X,70(''-''))')
      if(jbit.eq.0) then
        write(lunerr1,'('' No error message found'')')
      elseif(jbit.eq.ibPdfs4) then
        write(etxt2(min(jset,5))(5:7),'(I3)') jset
        write(lunerr1,'(1X,A37)') etxt2(min(jset,5))
        write(lunerr1,'(1X,A45)') commsg3((jword-1)*32+jbit)
      elseif(jbit.eq.ibInit4) then
        write(lundef1,'(1X,A45)') errmsg3((jword-1)*32+jbit)
        write(lundef1,'(1X,A45)') commsg3((jword-1)*32+jbit)
      elseif(jbit.eq.ibThrs4) then
        write(lundef1,'(1X,A45)') errmsg3((jword-1)*32+jbit)
        write(lundef1,'(1X,A45)') commsg3((jword-1)*32+jbit)
      else
        write(lunerr1,'(1X,A45)') errmsg3((jword-1)*32+jbit)
        write(lunerr1,'(1X,A45)') commsg3((jword-1)*32+jbit)
      endif
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,'(/''Error was detected in a call to '',A)')
     +       usrnam3(1:leng)
      endif

      stop

      end

C     ====================================
      subroutine sqcChekit(jset,ichk,jbit)
C     ====================================

C--   Check status bits as in sqcChkflg but no abort with error message.
C--   Called by each user function which should run quiet.
C--
C--   Input:  ichk(mbp0)            pattern of status bits to be checked
C--   Output: jbit                  status bit which failed: 0 = OK
C--
C--   Common: istat4(mbp0,mset0)    QCDNUM status passed through /qcsflg/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      include 'qstor7.inc'

      dimension ichk(mbp0)

C--   Initialize
      jbit = 0

C--   Loop over status words and check status bits
      do i = 1,mbp0
        jerr = imb_test1(ichk(i),istat4(i,jset))
        if(jerr.ne.0) then
          jword = i
          goto  500
        endif
      enddo

      return

 500  continue
C--   Search for first status bit which failed
      jbit = 0
      do i = 1,32
        if(imb_gbitn(ichk(jword),i)      .eq.1 .and.
     +     imb_gbitn(istat4(jword,jset),i).eq.0   ) then
           jbit = i
           goto 510
        endif
      enddo

 510  continue
      return
      end
      
C     ====================================
      subroutine sqcSetflg(iset,idel,ipdf)
C     ====================================

C--   Update QCDNUM status words.
C--   Called by each user subroutine.
C--
C--   Input:  iset(mbp0)           pattern of status bits to be set to 1
C--           idel(mbp0)           pattern of status bits to be set to 0
C--           ipdf = 0             update status word of all pdf sets
C--                # 0             update status word of ipdf set
C--
C--   Common: istat4(mbp0,mset0)   QCDNUM status passed through /qcsflg/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qsflg4.inc'
      include 'qstor7.inc'

      dimension iset(mbp0), idel(mbp0)
      
      if(ipdf.eq.0) then
        j1 = 1
        j2 = mset0
      else
        j1 = ipdf
        j2 = ipdf
      endif    

C--   Loop over status words
      do j = j1,j2
        do i = 1,mbp0
C--       Set bits to one
          istat4(i,j) = ior(iset(i),istat4(i,j))
C--       Set bits to zero
          istat4(i,j) = iand(not(idel(i)),istat4(i,j))
        enddo
      enddo

      return
      end

C==   ==================================================================
C==   Set and get status bits ==========================================
C==   ==================================================================

C     ==================================
      subroutine sqcSetbit(ibit,iword,n)
C     ==================================

C--   Set bit in sequence of n integer words to 1.
C--
C--   Input:   ibit     integer in the range [1,n*32]
C--            iword(n) array of n 32bit integers

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      dimension iword(n)

C--   Which word?
      nwd = (ibit-1)/32 + 1
C--   Check
      if(nwd.lt.1 .or. nwd.gt.n) goto 500

C--   Which bit?
      ibt = mod(ibit-1,32) + 1
C--   Check
      if(ibt.lt.1 .or. ibt.gt.32) goto 510

C--   Set bit
      call smb_sbit1(iword(nwd),ibt)

      return

C--   Error messages

 500  continue
      write(lunerr1,
     +    '(/'' sqcSETBIT: iwd .gt. maxwd '',2I15,'//
     +    '  '' ---> STOP'')') nwd, n
      write(lunerr1,*) ' Input ibit = ', ibit
      write(lunerr1,*) ' Input n    = ', n
      stop

 510  continue
      write(lunerr1,
     +    '(/'' sqcSETBIT: ibt not in range [1,32] '',I5,'//
     +    '  '' ---> STOP'')') ibt
      write(lunerr1,*) ' Input  ibit = ', ibit
      write(lunerr1,*) ' Input  n    = ', n
      write(lunerr1,*) ' Output ibt  = ', ibt
      stop

      end

C     ==================================
      subroutine sqcDelbit(ibit,iword,n)
C     ==================================

C--   Set bit in sequence of n integer words to 0.
C--
C--   Input:   ibit     integer in the range [1,n*32]
C--            iword(n) array of n 32bit integers

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      dimension iword(n)

C--   Which word?
      nwd = (ibit-1)/32 + 1
C--   Check
      if(nwd.lt.1 .or. nwd.gt.n) goto 500

C--   Which bit?
      ibt = mod(ibit-1,32) + 1
C--   Check
      if(ibt.lt.1 .or. ibt.gt.32) goto 510

C--   Set bit
      call smb_sbit0(iword(nwd),ibt)

      return

C--   Error messages

 500  continue
      write(lunerr1,
     +    '(/'' sqcDELBIT: iwd .gt. maxwd '',2I5,'//
     +    '  '' ---> STOP'')') nwd, n
      stop

 510  continue
      write(lunerr1,
     +    '(/'' sqcDELBIT: ibt not in range [1,32] '',I5,'//
     +    '  '' ---> STOP'')') ibt
      stop

      end
 
C     ========================================
      integer function iqcGetbit(ibit,iword,n)
C     ========================================

C--   Get bit in sequence of n integer words.
C--
C--   Input:   ibit     integer in the range [1,n*32]
C--            iword(n) array of n 32bit integers

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      dimension iword(n)

C--   Which word?
      nwd = (ibit-1)/32 + 1
C--   Check
      if(nwd.lt.1 .or. nwd.gt.n) goto 500

C--   Which bit?
      ibt = mod(ibit-1,32) + 1
C--   Check
      if(ibt.lt.1 .or. ibt.gt.32) goto 510

C--   Value of bit
      iqcGetbit = imb_gbitn(iword(nwd),ibt)

      return

C--   Error messages

 500  continue
      write(lunerr1,
     +    '(/'' iqcGETBIT: iwd .gt. maxwd '',2I5,'//
     +    '  '' ---> STOP'')') nwd, n
      stop

 510  continue
      write(lunerr1,
     +    '(/'' iqcGETBIT: ibt not in range [1,32] '',I5,'//
     +    '  '' ---> STOP'')') ibt
      stop

      end

C==   ==================================================================
C==   Standard error messages  =========================================
C==   ==================================================================

C     ====================================
      subroutine sqcErrMsg(subnam,message)
C     ====================================

C--   Print a message

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*(*) message
      character*(*) subnam

C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) message
C--   Extra message      
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     ===============================================
      subroutine sqcErrMsg2(subnam,message1,message2)
C     ===============================================

C--   Print a 2-line message

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*(*) message1, message2
      character*(*) subnam

C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      leng = imb_lenoc(message1)
      write(lunerr1,*) message1(1:leng)
      leng = imb_lenoc(message2)
      write(lunerr1,*) message2(1:leng)
C--   Extra message
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     =======================================
      subroutine sqcPdfMsg(subnam,idnam,idin)
C     =======================================

C--   Handles error code from iPdfTab
C--
C--   idin  (in) : +1990XX if iset not in range   (ierr = 1)
C--                +2XX099 if id   not in range   (ierr = 2)
C--                +3XX0YY if iset does not exist (ierr = 3)
C--                +4XX0YY if id   does not exist (ierr = 4)
C--                +5XX0YY if id   is not filled  (ierr = 5)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      character*(*) subnam
      character*(*) idnam
      character*10  idtxt, istxt, iptxt, ms0txt, mp0txt
      character*80  etxt
      character*22  txt

C--               1234567890123456789012
      data txt  /' = IPDFTAB(ISET,ID) : '/

C--   Decode idin
      ierr = idin/100000
      idgl = idin - 100000*ierr
      iset = idgl/1000
      ipdf = idgl-1000*iset

C--   Convert to character string
      call smb_itoch(idgl,idtxt,ldtxt)
      call smb_itoch(iset,istxt,lstxt)
      call smb_itoch(ipdf,iptxt,lptxt)
      call smb_itoch(mset0,ms0txt,lms0txt)
      call smb_itoch(mpdf0,mp0txt,lmp0txt)

C--   Error messages
      if(ierr.eq.1) then
        write(etxt,
     +  '(A,A,''ISET not in range [1,'',A,'']'')')
     +  idnam, txt, ms0txt(1:lms0txt)
        call sqcErrMsg(subnam,etxt)
      elseif(ierr.eq.2) then
        write(etxt,
     +  '(A,A,''ID not in range [0,'',A,'']'')')
     +  idnam, txt, mp0txt(1:lmp0txt)
        call sqcErrMsg(subnam,etxt)
      elseif(ierr.eq.3) then
        write(etxt,
     +  '(A,A,''ISET = '',A,'' does not exist'')')
     +  idnam, txt, istxt(1:lstxt)
        call sqcErrMsg(subnam,etxt)
      elseif(ierr.eq.4) then
        write(etxt,
     +  '(A,A,''ID = '',A,'' does not exist in ISET = '',A)')
     +  idnam, txt, iptxt(1:lptxt), istxt(1:lstxt)
        call sqcErrMsg(subnam,etxt)
      elseif(ierr.eq.5) then
        write(etxt,
     +  '(A,A,''Pdf set '',A,'' is empty'')')
     +  idnam, txt, istxt(1:lstxt)
        call sqcErrMsg(subnam,etxt)
      else
        stop 'sqcPdfMsg: unknown error code'
      endif

      return
      end

C     ========================================================
      subroutine sqcCutMsgD(subnam,vnam,cnam,dval,cut,noextra)
C     ========================================================

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qsnam3.inc'
      
      character*(*) subnam, vnam, cnam

      leng = imb_lenoc(subnam)
      lenv = imb_lenoc(vnam)
      lenc = imb_lenoc(cnam)
      
C--   Now for the error message      

      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')

      write(lunerr1,
     +   '(1X,A,'' = '',1PE11.3,'' fails '',A,'' cut '',1PE11.3)')
     +     vnam(1:lenv),dval,cnam(1:lenc),cut
      
C--   Extra message      
      leng = imb_lenoc(usrnam3)            
      if(leng.gt.0 .and. noextra.ne.1) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     ========================================================
      subroutine sqcCutMsgI(subnam,vnam,cnam,ival,cut,noextra)
C     ========================================================

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qsnam3.inc'
      
      character*(*) subnam, vnam, cnam

      leng = imb_lenoc(subnam)
      lenv = imb_lenoc(vnam)
      lenc = imb_lenoc(cnam)
      
C--   Now for the error message      

      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')

      write(lunerr1,
     +   '(1X,A,'' = '',I5,'' fails '',A,'' cut '',1PE11.3)')
     +     vnam(1:lenv),ival,cnam(1:lenc),cut
      
C--   Extra message      
      leng = imb_lenoc(usrnam3)            
      if(leng.gt.0 .and. noextra.ne.1) then
        write(lunerr1,*) ' '
        write(lunerr1,*)
     +  ' Error was detected in a call to ',usrnam3(1:leng)
      endif

      stop
      end

C     ============================
      subroutine sqcFstMsg(subnam)
C     ============================

C--   Complain that fast engine scope has changed

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      
      character*80  subnam

      if(iscopeslot6.ne.iscopeuse6) then
        call sqcErrMsg(subnam,
     +  'Pdf scope was changed after FASTINI or FASTCLR(0)')
      endif

      return
      end

C     ========================================
      subroutine sqcSetMsg(subnam,setnam,iset)
C     ========================================

C--   Complain that iset does not exist

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      
      character*80  subnam
      character*(*) setnam
      character*80  emsg
      character*10  etxt

      call smb_itoch(iset,etxt,ltxt)
      write(emsg,'(A,'' = '',A,'' : nonexistent or empty pdf set'')')
     +  setnam, etxt(1:ltxt)
      call sqcErrMsg(subnam,emsg)

      return
      end

C     ========================================
      subroutine sqcParMsg(subnam,setnam,iset)
C     ========================================

C--   Check if iset is in scope
C--   It is not checked that iset actually exists

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'
      
      character*80  subnam
      character*(*) setnam
      character*80  emsg
      character*10  etxt
      character*52  bla
      character*56  bla1
C--                       1         2         3         4         5
C--              1234567890123456789012345678901234567890123456789012
      data bla /' : pdf set is not accepted by the convolution engine'/
      data bla1/
C--                1         2         3         4         5
C--       12345678901234567890123456789012345678901234567890123456
     +   'Call IDSCOPE to set the scope of pdf input to the engine'/

C--   Check or nocheck thats the question
      if(Lscopechek6) then
        key = ikeyf7(iset)
        if(key.ne.iscopekey6) then
          call smb_itoch(iset,etxt,ltxt)
          write(emsg,'(A,'' = '',A,A)')
     +    setnam, etxt(1:ltxt), bla
          call sqcErrMsg2(subnam,emsg,bla1)
        endif
      endif

      return
      end

C     ========================================
      subroutine sqcNtbMsg(subnam,setnam,iset)
C     ========================================

C--   Complain that iset cannot hold all tables

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'
      
      character*80  subnam
      character*(*) setnam
      character*80  emsg
      character*10  ntxt, stxt, ftxt, ltxt

      ntab = iqcGetNumberOfTables(stor7,isetf7(iset),5)
      ifst = ifrst7(iset)
      ilst = ilast7(iset)
      call smb_itoch(ntab,ntxt,lntxt)
      call smb_itoch(iset,stxt,lstxt)
      call smb_itoch(ifst,ftxt,lftxt)
      call smb_itoch(ilst,ltxt,lltxt)
      write(emsg,
     + '(A,'' = '',A,'' exists but can only hold '',A,
     + '' pdf tables with id = '',A,'' to '',A)')
     + setnam,stxt(1:lstxt),ntxt(1:lntxt),ftxt(1:lftxt),ltxt(1:lltxt)
      call sqcErrMsg(subnam,emsg)

      return
      end

C     ========================================
      subroutine sqcMemMsg(subnam,nwords,ierr)
C     ========================================

C--   Print message according to sqcMakeTab error code
C--
C--   subnam  (in)  name of the calling routine
C--   nwords  (in)  number of words used (<0 not enough space)
C--   ierr    (in)  >= 0  OK no message
C--                   -1  attempt to create set with 0 tables
C--                   -2  not enough space
C--                   -3  iset limit mst0 exceeded

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      
      character*80 subnam 
      character*80 emsg
      character*10 etxt

C--   Try to create set without tables
      if(ierr.eq.-1) then
        call sqcErrMsg(subnam,'Attempt to create set with no tables')
C--   Not enough space
      elseif(ierr.eq.-2) then
        nwfirst = abs(nwords)+1
        call smb_itoch(nwfirst,etxt,ltxt)
        write(emsg,'(''Need at least '',A,
     &  '' words in memory --> increase NWF0 '',
     &  ''in qcdnum.inc'')') etxt(1:ltxt) 
        call sqcErrMsg(subnam,emsg)
C--   Iset count exceeded
      elseif(ierr.eq.-3) then
        call smb_itoch(mst0,etxt,ltxt)
        write(emsg,'(''Exceed max '',A,
     &  '' table sets --> increase MST0'',
     &  '' in qcdnum.inc'')') etxt(1:ltxt)
         call sqcErrMsg(subnam,emsg)
C--   No error message
      else
        return
      endif

      return
      
      end

C==   ==================================================================
C==   Handle table identifiers =========================================
C==   ==================================================================

C     ================================================================
      integer function iqcSjekId(subnam,idnam,w,idin,imi,ima,ifl,lint)
C     ================================================================

C--   Returns id in global format for internal and local workspace
C--   Produces a fatal error if it fails
C--
C--   subnam   (in)  subroutine name for error message
C--   idnam    (in)  id name for error messsage
C--   w        (in)  workspace
C--   idin     (in)  input identifier (global format)
C--   imi(2)   (in)  min accepted values of isign, itype
C--   ima(2)   (in)  max accepted values of isign, itype
C--   ifl(1)   (in)  1 = pass zero identifier
C--   ifl(2)   (in)  1 = check id is filled
C--   lint    (out)  .true. if id refers to stor7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'pstor8.inc'

      logical       lint
      logical       lqcWpartitioned
      logical       lqcIdExists, lqcIsFilled
      character*(*) subnam, idnam
      dimension     w(*), imi(*), ima(*), ifl(*)

      character*80 emsg
      character*10 idtxt, istxt, ittxt

      iqcSjekId = 0             !avoid compiler warning

C--   Pass zero identifier if asked to do so
      if(idin.eq.0 .and. ifl(1).eq.1) then
        iqcSjekId = 0
        lint      = .false.
        return
      endif

C--   Catch errors in the call to ipdftab
      if(imi(2).eq.5 .and. ima(2).eq.5 .and.
     +  idin.gt.-600000 .and. idin .le.-100000) then
        call sqcPdfMsg(subnam,idnam,abs(idin))
        return
      endif

C--   Catch error in the call to idspfun
      if(idin.eq.-99999) then
        write(emsg,
     + '(A,'' = IDSPFUN(PIJ,IORD,ISET) : wrong idspfun input'')') idnam
       call sqcErrMsg(subnam,emsg)
       return
      endif

C--   Encode idin
      call smb_itoch(idin,idtxt,lidtxt)

C--   Check that idin is a global identifier
      if(abs(idin).lt.1000 .or. abs(idin).gt.99699) then
        write(emsg,'(A,'' = '',A,'' is not a global identifier'')')
     +  idnam, idtxt(1:lidtxt)
        call sqcErrMsg(subnam,emsg)
        return
       endif

C--   Decompose idin
      if(idin.ge.0) then
        jsgn =  1
        lint = .false.
      else
        jsgn = -1
        lint = .true.
      endif
      jset = abs(idin)/1000
      jtab = abs(idin)-1000*jset
      jtyp = jtab/100

C--   Encode jset and jtyp
      call smb_itoch(jset,istxt,listxt)
      call smb_itoch(jtyp,ittxt,littxt)

C--   Check sign
      if(jsgn.lt.imi(1) .or. jsgn.gt.ima(1)) then
        if(jsgn.eq. 1) then
          write(emsg,'(A,'' = '',A,'' : Workspace table not allowed'')')
     +    idnam, idtxt(1:lidtxt)
        endif
        if(jsgn.eq.-1) then
          write(emsg,'(A,'' = '',A,'' : Internal table not allowed'')')
     +    idnam, idtxt(1:lidtxt)
        endif
        call sqcErrMsg(subnam,emsg)
        return
      endif

C--   Check if W partitioned
      if(.not.lint .and. .not.lqcWpartitioned(w)) then
        call sqcErrMsg(subnam,'Workspace W is not partitioned')
        return
      endif

C--   Check jtype
      if(jtyp.lt.imi(2) .or. jtyp.gt.ima(2)) then
        write(emsg,
     +  '(A,'' = '',A,'' : Table type = '',A,'' is not allowed'')')
     +  idnam, idtxt(1:lidtxt), ittxt(1:littxt)
        call sqcErrMsg(subnam,emsg)
        return
      endif

C--   Check if table exists
      ierr = 0
      if(lint) then
         if(.not.lqcIdExists(stor7,abs(idin))) ierr = 1
      else
         if(.not.lqcIdExists(w,abs(idin)))     ierr = 1
      endif
      if(ierr.eq.1) then
        write(emsg,'(A,'' = '',A,'' : Table does not exist'')')
     +  idnam, idtxt(1:lidtxt)
        call sqcErrMsg(subnam,emsg)
        return
      endif

C--   Check if table is filled (when requested)
      ierr = 0
      if(ifl(2).eq.1) then
        if(lint) then
          if(.not.lqcIsFilled(stor7,abs(idin))) ierr = 1
        else
          if(.not.lqcIsFilled(w,abs(idin)))     ierr = 1
        endif
      endif

C--   Error message
      if(ierr.eq.1) then
        write(emsg,'(A,'' = '',A,'' : Table is empty'')')
     +    idnam, idtxt(1:lidtxt)
        call sqcErrMsg(subnam,emsg)
        return
      endif

C--   Evolution parameter check for pdf tables
      if(Lscopechek6 .and. jtyp.eq.5) then
        if(lint) then
          key  = int(dparGetPar(stor7,jset,idipver8))
        else
          key  = int(dparGetPar(  w  ,jset,idipver8))
        endif
        if(key.ne.iscopekey6) then
          write(emsg,
     + '(A,'' = '',A,'' : Table not accepted by convolution engine'')')
     +    idnam, idtxt(1:lidtxt)
          call sqcErrMsg2(subnam,emsg,
     +   'Call IDSCOPE to set the scope of pdf input to the engine')
        endif
      endif

C--   Return global identifier
      iqcSjekId = abs(idin)

      return
      end

C     ======================================================
      subroutine sqcChkTyp12(subnam,cid1,cid2,id1,id2,mat12)
C     ======================================================

C--   Error message if mat12(ityp1,ityp2) = 0
C--
C--   subnam     (in)   Subroutine name for error message
C--   cid1       (in)   Id1 name for error message
C--   cid2       (in)   Id2 name for error message
C--   id1        (in)   Global identifier of weight table, type 1-4
C--   id2        (in)   Global identifier of weight table, type 1-4
C--   mat12(6,6) (in)   Error message if mat12(ityp1,ityp2) = 0

      implicit double precision (a-h,o-z)

      dimension mat12(6,6)

      character*(*) subnam,cid1,cid2
      character*80  emsg
      character*10  nid1,nid2

      ityp1 = iqcGetLocalId(id1)/100
      if(ityp1.lt.1 .or. ityp1.gt.6) stop 'sqcChkTyp12 : wrong ityp1'
      ityp2 = iqcGetLocalId(id2)/100
      if(ityp2.lt.1 .or. ityp2.gt.6) stop 'sqcChkTyp12 : wrong ityp2'
      if(mat12(ityp1,ityp2) .eq. 0) then
        call smb_itoch(id1,nid1,lid1)
        call smb_itoch(id2,nid2,lid2)
        write(emsg,'(A,'' = '',A,'' '',A,'' = '',A,
     +  '' : incompatible table types'')')
     +  cid1, nid1(1:lid1), cid2, nid2(1:lid2)
        call sqcErrMsg(subnam,emsg)
      endif

      return
      end

C     ======================================================
      subroutine sqcChkIoy12(subnam,cid1,cid2,w1,id1,w2,id2)
C     ======================================================

C--   Error if id1 (input) is single and id2 (output) is double
C--
C--   subnam     (in)   Subroutine name for error message
C--   cid1       (in)   Id1 name for error message
C--   cid2       (in)   Id2 name for error message
C--   w1, id1    (in)   Global identifier in w1
C--   w2, id2    (in)   Global identifier in w2

      character*(*) subnam,cid1,cid2
      character*80  emsg
      character*10  nid1,nid2

      dimension w1(*), w2(*)

      logical lqcIsDouble, Ldouble1, Ldouble2

      Ldouble1 = lqcIsDouble(w1,id1)
      Ldouble2 = lqcIsDouble(w2,id2)

      if(.not.Ldouble1 .and. Ldouble2) then
        call smb_itoch(id1,nid1,lid1)
        call smb_itoch(id2,nid2,lid2)
        write(emsg,'(A,'' = '',A,'' '',A,'' = '',A,
     +  '' : wrong mix of splitting and coefficient function table'')')
     +  cid1, nid1(1:lid1), cid2, nid2(1:lid2)
        call sqcErrMsg(subnam,emsg)
      endif

      return
      end

C==   ===============================================================
C==   Set and clear package subroutine names ========================
C==   ===============================================================

C-----------------------------------------------------------------------
CXXHDR    void setumsg(string name);
C-----------------------------------------------------------------------
CXXHFW  #define fsetumsgcpp FC_FUNC(setumsgcpp,SETUMSGCPP)
CXXHFW    void fsetumsgcpp(char*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setumsg(string name)
CXXWRP  {
CXXWRP    int ls = name.size();
CXXWRP    char *cname = new char[ls+1];
CXXWRP    strcpy(cname,name.c_str());
CXXWRP    fsetumsgcpp(cname, &ls);
CXXWRP    delete[] cname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================
      subroutine setUmsgCPP(name,ls)
C     ==============================

      implicit double precision (a-h,o-z)

      character*(100) name

      if(ls.gt.100) stop 'setUmsgCPP: input NAME size > 100 characters'

      call setUmsg(name(1:ls))

      return
      end
            
C     ========================
      subroutine setUmsg(name)
C     ========================

C--   Set name of user subroutine

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      
      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)

      character*(*) name
      
      character*80 subnam
      data subnam /'SETUMSG ( CHNAME )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Fill usrnam3 provided it is empty
      if(imb_lenoc(usrnam3).eq.0) then
        call smb_cfill(' ',usrnam3)
        len = min(imb_lenoc(name),80)
        usrnam3(1:len) = name(1:len)
      endif
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void clrumsg();
C-----------------------------------------------------------------------
CXXHFW  #define fclrumsg FC_FUNC(clrumsg,CLRUMSG)
CXXHFW    void fclrumsg();
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void clrumsg()
CXXWRP  {
CXXWRP    fclrumsg();
CXXWRP  }
C-----------------------------------------------------------------------
            
C     ==================
      subroutine clrUmsg
C     ==================

C--   Clear name of user subroutine

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      
      logical first
      save    first
      data    first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first)
      
      character*80 subnam
      data subnam /'CLRUMSG'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

      call smb_cfill(' ',usrnam3)

      return
      end


