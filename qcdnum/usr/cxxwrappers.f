
C--   This is the file cppwrappers.f with wrappers that set character
C--   arguments to fixed length, needed for the C++ interface to QCDNUM
C--
C--   subroutine qcardsCPP(usub,filename,ls,iprint)
C--   subroutine qcbookCPP(action,l1,key,l2)

C     =============================================
      subroutine qcardsCPP(usub,filename,ls,iprint)
C     =============================================

      implicit double precision (a-h,o-z)

      external usub
      character*(100) filename

      if(ls.gt.100) stop 'qcardsCPP: input file name > 100 characters'

      call qcards(usub,filename(1:ls),ls,iprint)

      return
      end

C     ======================================
      subroutine qcbookCPP(action,l1,key,l2)
C     ======================================

      implicit double precision (a-h,o-z)

      character*(100) action, key

      if(l1.gt.100) stop 'qcbookCPP: input ACTION size > 100 characters'
      if(l2.gt.100) stop 'qcbookCPP: input KEY size > 100 characters'

      call qcbook(action(1:l1),key(1:l2))

      return
      end



