
C--  This is sqcPdfFill.f  with filling routines other than EVOLFG

C--   subroutine sqcPdfCpy(k7set1,k7set2)
C--   subroutine sqcExtPdf(func,ig0,n,del,nfheavy)
C--   subroutine sqcUsrPdf(func,ig0,n,del,nfheavy)

C     ===================================
      subroutine sqcPdfCpy(k7set1,k7set2)
C     ===================================

C--   Copy jset1 [1-mset0] to jset2 [1-mset0]

C--   k7set1 (in)    : source set id in stor7
C--   k7set2 (in)    : target set id in stor7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

C--   Copy parameter list
      call sparParAtoB(stor7,k7set1,stor7,k7set2)
C--   Copy type 5 + satellite tables, set iset pointer and validate
      idg1 = 1000*k7set1+500
      idg2 = 1000*k7set2+500
      ntab = iqcGetNumberOfTables(stor7,k7set1,5)
      do i = 1,ntab
        call sqcCopyType5(stor7,idg1+i,stor7,idg2+i)
        call sqcValidate(stor7,idg2+i)
      enddo
      ntb2 = iqcGetNumberOfTables(stor7,k7set2,5)
      do i = ntab+1,ntb2
        call sqcInvalidate(stor7,idg2+i)
      enddo
C--   Copy type 6
      idg1 = 1000*k7set1+600
      idg2 = 1000*k7set2+600
      ntab = iqcGetNumberOfTables(stor7,k7set1,6)
      do i = 1,ntab
        call sqcCopyType6(stor7,idg1+i,stor7,idg2+i)
      enddo
C--   Copy type 7
      idg1 = 1000*k7set1+700
      idg2 = 1000*k7set2+700
      ntab = iqcGetNumberOfTables(stor7,k7set1,7)
      do i = 1,ntab
        call sqcCopyType7(stor7,idg1+i,stor7,idg2+i)
      enddo

      return
      end

C     ============================================
      subroutine sqcExtPdf(func,ig0,n,del,nfheavy)
C     ============================================

C--   Fill set of basis pdf tables.

C--   func     (external): user supplied function giving gluon, qqbar
C--   ig0      (in)      : ig0 global id of gluon table in stor7
C--   n        (in)      : number of extra tables beyond gluon, qqbar
C--   del      (in)      : threshold offset mu2h*(1+-del)
C--   nfheavy (out)      : max nf encountered

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension qqbar(-6:mpdf0),epm(0:12)
      
      external func

C--   Dummy call to func
      dummy = func(0,1.D-1,1.D1,.true.)
C--   Global identifiers
      ig1  = ig0+1
      ig2  = ig0+2
      ig12 = ig0+12
C--   Address increment
      inc = iqcG5ijk(stor7,1,1,ig2)-iqcG5ijk(stor7,1,1,ig1)
C--   Loop over iz and it
      do iz = izmic5,izmac5
        it  = itfiz5( iz)
        nf  = itfiz5(-iz)
        qi  = exp(tgrid2(it))
C--     Detect threshold crossing: isign =  ... 0, 0, 0, -1, +1, 0, 0 ... 
        if(iz.ne.1 .and. iz.ne.nzz5) then
          nfizm = itfiz5(-(iz-1))
          nfizp = itfiz5(-(iz+1))
          isign = 2*nf-nfizm-nfizp
        else
          isign = 0
        endif  
        qi = qi * (1.D0+isign*del)        
        do iy = 1,iymac2
          xi = exp(-ygrid2(iy))
C--       Get qbar, gluon, q, extra
          do i = -6,6+n
            qqbar(i) = func(i,xi,qi,.false.)
          enddo
C--       Convert to epm
          call sqcQQBtoEPM(qqbar(-6),epm(0),nf)
C--       Fill e+- tables
          ia  = iqcG5ijk(stor7,iy,iz,ig0)-inc
          do i = 0,12
            ia        = ia+inc
            stor7(ia) = epm(i)
          enddo
C--       Fill extra tables
          ia  = iqcG5ijk(stor7,iy,iz,ig12)
          do i = 1,n
            ia        = ia+inc
            stor7(ia) = qqbar(6+i)
          enddo
C--     End of loop over y
        enddo  
C--   End of loop over t
      enddo

      nfheavy  = itfiz5(-izmac5)

      return
      end

C     ============================================
      subroutine sqcUsrPdf(func,ig0,n,del,nfheavy)
C     ============================================

C--   Fill set of user pdf tables.

C--   func     (external): user supplied function giving gluon, userpdf
C--   ig0      (in)      : ig0 global id of gluon table in stor7
C--   n        (in)      : number of tables beyond gluon
C--   del      (in)      : threshold offset mu2h*(1+-del)
C--   nfheavy (out)      : max nf encountered

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'


      external func

C--   Dummy call to func
      dummy = func(0,1.D-1,1.D1,.true.)
C--   Global identifiers
      ig1  = ig0+1
      ig2  = ig0+2
C--   Address increment
      inc = iqcG5ijk(stor7,1,1,ig2)-iqcG5ijk(stor7,1,1,ig1)
C--   Loop over iz and it
      do iz = izmic5,izmac5
        it  = itfiz5( iz)
        nf  = itfiz5(-iz)
        qi  = exp(tgrid2(it))
C--     Detect threshold crossing: isign =  ... 0, 0, 0, -1, +1, 0, 0 ...
        if(iz.ne.1 .and. iz.ne.nzz5) then
          nfizm = itfiz5(-(iz-1))
          nfizp = itfiz5(-(iz+1))
          isign = 2*nf-nfizm-nfizp
        else
          isign = 0
        endif
        qi = qi * (1.D0+isign*del)
        do iy = 1,iymac2
          xi = exp(-ygrid2(iy))
C--       Fill tables
          ia  = iqcG5ijk(stor7,iy,iz,ig0)-inc
          do i = 0,n
            ia        = ia+inc
            stor7(ia) = func(i,xi,qi,.false.)
          enddo
C--     End of loop over y
        enddo
C--   End of loop over t
      enddo

      nfheavy  = itfiz5(-izmac5)

      return
      end
