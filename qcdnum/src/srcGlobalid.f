
C--   This is the file srdGlobalId with conversions to global addressing

C--   integer function iqcIdAtab(id,iset)
C--   integer function iqcIaAtab(iz,id,iset)

C--   integer function iqcIdStab(id,iset)
C--   integer function iqcIaStab(iz,id,iset)

C--   integer function iqcPdfIjkl(iy,it,id,iset)
C--   integer function iqcIdPdfLtoG(iset,id)
C--   subroutine sqcIdPdfGtoL(idgin,isetout,idout)

C--   integer function iqcIbufGlobal(ibuf)
C--   subroutine sqcInvalidateBuf(ibuf)


C===================================================================
C==== Alpha_s table addressing =====================================
C===================================================================

C     ===================================
      integer function iqcIdAtab(id,iset)
C     ===================================

C--   Get global identifier of alfas table in stor7
C--
C--   id   (in) :  alfas table index [-mord0,mord0]
C--   iset (in) : -1 not used here (fast tables)
C--                0 base tables
C--               >0 pdf set id [1,mset0]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      if(iset.lt.0 .or. iset.gt.mset0)  stop 'iqcIdAtab wrong iset'
      if(id.lt.-mord0 .or. id.gt.mord0) stop 'iqcIdAtab wrong id'

      iqcIdAtab = 1000*isetf7(iset) + 601 + id + mord0

      return
      end

C     ======================================
      integer function iqcIaAtab(iz,id,iset)
C     ======================================

C--   Get address of entry alfas(iz,id,iset) in stor7

C--   iz   (in) :  iz index    [1,nzz2]
C--   id   (in) :  table index [-mord0,mord0]
C--   iset (in) : -1 not used here (fast tables)
C--                0 base set
C--               >0 pdf set id [1,mset0]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      if(iz.lt.1 .or. iz.gt.nzz2) stop 'iqcIaAtab wrong iz'

      jd        = iqcIdAtab(id,iset)                      !get global id
      iqcIaAtab = iqcG6ij(stor7,iz,jd)                      !get address

      return
      end


C=======================================================================
C==   Subgrid pointer tables in the internal store =====================
C=======================================================================

C     ===================================
      integer function iqcIdStab(id,iset)
C     ===================================

C--   Get global identifier of subgrid pointer table in stor7
C--
C--   id   (in) :  subgrid pointer table index [1,nsubt0]
C--   iset (in) : -1 not used here (fast tables)
C--                0 base tables
C--               >0 pdf set id [1,mset0]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      if(iset.lt.0 .or. iset.gt.mset0)  stop 'iqcIdStab wrong iset'
      if(id.lt.1 .or. id.gt.nsubt0)     stop 'iqcIdStab wrong id'

      iqcIdStab = 1000*isetf7(iset) + 700 + id

      return
      end

C     ======================================
      integer function iqcIaStab(iz,id,iset)
C     ======================================

C--   Get address of entry isubgrid(iz,id,iset) in stor7

C--   iz   (in) :  iz index    [1,nzz2]
C--   id   (in) :  subgrid pointer table index [1,nsubt0]
C--   iset (in) : -1 not used here (fast tables)
C--                0 base set
C--               >0 pdf set id [1,mset0]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      if(iz.lt.0 .or. iz.gt.nzz2) stop 'iqcIaStab wrong iz'

      jd        = iqcIdStab(id,iset)                      !get global id
      iqcIaStab = iqcG7ij(stor7,iz,jd)                      !get address

      return
      end

C=======================================================================
C==   Pdf address in the internal store ================================
C=======================================================================

C     ==========================================
      integer function iqcPdfIjkl(iy,it,id,iset)
C     ==========================================

C--   Return address of pdf in internal memory
C--   iy   (in) :  y index
C--   it   (in) :  t index
C--   id   (in) :  pdf table identifier [ifrst7,ilast7] usually [0,12]
C--   iset (in) : -1 fast tables
C--                0 base tables
C--               >0 pdf set id [1,mset0]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      if(iset.lt.-1 .or. iset.gt.mset0) stop 'iqcPdfIjk wrong iset'

      if(id.ge.0) then
        if(id.lt.ifrst7(iset) .or. id.gt.ilast7(iset)) then
          write(6,*) 'iqcPdfIjk wrong id = ',id
          stop
        endif
        jd = 1000*isetf7(iset) + id - ifrst7(iset) + 501
      else
        if(-id.lt.ifrst7(0) .or. -id.gt.ilast7(0)) then
          write(6,*) 'iqcPdfIjk wrong id = ',id
          stop
        endif
        jd = 1000*isetf7(0) - id - ifrst7(0) + 501
      endif

      iqcPdfIjkl = iqcG5ijk(stor7,iy,it,jd)

      return
      end

C     ========================================
      integer function iqcIdPdfLtoG(iset,id)
C     ========================================

C--   Convert iset [0,mset0] and id [ifrst7,ilast7] to global pdf id
C--
C--   iset (in) : -1 fast tables
C--                0 base tables
C--               >0 pdf set id [1,mset0]
C--   id   (in) : pdf table identifier [ifrst7,ilast7] usually [0,12]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      if(iset.lt.-1 .or. iset.gt.mset0) then
        write(6,*) 'iqcIdPdfLtoG wrong iset = ',iset
        stop 'iqcIdPdfLtoG wrong iset'
      endif

      if(id.ge.0) then
        if(id.lt.ifrst7(iset) .or. id.gt.ilast7(iset)) then
         write(6,*) 'iqcIdPdfLtoG wrong id = ',id
         stop
        endif
        iqcIdPdfLtoG = 1000*isetf7(iset) + id - ifrst7(iset) + 501
      else
        if(-id.lt.ifrst7(0) .or. -id.gt.ilast7(0)) then
         write(6,*) 'iqcIdPdfLtoG wrong id = ',id
         stop
        endif
        iqcIdPdfLtoG = 1000*isetf7(0) - id - ifrst7(0) + 501
      endif

      return
      end

C     ============================================
      subroutine sqcIdPdfGtoL(idgin,isetout,idout)
C     ============================================

C--   Convert global pdf id to iset [-1,mset0] and id [ifrst7,ilast7]
C--   Returns -2 if isetout not found

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      jset = abs(idgin)/1000
      if(jset.lt.1 .or. jset.gt.mst0) stop 'sqcIdPdfGtoL wrong jset'
      do i = -1,mset0
        if(jset.eq.isetf7(i)) then
          isetout = i
          idout   = abs(idgin) - 1000*isetf7(i) + ifrst7(i) - 501
          if(isetout.lt.-1 .or. isetout.gt.mset0)
     +    stop 'sqcIdPdfGtoL wrong isetout'
          if(idout.lt.ifrst7(i) .or. idout.gt.ilast7(i))
     +    stop 'sqcIdPdfGtoL wrong idout'
          return
        endif
      enddo

      stop 'sqcIdPdfGtoL isetout not found'
      isetout = -2
      idout   = -2

      return
      end

C=======================================================================
C==  Fast buffer addressing ============================================
C=======================================================================

C     ====================================
      integer function iqcIbufGlobal(ibuf)
C     ====================================

C--   Get the global identifier of fast buffer ibuf in stor7
C--
C--   ibuf (in) : Buffer index [1,mbf0]
C--
C--   Note that scratch buffers have iset = 0

      implicit double precision (a-h,o-z)

      iqcIbufGlobal = iqcIdPdfLtoG(-1,ibuf)

      return
      end

C     =================================
      subroutine sqcInvalidateBuf(ibuf)
C     =================================

C--   Invalidate fast buffer ibuf
C--
C--   ibuf (in) : buffer index in local format [1,mbf0]
C--
C--   Remark: ibuf = 0 invalidates all fast buffers

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'

      if(ibuf.eq.0) then
        ibmi = ibuf
        ibma = ibuf
      else
        ibmi = 1
        ibma = mbf0
      endif

C--   Global identifiers
      igmi = iqcIbufGlobal(ibmi)
      igma = iqcIbufGlobal(ibma)

C--   Invalidate tables
      do idg = igmi,igma
        call sqcInvalidate(stor7,idg)
      enddo

      return
      end


