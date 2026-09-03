
C--   This is the file srcPdfNN.f containing the nxn pdf output routines

C--   double precision function dqcEvPdfij(ww,id,iy,it)
C--   subroutine sqcInterpList(subnam,w,igl,x,q,f,n,jchk)
C--   subroutine sqcEvPCopy(w,jdf,def,n,id0,nfmax)
C--   subroutine sqcEvTable(w,idf,xx,nxx,qq,nqq,fff)

C     =================================================
      double precision function dqcEvPdfij(ww,id,iy,it)
C     =================================================

C--   Get pdf at gridpoint
C--
C--   ww     (in)  store (not aligned)
C--   id     (in)  pdf table identifier in global format
C--   iy     (in)  y-grid index
C--   it     (in)  t-grid index (same as iq) can be < 0

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension ww(*)

      iz         = izfit5(it)
      ia         = iqcG5ijk(ww,iy,iz,id)
      dqcEvPdfij = ww(ia)

      return
      end

c     ===================================================
      subroutine sqcInterpList(subnam,w,igl,x,q,f,n,jchk)
C     ===================================================

C--   Interpolate arbitrary long list of pdfs (called by EVPLIST)
C--   Uses interpolation routine of fast convolution engine
C--
C--   submnam   (in) : subroutine name (for error message)
C--   w         (in) : workspace
C--   igl       (in) : global pdf identifier
C--   x         (in) : list of x points
C--   q         (in) : list of Q2 points
C--   f        (out) : list of interpolated pdfs
C--   n         (in) : number of items in x, q and f
C--   jchk      (in) : 1 = check all (x,q) inside grid or cuts

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*),x(*),q(*),f(*)

      character*(*) subnam

      dimension xx(mpt0),qq(mpt0)

C--   Fill output array f in batches of mpt0 words
      margin = 0
      ipt    = 0
      jj     = 0
      do i = 1,n
        ipt     = ipt+1
        xx(ipt) = x(i)
        qq(ipt) = q(i)
        if(ipt.eq.mpt0) then
C--       Setup interpolation mesh
          call sqcSetMark(xx,qq,mpt0,margin,ierr)
C--       At least one x,qmu2 outside grid
          if(jchk.eq.1 .and. ierr.eq.1) then
            call sqcErrMsg(subnam,'At least one x, mu2 outside cuts')
          endif
          call sqcFastFxq(w,igl,f(jj*mpt0+1),mpt0) !interpolate
          ipt = 0
          jj  = jj+1
        endif
      enddo
C--   Flush remaining ipt points
      if(ipt.ne.0) then
C--     Setup interpolation mesh
        call sqcSetMark(xx,qq,ipt,margin,ierr)
C--     At least one x,qmu2 outside grid
        if(jchk.eq.1 .and. ierr.eq.1) then
          call sqcErrMsg(subnam,'At least one x, mu2 outside cuts')
        endif
        call sqcFastFxq(w,igl,f(jj*mpt0+1),ipt) !interpolate
      endif

      return
      end

C     =================================================
      subroutine sqcEvPCopy(w,jdf,def,n,id0,nfmax,ierr)
C     =================================================

C--   Copy a set of pdfs in w to QCDNUM internal memory

C--   w               (in) : local workspace
C--   jdf(0:12+n)     (in) : pdf identifiers in w in global format
C--   def(-6:6,12)    (in) : contribution of flavour i to jdf(j)
C--   n               (in) : number of extra tables beyond 13
C--   id0             (in) : gluon id in stor7 (global format)
C--   nfmax           (in) : max number of flavours
C--   ierr           (out) :  0 = OK
C--                          >0 = cannot invert (submatrix of) def 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      dimension w(*)
      dimension itypes(6)
      dimension jdf(0:12+n), def(-6:6,12), dinv(12,12,3:6)
      dimension coef00(3:6), coefij(3:6)

C--   Initialise
      ierr = 0
      call smb_IFill(itypes,6,0)
      call smb_VFill(coef00,4,0.D0)

C--   dinv(i,j,nf) is the contribution of pdf(j) in w to pdf(i) in stor7
      do nf = 3,nfmax
        call sqcGetDinv(def,dinv(1,1,nf),nf,jerr)
        if(jerr.ne.0) then
          ierr = nf
          return
        endif
      enddo

C--   Put jdf(0)=gluon from w into stor7
      caLL sqcPdfCopy(w,jdf(0),stor7,id0,coef00,0)  !copy gluon to stor7

C--   Accumulate weighted sum cj*jdf(j) from w into stor7
C--   Copy eplus
      do i = 1,nfmax                   !loop over quark pdf i in stor7
        idi = id0+i                    !get global idi in stor7
        do j = 1,2*nfmax               !loop over quark pdf j in w
          do nf = 3,nfmax
            coefij(nf) = dinv(i,j,nf)  !contribution of pdf j to pdf i
          enddo
          iadd = min(j-1,1)            !iadd=0 for j=1, iadd=1 otherwise
          call sqcPdfCopy(w,jdf(j),stor7,idi,coefij,iadd)  !weighted sum
        enddo
      enddo
C--   Set rest of eplus to singlet
      id1 = id0+1                      !get global id1 in stor7
      do i = nfmax+1,6
        idi = id0+i                    !get global idi in stor7
        call sqcPdfCopy(stor7,id1,stor7,idi,coef00,0)   !copy id1 to idi
      enddo
C--   Copy eminus
      do i = 7,6+nfmax                 !loop over quark pdf i in stor7
        idi = id0+i                    !get global idi in stor7
        do j = 1,2*nfmax               !loop over quark pdf j in w
          do nf = 3,nfmax
            coefij(nf) = dinv(i,j,nf)  !contribution of pdf j to pdf i
          enddo
          iadd = min(j-1,1)            !iadd=0 for j=1, iadd=1 otherwise
          call sqcPdfCopy(w,jdf(j),stor7,idi,coefij,iadd)  !weighted sum
        enddo
      enddo
C--   Set rest of eminus to valence
      id7 = id0+7                      !get global id7 in stor7
      do i = 7+nfmax,12
        idi = id0+i                    !get global idi in stor7
        call sqcPdfCopy(stor7,id7,stor7,idi,coef00,0)   !copy id7 to idi
      enddo
C--   Now copy the n extra tables
      do i = 13,12+n
        idi = id0+i                    !get global idi in stor7
        call sqcPdfCopy(w,jdf(i),stor7,idi,coef00,0)  !copy pdf j to idi
      enddo

      return
      end

C     ==============================================
      subroutine sqcEvTable(w,idf,xx,nxx,qq,nqq,fff)
C     ==============================================

C--   Return pdfs interpolated on a x-mu2 grid 
C--
C--   w            (in)   local workspace
C--   idf          (in)   pdf identifier in global format
C--   xx(i)        (in)   table of x-values
C--   nxx          (in)   number of x-values
C--   qq(i)        (in)   table of mu2 values
C--   nqq          (in)   number of mu2 values
C--   fff(nxx*nqq) (out)  interpolated pdfs (linear store)
C--
C--   Remark: xx(i) and qq(i) are all supposed to be within the grid
C--   Remark: result fff(nx*nq) is stored linearly 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension w(*), xx(*), qq(*), fff(*)
      dimension wy(6),wz(6)

C--   Loop over x and q and interpolate idf, store result in fff
      margin = 0
      ii     = 0
      do iq = 1,nqq
        t = log(qq(iq))
        do ix = 1,nxx
          ii   = ii + 1
          y    = -log(xx(ix))
          call sqcZmesh(y,t,margin,iy1,iy2,iz1,iz2,it1)
          ny   = iy2-iy1+1
          nz   = iz2-iz1+1
          iag  = iqcG5ijk(w,iy1,iz1,idf)
          call sqcIntWgt(iy1,ny,it1,nz,y,t,wy,wz)
          fff(ii) = dqcPdfPol(w,iag,ny,nz,wy,wz)
        enddo
      enddo
      
      return
      end
