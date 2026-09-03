
C--   This is the file srcPdfInOut.f containing pdfs in/out routines

C--   double precision function dqcOneQpm(id,ia,nf,nfmax)
C--   subroutine sqcAllQpm(ia,nf,nfmax,qout)
C--   subroutine sqcEfromQQ(qvec,evec,nf,nfmax)
C--   subroutine sqcElistQQ(qvec,wt,id,n,nf,nfmax)
C--   subroutine sqcElistFF(qvec,isel,wt,id,n,nf)
C--   subroutine sqcGetDinv(def,dinv,nf,ierr)
C--
C--   double precision function dqcXSplne(idg,y,it)
C--   double precision function dqcBvalyt(idg,yy,tt)
C--   double precision function dqcBvalij(idg,iy,it)
C--   double precision function dqcFvalyt(idg,id,yy,tt)
C--   double precision function dqcFvalij(idg,id,iy,it)
C--   subroutine sqcAllFyt(idg,yy,tt,pdf,n)
C--   subroutine sqcAllFij(idg,iy,it,pdf,n)
C--   double precision function dqcFsumyt(idg,qvec,isel,yy,tt)
C--   double precision function dqcFsumij(idg,qvec,isel,iy,it)
C--   double precision function dqcSplChk(idg,it)

C=======================================================================
C==   Get q+- from e+-  ================================================
C=======================================================================

C--   Rollup algorithm, efficient to get all flavour pdfs in one call
C--
C--   Call Sn = Sum_n q+- and then rollup basis e+- starting with e6
C--   For example for 6 flavours, 4 active we have,
C--
C--   Basis pdfs              ->  Flavour pdfs
C--   ----------                  ------------
C--   e1 = S4                 ->  q1 = d = [S2 - e2]/2
C--   e2 = u - d              ->  q2 = u = [S2 + e2]/2
C--   e3 = S2 - 2s = S3 - 3s  ->  q3 = s = [S3 - e3]/3  and  S2 = S3 - s
C--   e4 = S3 - 3c = S4 - 4c  ->  q4 = c = [S4 - e4]/4  and  S3 = S4 - c
C--
C--   e5 = S4 - 4b            ->  q5 = b = [S4 - e5]/4
C--   e6 = S4 - 5t            ->  q6 = t = [S4 - e6]/5

C     ===================================================
      double precision function dqcOneQpm(id,ia,nf,nfmax)
C     ===================================================

C--   Roll-up basis pdfs to get one flavour pdf
C--   Not really faster than doing matrix multiplication
C--
C--   id    (in) : [1-6] identifier of d,u,s,c,b,t
C--   ia(6) (in) : addresses of basis pdfs e+- at (iy,iz)
C--   nf    (in) : active number of flavours nf(iz)
C--   nfmax (in) : max number of flavours nf(nz)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ia(6)

      dqcOneQpm = 0.D0
      qj        = 0.D0

      if(id.gt.nfmax) then
        dqcOneQpm = 0.D0
      elseif(id.gt.nf+1) then
        si = stor7(ia(1))
        ei = stor7(ia(id))
        dqcOneQpm = (si-ei)/(id-1)
      else
        j1 = nf
        j2 = max(id,2)
        sj = stor7(ia(1))
        do j = j1,j2,-1
          ej = stor7(ia(j))
          qj = (sj-ej)/j
          sj = sj-qj
        enddo
        if(id.eq.1) then
          dqcOneQpm = sj
        else
          dqcOneQpm = qj
        endif
      endif

      return
      end


C     ======================================
      subroutine sqcAllQpm(ia,nf,nfmax,qout)
C     ======================================

C--   Roll-up basis pdfs to get all flavour pdfs
C--   Faster than doing matrix multiplication
C--
C--   ia(6)    (in) : addresses of basis pdfs e+- at (iy,iz)
C--   nf       (in) : active number of flavours nf(iz)
C--   nfmax    (in) : max number of flavours nf(nz)
C--   qout(6) (out) : flavour pdfs q+-

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ia(6), qout(6)

      do i = 6,nfmax+1,-1
        qout(i) = 0.D0
      enddo
      do i = nfmax,nf+1,-1
        qout(i) = stor7(ia(i))
      enddo
      si = stor7(ia(1))
      do i = nf,2,-1
        ei      = stor7(ia(i))
        qout(i) = (si-ei)/i
        si      = si-qout(i)
      enddo
      qout(1) = si

      return
      end

C     =========================================
      subroutine sqcEfromQQ(qvec,evec,nf,nfmax)
C     =========================================

C--   Transform quark coefficients from flavor basis to si/ns basis
C--   Gluon is ignored
C--
C--          1  2  3  4  5  6  7  8  9 10 11 12 13
C--   qvec  tb bb cb sb ub db  g  d  u  s  c  b  t
C--   pdfid -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--
C--   |epm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  si 2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-
C--
C--   qvec(-6:6)   (in)  Coefficients in flavour space
C--                      qvec(0) is ignored --> quarks only!
C--   evec(12)     (out) Coefficients in si/ns space
C--   nf           (in)  Number of active flavours [3-6]
C--   nmfax        (in)  Max number of flavours

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpdfs7.inc'

      dimension qvec(-6:6),evec(12)

      mf = max(nf,nfmax)

C--   Transform
      do i = 1,12
        evec(i) = 0.D0
      enddo
      do i = 1,mf
        diplu = 0.D0
        dimin = 0.D0
        do j = 1,mf
          diplu = diplu + qvec( j)*tmatqen7(7+j,i+1,nf)
          diplu = diplu + qvec(-j)*tmatqen7(7-j,i+1,nf)
          dimin = dimin + qvec( j)*tmatqen7(7+j,i+7,nf)
          dimin = dimin + qvec(-j)*tmatqen7(7-j,i+7,nf)
        enddo
        evec(i  ) = diplu
        evec(i+6) = dimin
      enddo

      return
      end

C     ============================================
      subroutine sqcElistQQ(qvec,wt,id,n,nf,nfmax)
C     ============================================

C--   Transform quark coefficients from flavor basis to si/ns basis
C--   and produce a list of basis ids with their associated weights
C--   This routine ignores the gluon in qvec(0) and also does not
C--   include extra pdfs
C--
C--         -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   qvec  tb bb cb sb ub db  g  d  u  s  c  b  t
C--
C--   |qpm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--   |epm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  si 2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-
C--
C--   qvec(-6:6)  (in) :  Coefficients in flavour space
C--   wt(12)     (out) :  List of coefficients in si/ns space
C--   id(12)     (out) :  List of basis function identifiers [0,12]
C--   n          (out) :  Number of items in wt and id
C--   nf          (in) :  Number of active flavours [3-6]
C--   nfmax       (in) :  Max number of flavours
C--
C--   Remark: qvec(0) = gluon is ignored

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qpdfs7.inc'

      logical lmb_ne

      dimension qvec(-6:6),wt(12),id(12)

      mf = max(nf,nfmax)

C--   Transform mf flavours
      n = 0
      do i = 1,mf
        diplu = 0.D0
        dimin = 0.D0
        do j = 1,mf
          diplu = diplu + qvec( j)*tmatqen7(7+j,i+1,nf)
          diplu = diplu + qvec(-j)*tmatqen7(7-j,i+1,nf)
          dimin = dimin + qvec( j)*tmatqen7(7+j,i+7,nf)
          dimin = dimin + qvec(-j)*tmatqen7(7-j,i+7,nf)
        enddo
C--     Store non-zero coefficients and their identifiers
        if(lmb_ne(diplu,0.D0,aepsi6)) then
          n     = n+1
          wt(n) = diplu
          id(n) = i
        endif
        if(lmb_ne(dimin,0.D0,aepsi6)) then
          n     = n+1
          wt(n) = dimin
          id(n) = i+6
        endif
      enddo

      return
      end

C     ===========================================
      subroutine sqcElistFF(qvec,isel,wt,id,n,nf)
C     ===========================================

C--   Transform quark coefficients from flavor basis to si/ns basis
C--   and produce a list of basis ids with their associated weights
C--
C--         -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   qvec  tb bb cb sb ub db  g  d  u  s  c  b  t
C--
C--   |qpm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--   |epm>  0  1  2  3  4  5  6  7  8  9 10 11 12   13 ...
C--          g  si 2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-  f1 ...
C--
C--   qvec(-6:6)   (in) :  Coefficients in flavour space
C--   isel         (in) :  Selection flag
C--   wt(12)      (out) :  List of coefficients in si/ns space
C--   id(12)      (out) :  List of basis function identifiers [0,12+i]
C--   n           (out) :  Number of items in wt and id (may be zero)
C--   nf           (in) :  Number of active flavours [3-6]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qpdfs7.inc'

      logical lmb_ne

      dimension qvec(-6:6),evec(12),wt(12),id(12)

      dimension mask(0:12,9)
C--              g s           v
C--              0 1 2 3 4 5 6 7 8 9 0 1 2
      data mask /0,1,1,1,1,1,1,1,1,1,1,1,1,     !1=all quarks
     +           0,1,0,0,0,0,0,0,0,0,0,0,0,     !2=singlet
     +           0,0,1,1,1,1,1,1,1,1,1,1,1,     !3=all ns
     +           0,0,1,1,1,1,1,0,0,0,0,0,0,     !4=ns+
     +           0,0,0,0,0,0,0,1,1,1,1,1,1,     !5=v and ns-
     +           0,0,0,0,0,0,0,0,1,1,1,1,1,     !6=ns-
     +           0,0,0,0,0,0,0,1,0,0,0,0,0,     !7=valence
     +           1,0,0,0,0,0,0,0,0,0,0,0,0,     !8=gluon
     +           0,1,1,1,1,1,1,1,1,1,1,1,1  /   !9=all quarks

      if(isel.eq.0) then                        !gluon
        n     = 1
        wt(1) = 1.D0
        id(1) = 0

      elseif(isel.ge.1 .and. isel.le.9) then    !quarks or weigted gluon
        if(isel.ne.9) then
          mf = nf
        else
          mf = 6
        endif
        call sqcEfromQQ(qvec,evec,nf,mf)        !get basis coefficients
        if(isel.eq.8) then                      !weighted gluon
          n     = 1
          wt(1) = evec(1)
          id(1) = 0
        elseif(isel.eq.2)               then    !weighted singlet
          n     = 1
          wt(1) = evec(1)
          id(1) = 1
        elseif(isel.eq.7)               then    !weighted valence
          n     = 1
          wt(1) = evec(7)
          id(1) = 7
        else                                    !weighted quarks
          n = 0
          do i = 1,12
            if(lmb_ne(evec(i)*mask(i,isel),0.D0,aepsi6)) then
              n     = n+1
              wt(n) = evec(i)
              id(n) = i
            endif
          enddo
        endif

      elseif(isel.ge.13) then                   !extra pdf
        n     = 1
        wt(1) = 1.D0
        id(1) = isel

      else                                      !error
        stop 'sqcElistFF: wrong value of ISEL'
      endif

      return
      end

C     ============================================
      subroutine sqcGetDinv( def, dinv, nf, ierr )
C     ============================================

C--   The basis pdfs |e> are given by a linear combination of the 2nf
C--   user input pdfs |f> :
C--
C--               |ei> = Sum_j dinv(i,j) |fj>  j = 1,...2nf
C--
C--   This routine calculates the coefficients dinv(i,j)
C--
C--   def(-6:6,12)   (in) : def(i,j)   contribution of flavour i to |fj>
C--   dinv(12,12)   (out) : dinv(i,j)  contribution of |fj> to |ei>
C--   nf             (in) : number of active flavours [3-6]

      implicit double precision (a-h,o-z)
      dimension def(-6:6,12), dinv(12,12)
      dimension d(12), work(12,12), iw(12)

      do i = 1,2*nf
C--     Get d coefficients
        call sqcEfromQQ( def(-6,i), d, nf, nf )
C--     Copy d to work array
        do j = 1,nf
            work(i,j   ) = d(j  )
            work(i,j+nf) = d(j+6)
        enddo
      enddo
C--   Invert work array
      call smb_dminv( 2*nf, work, 12, iw, ierr )
      if( ierr.ne.0 ) return   !singular matrix
C--   Initialize dinv
      do i = 1,12
        do j = 1,12
          dinv(i,j) = 0.D0
        enddo
      enddo
C--   Copy work array to dinv
      do i = 1,nf
        do j = 1,2*nf
          dinv(i  ,j) = work(i   ,j)
          dinv(i+6,j) = work(i+nf,j)
        enddo
      enddo

      return
      end


C=======================================================================
C==   Pdf output routines ==============================================
C=======================================================================

C     =============================================
      double precision function dqcXSplne(idg,y,it)
C     =============================================

C--   Spline interpolation in x

C--   idg    (in)  : Global pdf id in stor7
C--   y      (in)  : value of y
C--   it     (in)  : t grid point

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      
      logical lmb_eq
      
      dimension acoef(mxx0)

C--   Catch y = 0  (x = 1)
      if(lmb_eq(y,0.D0,aepsi6)) then
        dqcXSplne = 0.D0
        return
      endif  
C--   Spline storage index
      idk = ioy2-1      
C--   Find y grid index in main grid G0
      iy = iqcFindIy(y)
C--   z grid index
      iz = izfit5(it)
C--   Convert to A coefficients
      call sqcGetSplA(stor7,idg,iy,iz,ig,iyg,acoef)
C--   We did hit the end-point: iyg --> iyg-1
      iyg = min(iyg,nyy2(ig)-1)
C--   Spline index range in y
      call sqcByjLim(idk,iyg+1,minby,maxby)  !iy --> iy+1
C--   Loop over spline y
      val = 0.D0
      do jy = minby,maxby
        yjm1  = (jy-1)*dely2(ig)
        spy = dqcBsplyy(idk,1,(y-yjm1)/dely2(ig))
        val = val + acoef(jy)*spy
      enddo

      dqcXSplne = val

      return
      end
      
C     ==============================================
      double precision function dqcBvalyt(idg,yy,tt)
C     ==============================================

C--   Value of an internal basis pdf at y and t
C--
C--   id  =  0   1   2   3   4   5   6   7   8   9  10  11  12  13 ...
C--          g  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2- e3- e4- e5- e6- extra
C--
C--   idg   (in) : stor7 pdf identifier in global format
C--   yy,tt (in) : interpolation point
      
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qpars6.inc'
      
      logical lmb_eq
      dimension wy(6),wz(6)

C--   Catch y = 0  (x = 1)      
      if(lmb_eq(yy,0.D0,aepsi6)) then
        dqcBvalyt = 0.D0
      else
        call sqcZmesh(yy,tt,0,iy1,iy2,iz1,iz2,it1)
        ny  = iy2-iy1+1
        nz  = iz2-iz1+1
        call sqcIntWgt(iy1,ny,it1,nz,yy,tt,wy,wz)
        iag = iqcG5ijk(stor7,iy1,iz1,idg)
        dqcBvalyt = dqcPdfPol(stor7,iag,ny,nz,wy,wz)
      endif

      return
      end
      
C     ==============================================
      double precision function dqcBvalij(idg,iy,it)
C     ==============================================

C--   Value of an internal pdf at iy and it
C--
C--   id  =  0   1   2   3   4   5   6   7   8   9  10  11  12  13 ...
C--          g  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2- e3- e4- e5- e6- extra
C--
C--   idg  ( in) : stor7 pdf identifier in global format
C--   ksetw (in) : stor7 set identifier
C--   iy,it (in) : grid point
C--
C--   Remark: it > 0 then follow the 4,5,6 convention at iqc,b,t
C--              < 0 then follow the 3,4,5 convention at iqc,b,t
      
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      iz        = izfit5(it)
      ia        = iqcG5ijk(stor7,iy,iz,idg)
      dqcBvalij = stor7(ia)

      return
      end

C     =================================================
      double precision function dqcFvalyt(idg,id,yy,tt)
C     =================================================

C--   Return gluon or q, qbar
C--
C--   |id|   =  0  1  2  3  4  5  6  7 ...
C--   q/qb   =  g  d  u  s  c  b  t f1 ...
C--
C--   idg    (in) : gluon pdf id in stor7 (global format)
C--   ksetw  (in) : pdf set id in stor7
C--   id     (in) : local pdf identifier [-6,6+n]
C--   y,t    (in) : interpolation point

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qpars6.inc'

      logical lmb_eq
      dimension wy(6),wz(6)
      dimension qvec(-6:6),ide(12),wte(12)

C--   Initialise
      call smb_VFill(qvec,13,0.D0)

C--   Catch y = 0  (x = 1)
      if(lmb_eq(yy,0.D0,aepsi6)) then
        dqcFvalyt = 0.D0
      else
        call sqcZmesh(yy,tt,0,iy1,iy2,iz1,iz2,it1)
        ny  = iy2-iy1+1
        nz  = iz2-iz1+1
        call sqcIntWgt(iy1,ny,it1,nz,yy,tt,wy,wz)

        iag   = iqcG5ijk(stor7,iy1,iz1,idg)
        idabs = abs(id)
        it    = iqcItfrmt(tt)
        if(it.eq.0)  stop 'sqcQQByt: t out of range ---> STOP'
        iz    = izfit5(it)
        nf    = itfiz5(-iz)
        if(nfix5.eq.1) then
C--       intrinsic heavy flavours
          nfmax = itfiz5(-izmac5)
        else
C--       no intrinsic heavy flavours
          nfmax =  nf
        endif

        if(id.eq.0)     then
C--       Output gluon
          dqcFvalyt = dqcPdfPol(stor7,iag,ny,nz,wy,wz)
        elseif(id.gt.6) then
C--       Output extra pdf
          dqcFvalyt = dqcPdfPol(stor7,iag+(id+6)*incid7,ny,nz,wy,wz)
        elseif(idabs.gt.nfmax) then
C--       Output nonexisting flavour
          dqcFvalyt = 0.D0
        elseif(idabs.gt.nf) then
C--       Output non-active flavour
          eplus = dqcPdfPol(stor7,iag+ idabs   *incid7,ny,nz,wy,wz)
          eminu = dqcPdfPol(stor7,iag+(idabs+6)*incid7,ny,nz,wy,wz)
          if(id.gt.0) then
            dqcFvalyt = 0.5D0*(eplus+eminu)
          else
            dqcFvalyt = 0.5D0*(eplus-eminu)
          endif
        else
C--       Output quark
C--       Get list of weights and basis ids [1,12]
          qvec(id) = 1.D0
          call sqcElistQQ(qvec,wte,ide,n,nf,nf)
          qvec(id) = 0.D0
          sum      = 0.D0
C--       Go for weighted sum
          do i = 1,n
            ia  = iag + ide(i) * incid7
            sum = sum + wte(i) * dqcPdfPol(stor7,ia,ny,nz,wy,wz)
          enddo
          dqcFvalyt = sum
        endif
      endif

      return
      end

C     =================================================
      double precision function dqcFvalij(idg,id,iy,it)
C     =================================================

C--   Return gluon or 0.5(qplus +- qminus)
C--
C--   |id|   =  0  1  2  3  4  5  6  7  ...
C--   q/qb   =  g  d  u  s  c  b  t f1  ...
C--
C--   idg    (in) : gluon pdf id in stor7 (global format)
C--   id     (in) : local pdf identifier [-6,6+n]
C--   iy,it  (in) : grid point

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qpars6.inc'

      dimension qvec(-6:6),ide(12),wte(12)

      call smb_VFill(qvec,13,0.D0)

      iz    = izfit5(it)
      iag   = iqcG5ijk(stor7,iy,iz,idg)
      idabs = abs(id)
      nf    = itfiz5(-iz)
      if(nfix5.eq.1) then
C--     intrinsic heavy flavours
        nfmax = itfiz5(-izmac5)
      else
C--    no intrinsic heavy flavours
        nfmax =  nf
      endif

      if(id.eq.0)     then
C--     Output gluon
        dqcFvalij = stor7(iag)
      elseif(id.gt.6) then
C--     Output extra pdf
        dqcFvalij = stor7(iag+(id+6)*incid7)
      elseif(idabs.gt.nfmax) then
C--     Output non-pexisting flavour
        dqcFvalij = 0.D0
      elseif(idabs.gt.nf) then
C--     Output non-active flavour
        eplus = stor7(iag+ idabs   *incid7)
        eminu = stor7(iag+(idabs+6)*incid7)
        if(id.gt.0) then
          dqcFvalij = 0.5D0*(eplus+eminu)
        else
          dqcFvalij = 0.5D0*(eplus-eminu)
        endif
      else
C--     Output quark
C--     Get list of weights and basis ids [1,12]
        qvec(id) = 1.D0
        call sqcElistQQ(qvec,wte,ide,n,nf,nf)
        qvec(id) = 0.D0
        sum      = 0.D0
C--     Go for weighted sum
        do i = 1,n
          ia  = iag + ide(i) * incid7
          sum = sum + wte(i) * stor7(ia)
        enddo
        dqcFvalij = sum
      endif

      return
      end

C     =====================================
      subroutine sqcAllFyt(idg,yy,tt,pdf,n)
C     =====================================

C--   Return all flavour pdfs
C--
C--   |id|   =  0  1  2  3  4  5  6  7  ...
C--   q/qb   =  g  d  u  s  c  b  t f1  ...
C--
C--   idg         (in)  : gluon pdf id in stor7 (global format)
C--   yy,tt       (in)  : interpolation point
C--   pdf(-6:6+n) (out) : output pdfs
C--   n           (in)  : number of extra pdfs

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qpars6.inc'

      logical lmb_eq
      dimension wy(6),wz(6),pdf(-6:6+n)

C--   Initialise
      do i = -6,6+n
        pdf(i) = 0.D0
      enddo
C--   Catch y = 0  (x = 1)
      if(lmb_eq(yy,0.D0,aepsi6)) return
C--   Figure out number of flavors
      it = iqcItfrmt(tt)
      if(it.eq.0)  stop 'sqcAllFyt: t out of range ---> STOP'
      iz = izfit5( it)
      nf = itfiz5(-iz)
      if(nfix5.eq.1) then
C--     intrinsic heavy flavours
        nfmax = itfiz5(-izmac5)
      else
C--     no intrinsic heavy flavours
        nfmax =  nf
      endif
C--   Setup interpolation
      call sqcZmesh(yy,tt,0,iy1,iy2,iz1,iz2,it1)
      ny  = iy2-iy1+1
      nz  = iz2-iz1+1
      call sqcIntWgt(iy1,ny,it1,nz,yy,tt,wy,wz)

      iag = iqcG5ijk(stor7,iy1,iz1,idg)

C--   Output gluon
      pdf(0) = dqcPdfPol(stor7,iag,ny,nz,wy,wz)
C--   Output extra pdfs
      do id = 7,6+n
        pdf(id) = dqcPdfPol(stor7,iag+(id+6)*incid7,ny,nz,wy,wz)
      enddo
C--   Output non-existing flavours
      do id = nfmax+1,6
        pdf( id) = 0.D0
        pdf(-id) = 0.D0
      enddo
C--   Output non-active flavours
      do id = nf+1,nfmax
        eplus = dqcPdfPol(stor7,iag+ id   *incid7,ny,nz,wy,wz)
        eminu = dqcPdfPol(stor7,iag+(id+6)*incid7,ny,nz,wy,wz)
        pdf( id) = 0.5D0*(eplus+eminu)
        pdf(-id) = 0.5D0*(eplus-eminu)
      enddo
C--   Output quarks
      do id = 1,nf
        iap    = iag
        iam    = iag+6*incid7
        sumpl  = 0.D0
        sumin  = 0.D0
        do j = 1,nf
          iap   = iap + incid7
          eplus = dqcPdfPol(stor7,iap,ny,nz,wy,wz)
          sumpl = sumpl + umatqe7(id,j,nf)*eplus
          iam   = iam + incid7
          eminu = dqcPdfPol(stor7,iam,ny,nz,wy,wz)
          sumin = sumin + umatqe7(id,j,nf)*eminu
        enddo
        pdf( id) = 0.5D0 * (sumpl+sumin)
        pdf(-id) = 0.5D0 * (sumpl-sumin)
      enddo

      return
      end

C     =====================================
      subroutine sqcAllFij(idg,iy,it,pdf,n)
C     =====================================

C--   Return all flavour pdfs
C--
C--   |id|   =  0  1  2  3  4  5  6
C--   q/qb   =  g  d  u  s  c  b  t
C--
C--   idg       (in)  : gluon pdf id in stor7 (global format)
C--   iy,it     (in)  : grid point
C--   pdf(-6,6) (out) : output pdfs
C--   n         (in)  : number of extra pdfs

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension pdf(-6:6+n)

      iz    = izfit5(it)
      nf    = itfiz5(-iz)
      if(nfix5.eq.1) then
C--     intrinsic heavy flavours
        nfmax = itfiz5(-izmac5)
      else
C--     no intrinsic heavy flavours
        nfmax =  nf
      endif
      iag   = iqcG5ijk(stor7,iy,iz,idg)

C--   Output gluon
      pdf(0) = stor7(iag)
C--   Output extra pdfs
      do id = 7,6+n
        pdf(id) = stor7(iag+(id+6)*incid7)
      enddo
C--   Output non-existing flavours
      do id = nfmax+1,6
        pdf( id) = 0.D0
        pdf(-id) = 0.D0
      enddo
C--   Output non-active flavours
      do id = nf+1,nfmax
        eplus = stor7(iag+ id   *incid7)
        eminu = stor7(iag+(id+6)*incid7)
        pdf( id) = 0.5D0*(eplus+eminu)
        pdf(-id) = 0.5D0*(eplus-eminu)
      enddo
C--   Output quarks
      do id = 1,nf
        iap    = iag
        iam    = iag+6*incid7
        sumpl  = 0.D0
        sumin  = 0.D0
        do j = 1,nf
          iap   = iap + incid7
          eplus = stor7(iap)
          sumpl = sumpl + umatqe7(id,j,nf)*eplus
          iam   = iam + incid7
          eminu = stor7(iam)
          sumin = sumin + umatqe7(id,j,nf)*eminu
        enddo
        pdf( id) = 0.5D0 * (sumpl+sumin)
        pdf(-id) = 0.5D0 * (sumpl-sumin)
      enddo

      return
      end

C     ========================================================
      double precision function dqcFsumyt(idg,qvec,isel,yy,tt)
C     ========================================================

C--   Return linear combination of q and qbar
C--   id     =  -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   def    =  tb bb cb sb ub db  g  d  u  s  c  b  t
C--
C--   idg         (in) : gluon id in stor7 (global format)
C--   qvec(-6:6)  (in) : coefficients of the linear combination
C--   isel        (in) : selection flag
C--   y,t         (in) : interpolation point
C--
C--   NB: qvec(0) corresponds to the gluon and is ignored

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qpars6.inc'

      logical lmb_eq

      dimension wy(6),wz(6),qvec(-6:6),wte(12),ide(12)

C--   Catch y = 0  (x = 1)
      if(lmb_eq(yy,0.D0,aepsi6)) then
        sum = 0.D0
      else
C--     Figure out number of flavors
        it = iqcItfrmt(tt)
        if(it.eq.0)  stop 'sqcSumQQByt: t out of range ---> STOP'
        iz = izfit5( it)
        nf = itfiz5(-iz)
        if(nfix5.eq.1) then
C--       intrinsic heavy flavours
          nfmax = itfiz5(-izmac5)
        else
C--       no intrinsic heavy flavours
          nfmax =  nf
        endif
C--     Setup interpolation
        call sqcZmesh(yy,tt,0,iy1,iy2,iz1,iz2,it1)
        ny  = iy2-iy1+1
        nz  = iz2-iz1+1
        call sqcIntWgt(iy1,ny,it1,nz,yy,tt,wy,wz)
C--     Setup weighted sum of quark basis functions
        call sqcElistFF(qvec,isel,wte,ide,n,nf)         !n might be zero
C--     Base address of gluon
        iag = iqcG5ijk(stor7,iy1,iz1,idg)
C--     Go for weighted sum
        sum = 0.D0
        do i = 1,n                                      !n might be zero
          ia  = iag + ide(i) * incid7
          sum = sum + wte(i) * dqcPdfPol(stor7,ia,ny,nz,wy,wz)
        enddo
C--     Add intrinsic heavy quarks for isel = 9
        if(isel.eq.9) then
          do id = nf+1,nfmax
            eplus = dqcPdfPol(stor7,iag+id*incid7,ny,nz,wy,wz)
            sum   = sum+eplus
          enddo
        endif
      endif

      dqcFsumyt = sum

      return
      end

C     ========================================================
      double precision function dqcFsumij(idg,qvec,isel,iy,it)
C     ========================================================

C--   Return linear combination of q and qbar
C--   id     =  -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   def    =  tb bb cb sb ub db  g  d  u  s  c  b  t
C--
C--   idg         (in) : gluon id in stor7 (global format)
C--   qvec(-6:6)  (in) : coefficients of the linear combination
C--   isel        (in) : selection flag
C--   iy,it       (in) : grid point
C--
C--   NB: qvec(0) corresponds to the gluon and is ignored

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension qvec(-6:6),wte(12),ide(12)

      iz = izfit5( it)
      nf = itfiz5(-iz)
      if(nfix5.eq.1) then
C--     intrinsic heavy flavours
        nfmax = itfiz5(-izmac5)
      else
C--     no intrinsic heavy flavours
        nfmax =  nf
      endif
C--   Setup weighted sum of quark basis functions
      call sqcElistFF(qvec,isel,wte,ide,n,nf)
C--   Base address of gluon
      iag = iqcG5ijk(stor7,iy,iz,idg)
C--   Go for weighted sum
      sum = 0.D0
      do i = 1,n
        ia  = iag + ide(i) * incid7
        sum = sum + wte(i) * stor7(ia)
      enddo
C--   Add intrinsic heavy quarks for isel = 9
      if(isel.eq.9) then
        do id = nf+1,nfmax
          eplus = stor7(iag+id*incid7)
          sum   = sum+eplus
        enddo
      endif

      dqcFsumij = sum

      return
      end
      
C     ===========================================
      double precision function dqcSplChk(idg,it)
C     ===========================================

C--   Returns epsi = || quad-lin || interpolation at midpoints
C--
C--   idg   (in) : stor7 pdf identifier in global format
C--   it    (in) : mu2 grid point

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      
      dimension acoef(mxx0), epsi(mxx0)

C--   Initialize
      dqcSplChk = 0.D0
C--   Linear interpolation
      if(ioy2.ne.3) return
C--   Find z-bin
      iz = izfit5(it)

C--   Loop over subgrids
      do jg = 1,nyg2
        iy = iyma2(jg)
        call sqcGetSplA(stor7,idg,iy,iz,ig,ny,acoef)
C--     Debug checks
        if(ig.ne.jg)        stop 'dqcSplChk: ig not jg'
        if(ny.ne.nyy2(jg))  stop 'dqcSplChk: ny not nyy2(jg)'
C--     Now get vector of deviations 
        nyma = iqcIyMaxG(iymac2,ig)       
        call sqcDHalf(ioy2,acoef,epsi,nyma)
C--     Max deviation
        do iy = 1,nyma
          dqcSplChk = max(dqcSplChk,abs(epsi(iy)))
        enddo  
      enddo
      
      return
      end
