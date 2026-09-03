
C--   This is the file srcsplint.f with free-runnung workhorses

C--   sspYnMake(istep, ynodes, nys, ierr)
C--   sspTnMake(istep, tnodes, nts, ierr)
C--   sspYnUser(xarr, nx, ynodes, nys, ierr)
C--   sspTnUser(qarr, nq, tnodes, nts, ierr)
C--   ispS1make(w, unodes, nu, isign)
C--   sspS1fill(w, iasp, fvals)
C--   dspS1fun(w, ia, u)
C--   sspGetIaOneD(w, ia, iat, iau, nus, iaf, iab, iac, iad)

C--   dspBintYi(w, iasp, iy, y)
C--   dspSpIntY(w, iasp, ymin, ymax)
C--   dspBintTi(w, iasp, it, t)
C--   dspSpIntT(w, iasp, tmin, tmax)

C--   ispS2make(w, unodes, nu, vnodes, nv)
C--   sspS2fill(w, iasp2, fvals)
C--   sspRangeYT(w, ia, rscut)
C--   dspS2fun(w, ia, u, v)
C--   sspGetIaTwoD(w, ia, iat, iau, nus, iav, nvs, iaFF, iaCC)

C--   ispN2make(w, unodes, nu, vnodes, nv)

C--   dspBintYij(w, iasp, y1, y2, dt)
C--   dspBintTij(w, iasp, t1, t2, dy)
C--   dspBintYTij(w, iasp, y1, y2, t1, t2)

C--   dspBintYYTT(w, iasp, iy, it, y1, y2, t1, t2, rs, np)
C--   dspGausFun(t)
C--   dspSpIntYT(w, iasp, y1, y2, t1, t2, rs, np)
C--   sspSnipSnip(y1, y2, t1, t2, sc, ta, tb)
C--   sspGetCoefs(w, iasp, iy, it, y1, t1)
C--   dspDerSp2(w, iasp, n, m, iy, it, y, t)

C--   ispSplineType(w, ia)
C--   ispReadOnly(w, ia)
C--   sspSpLims(w, ia, nu, umi, uma, nv, vmi, vma, ndim, nb)

C--   ispIaFromI(w, iasp, i)
C--   ispIyFromY(w, ias2, y)
C--   ispItFromT(w, ias2, t)
C--   lspIsaFbin(w, ias2, iy, it)
C--   sspBinLims(w, ias2, iy, it, y1, y2, t1, t2)
C--   ispCrossSc(y1, y2, t1, t2, sc)
C--   dspRsMax(w, ias2, sc)

C--   dspPol3(x, coef, n)
C--   ispGetBin(x, xarray, n)
C--   sspEplus(x, Eplus)
C--   sspEminu(x, Eminu)


C=======================================================================
C===  1-dim splines  ===================================================
C=======================================================================

C     ==============================================
      subroutine sspYnMake(istep, ynodes, nys, ierr)
C     ==============================================

C--   Automake xnodes and convert to ynodes

C--   istep   (in): x-grid step size
C--   ynodes (out): array of ynodes in ascending order
C--   nys    (out): number of ynodes
C--   ierr   (out): 0 = OK, 1 = array size exceeded

      implicit double precision(a-h,o-z)

      include 'splint.inc'

      dimension ynodes(*)

C--   Sample evolution grid
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)

C--   Fill y-nodes
      ynodes(1) = 0.D0                 !boundary x = 1
      nys       = 1
      do ix = nx,2,-istep
        if(nys.gt.maxn0-2) goto 500
        nys         =  nys+1
        ynodes(nys) = -log(xfrmix(ix))
      enddo
      nys         =  nys+1
      ynodes(nys) = -log(xfrmix(1))    !boundary x = xmin
      ierr        =  0
      return

 500  continue
      ierr = 1

      return
      end

C     ==============================================
      subroutine sspTnMake(istep, tnodes, nts, ierr)
C     ==============================================

C--   Automake qnodes and convert to tnodes

C--   istep   (in): q-grid step size
C--   tnodes (out): array of tnodes in ascending order
C--   nts    (out): number of tnodes
C--   ierr   (out): 0 = OK, 1 = array size exceeded

      implicit double precision(a-h,o-z)

      include 'splint.inc'

      dimension tnodes(*)

C--   Sample evolution grid
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)

C--   Fill t-nodes
      tnodes(1) = log(qfrmiq(1))        !lower boundary qmin
      nts = 1
      do iq = 2,nq-1,istep
        if(nts.gt.maxn0-2) goto 500
        nts         = nts+1
        tnodes(nts) = log(qfrmiq(iq))
      enddo
      nts         = nts+1
      tnodes(nts) = log(qfrmiq(nq))     !upper boundary qmax
      ierr        = 0
      return

 500  continue
      ierr = 1

      return
      end

C     =================================================
      subroutine sspYnUser(xarr, nx, ynodes, nys, ierr)
C     =================================================

C--   Weed and sort xnodes and convert to ynodes

C--   xarr    (in): x-node user input
C--   nx      (in): number of user x-nodes
C--   ynodes (out): array of ynodes in ascending order
C--   nys    (out): number of ynodes
C--   ierr   (out): 0 = OK, 1 = array size exceeded

      implicit double precision(a-h,o-z)

      include 'splint.inc'

      dimension xarr(*), ynodes(*)
      real      rix(maxn0)

C--   Discard points outside grid
      ny = 0
      do i = 1,nx
        x    = xarr(i)
        ix   = ixfrmx(x)
        if(ix.ne.0) then
          if(ny.gt.maxn0-1) goto 500
          ny       = ny+1
          rix(ny)  = real(ix)
        endif
      enddo
C--   Sort and discard equal entries
      call smb_asort ( rix, ny, nys )
C--   Fill y-nodes
      do iy = 1,nys
        ix         =  int(rix(nys+1-iy))
        x          =  xfrmix(ix)
        ynodes(iy) = -log(x)
      enddo
      ierr = 0
      return

 500  continue
      ierr = 1
    
      return
      end

C     =================================================
      subroutine sspTnUser(qarr, nq, tnodes, nts, ierr)
C     =================================================

C--   Weed and sort qnodes and convert to tnodes

C--   qarr    (in): q-node user input
C--   nq      (in): number of user q-nodes
C--   tnodes (out): array of tnodes in ascending order
C--   nts    (out): number of tnodes
C--   ierr   (out): 0 = OK, 1 = array size exceeded

      implicit double precision(a-h,o-z)

      include 'splint.inc'

      dimension qarr(*), tnodes(*)
      real      riq(maxn0)

C--   Discard points outside grid
      nt = 0
      do i = 1,nq
        q    = qarr(i)
        iq   = iqfrmq(q)
        if(iq.ne.0) then
          if(nt.gt.maxn0-1) goto 500
          nt       = nt+1
          riq(nt)  = real(iq)
        endif
      enddo
C--   Sort and discard equal entries
      call smb_asort ( riq, nt, nts )
C--   Fill t-nodes
      do it = 1,nts
        iq         =  int(riq(it))
        q          =  qfrmiq(iq)
        tnodes(it) =  log(q)
      enddo
      ierr = 0
      return

 500  continue
      ierr = 1
    
      return
      end

C     ================================================
      integer function ispS1make(w, unodes, nu, isign)
C     ================================================

C--   Create 1-dimensional spline object

C--   w       (in): workspace
C--   unodes  (in): array of input u-nodes
C--   nu      (in): number of unodes
C--   isign   (in): log transformation sign -1 for x and +1 for q

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*), unodes(*), imi(1), ima(1)

C--   Create new set
      iasp  = iws_newset(w)
C--   Create user table
      imi     = 1
      ima     = nusr0
      iausr   = iws_wtable(w,imi,ima,1)
      iuser   = iws_BeginTbody(w,iausr)
C--   Create and fill u-node table
      imi   = 1
      ima   = nu
      itabu = iws_wtable(w,imi,ima,1)
      ia    = iws_BeginTbody(w,itabu)-1
      do iu = 1,nu
        w(ia+iu) = unodes(iu)
      enddo
C--   Create coefficient arrays
      imi   = 1
      ima   = nu
      itabf = iws_wtable(w,imi,ima,1)
      itabb = iws_wtable(w,imi,ima,1)
      itabc = iws_wtable(w,imi,ima,1)
      itabd = iws_wtable(w,imi,ima,1)
C--   Store local addresses in tag field
      ia = iws_IaFirstTag(w,iasp)
      w(ia+ImarkS0) = dble(MarkSp0)
      w(ia+IdimSp0) = dble(isign)
      w(ia+NedegU0) = dble(3)
      w(ia+NedegV0) = dble(3)
      w(ia+IsUtab0) = dble(itabu-iasp)
      w(ia+NwUtab0) = dble(nu)
      w(ia+IsVtab0) = dble(0)
      w(ia+NwVtab0) = dble(0)
      w(ia+Nodesa0) = dble(nu)
      w(ia+IsUser0) = dble(iuser-iasp)
      w(ia+IsFtab0) = dble(itabf-iasp)
      w(ia+IsBtab0) = dble(itabb-iasp)
      w(ia+IsCtab0) = dble(itabc-iasp)
      w(ia+IsDtab0) = dble(itabd-iasp)
C--   Store address of first spline
      iaR = iws_iaRoot()
      iat = iws_IaFirstTag(w,iaR)
      if(int(w(iat+IaSfirst0)).eq.0) w(iat+IaSfirst0) = dble(iasp)
C--   Return address
      ispS1make     = iasp

      return
      end

C     ====================================
      subroutine sspS1fill(w, iasp, fvals)
C     ====================================

C--   Fill and construct 1-dimensional spline

C--   w       (in): workspace
C--   fvals   (in): array of input f(u_i)

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*), fvals(*)

C--   Get table-body addresses
      call sspGetIaOneD(w,iasp,iat,iau,nus,iaf,iab,iac,iad)
C--   Fill f
      ia    = iaf-1
      do iu = 1,nus
        w(ia+iu) = fvals(iu)
      enddo
C--   Make spline coefficients
      call smb_spline(nus,w(iau),w(iaf),w(iab),w(iac),w(iad))

      return
      end

C     ============================================
      double precision function dspS1fun(w, ia, u)
C     ============================================

C--   Compute 1-dim spline function

      implicit double  precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

      save iarem, iau, iaf, iab, iac, iad, nus
      data iarem/0/

      dimension coef(0:3)

C--   Get addresses
      if(ia.ne.iarem) then
        call sspGetIaOneD(w,ia,iat,iau,nus,iaf,iab,iac,iad)
        iarem = ia
      endif
C--   Get node bin (< 0 = outside range)
      iu   = ispGetBin(u,w(iau),nus)
      ium1 = abs(iu)-1
      ndeg = 3
      if(iu.lt.0) ndeg = int(w(iat+NedegU0))
C--   Get coefficients
      coef(0) = w(iaf+ium1)
      coef(1) = w(iab+ium1)
      coef(2) = w(iac+ium1)
      coef(3) = w(iad+ium1)
C--   Spline function
      du = u - w(iau+ium1)
      dspS1fun = dspPol3(du,coef,ndeg)

      return
      end

C     =================================================================
      subroutine sspGetIaOneD(w, ia, iat, iau, nus, iaf, iab, iac, iad)
C     =================================================================

C--   Get table body addresses of 1-dim spline

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

      iat   = iws_IaFirstTag(w,ia)
      nus   = int(w(iat+NwUtab0))

      itabu = int(w(iat+IsUtab0)) + ia
      itabf = int(w(iat+IsFtab0)) + ia
      itabb = int(w(iat+IsBtab0)) + ia
      itabc = int(w(iat+IsCtab0)) + ia
      itabd = int(w(iat+IsDtab0)) + ia

      iau   = iws_BeginTbody(w,itabu)
      iaf   = iws_BeginTbody(w,itabf)
      iab   = iws_BeginTbody(w,itabb)
      iac   = iws_BeginTbody(w,itabc)
      iad   = iws_BeginTbody(w,itabd)

      return
      end

C=======================================================================
C===  1-dim integration  ===============================================
C=======================================================================

C     ===================================================
      double precision function dspBintYi(w, iasp, iy, y)
C     ===================================================

C--   Compute Int_yi^y exp(-u) Pol3(u) du

C--   w     (in): workspace
C--   iasp  (in): address of 1-dim x-spline
C--   iy    (in): node-point index in y must be below upper index limit
C--   y     (in): y-coordinate inside node-bin iy

      implicit double precision (a-h,o-z)
      logical lmb_le, lmb_gt

      include 'splint.inc'

      dimension w(*), eminu(0:3), coef(0:3)

C--   Get table addresses
      call sspGetIaOneD(w,iasp,iat,iau,nus,iaf,iab,iac,iad)
C--   Lower and upper bin limits
      iym1 = iy-1
      ymin = w(iau+iym1)
      ymax = w(iau+iy)
C--   Range check
      if(lmb_le(y,ymin,-deps0) .or. lmb_gt(y,ymax,-deps0)) then
        dspBintYi = 0.D0
        return
      endif
C--   Get E-integrals
      dy = y-ymin
      call sspEminu(dy,eminu)
C--   Get spline coefficients
      coef(0) = w(iaf+iym1)
      coef(1) = w(iab+iym1)
      coef(2) = w(iac+iym1)
      coef(3) = w(iad+iym1)
C--   Compute integral
      val = 0.D0
      do i = 0,3
        val = val + coef(i)*eminu(i)
      enddo
C--   Done (almost)
      dspBintYi = exp(-ymin)*val

      return
      end

C     ========================================================
      double precision function dspSpIntY(w, iasp, ymin, ymax)
C     ========================================================

C--   Compute Int_ymin^ymax exp(-u) P(u) du

C--   w     (in): workspace
C--   iasp  (in): address of 1-dim x-spline
C--   ymin  (in): lower limit inside range of spline < ymax
C--   ymax  (in): upper limit inside range of spline > ymin

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

C--   Get node table address
      ia    = iws_IaFirstTag(w,iasp)
      nus   = int(w(ia+NwUtab0))
      itabu = int(w(ia+IsUtab0)) + iasp
      iau   = iws_BeginTbody(w,itabu)
C--   Find node index limits
      iymi  = ispGetBin(ymin,w(iau),nus)
      iyma  = ispGetBin(ymax,w(iau),nus)
C--   Catch out of range
      if(iymi.le.0 .or. iyma.le.0) stop
     +   ' SPLINT::dspSpIntY: problem with limits out of range'
C--   Reset bin-number when ymax = upper spline limit
      if(iyma.eq.nus) iyma = iyma-1
C--   Compute integral
      val = 0.D0
      do iy = iymi,iyma-1
        yma = w(iau+iy)
        val = val + dspBintYi(w,iasp,iy,yma)
      enddo
      val = val + dspBintYi(w,iasp,iyma,ymax)
      val = val - dspBintYi(w,iasp,iymi,ymin)

      dspSpIntY = val

      return
      end

C     ===================================================
      double precision function dspBintTi(w, iasp, it, t)
C     ===================================================

C--   Compute Int_ti^t exp(u) Pol3(u) du

C--   w     (in): workspace
C--   iasp  (in): address of 1-dim q-spline
C--   it    (in): node-point index in t must be below upper index limit
C--   t     (in): t-coordinate inside node-bin it

      implicit double precision (a-h,o-z)
      logical lmb_le, lmb_gt

      include 'splint.inc'

      dimension w(*), eplus(0:3), coef(0:3)

C--   Get table addresses
      call sspGetIaOneD(w,iasp,iat,iau,nus,iaf,iab,iac,iad)
C--   Lower and upper bin limits
      itm1 = it-1
      tmin = w(iau+itm1)
      tmax = w(iau+it)
C--   Range check
      if(lmb_le(t,tmin,-deps0) .or. lmb_gt(t,tmax,-deps0)) then
        dspBintTi = 0.D0
        return
      endif
C--   Get E-integrals
      dt = t-tmin
      call sspEplus(dt,eplus)
C--   Get spline coefficients
      coef(0) = w(iaf+itm1)
      coef(1) = w(iab+itm1)
      coef(2) = w(iac+itm1)
      coef(3) = w(iad+itm1)
C--   Compute integral
      val = 0.D0
      do i = 0,3
        val = val + coef(i)*eplus(i)
      enddo
C--   Done (almost)
      dspBintTi = exp(tmin)*val

      return
      end

C     ========================================================
      double precision function dspSpIntT(w, iasp, tmin, tmax)
C     ========================================================

C--   Compute Int_tmin^tmax exp(u) P(u) du

C--   w     (in): workspace
C--   iasp  (in): address of 1-dim q-spline
C--   tmin  (in): lower limit inside range of spline < tmax
C--   tmax  (in): upper limit inside range of spline > tmin

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

C--   Get node table address
      ia    = iws_IaFirstTag(w,iasp)
      nus   = int(w(ia+NwUtab0))
      itabu = int(w(ia+IsUtab0)) + iasp
      iau   = iws_BeginTbody(w,itabu)
C--   Find node index limits
      itmi  = ispGetBin(tmin,w(iau),nus)
      itma  = ispGetBin(tmax,w(iau),nus)
C--   Catch out of range
      if(itmi.le.0 .or. itma.le.0) stop
     +   ' SPLINT::dspSpIntT: problem with limits out of range'
C--   Reset bin-number when tmax = upper spline limit
      if(itma.eq.nus) itma = itma-1
C--   Compute integral
      val = 0.D0
      do it = itmi,itma-1
        tma = w(iau+it)
        val = val + dspBintTi(w,iasp,it,tma)
      enddo
      val = val + dspBintTi(w,iasp,itma,tmax)
      val = val - dspBintTi(w,iasp,itmi,tmin)

      dspSpIntT = val

      return
      end

C=======================================================================
C===  2-dim splines  ===================================================
C=======================================================================

C     =====================================================
      integer function ispS2make(w, unodes, nu, vnodes, nv)
C     =====================================================

C--   Create 2-dimensional spline object

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*), unodes(*), vnodes(*)
      dimension imin(4), imax(4)

C--   Create new set
      iasp2   = iws_newset(w)
C--   Create user table
      imin(1) = 1
      imax(1) = nusr0
      iausr   = iws_wtable(w,imin,imax,1)
      iuser   = iws_BeginTbody(w,iausr)
C--   Create and fill u-node table (2nd dim = range table filled later)
      imin(1) = 1
      imax(1) = nu
      imin(2) = 1
      imax(2) = 2
      itabu   = iws_wtable(w,imin,imax,2)
      ia      = iws_BeginTbody(w,itabu)-1
      do i = 1,nu
        w(ia+i) = unodes(i)
      enddo
C--   Create and fill v-node table (2nd dim = range table filled later)
      imin(1) = 1
      imax(1) = nv
      imin(2) = 1
      imax(2) = 2
      itabv   = iws_wtable(w,imin,imax,2)
      ia      = iws_BeginTbody(w,itabv)-1
      do i = 1,nv
        w(ia+i) = vnodes(i)
      enddo
C--   Create coefficient table
      imin(1) = 1
      imax(1) = nu
      imin(2) = 1
      imax(2) = nv
      imin(3) = 0
      imax(3) = 3
      imin(4) = 0
      imax(4) = 3
      ifbcd   = iws_wtable(w,imin,imax,4)
C--   Create coefficient array
      imin(1) = 0
      imax(1) = 3
      imin(2) = 0
      imax(2) = 3
      iccij   = iws_wtable(w,imin,imax,2)
C--   Store local addresses in tag field
      ia = iws_IaFirstTag(w,iasp2)
      w(ia+ImarkS0) = dble(MarkSp0)
      w(ia+IdimSp0) = dble(2)
      w(ia+NedegU0) = dble(3)
      w(ia+NedegV0) = dble(3)
      w(ia+IsUtab0) = dble(itabu-iasp2)
      w(ia+NwUtab0) = dble(nu)
      w(ia+IsVtab0) = dble(itabv-iasp2)
      w(ia+NwVtab0) = dble(nv)
      w(ia+Nodesa0) = dble(nu*nv)
      w(ia+IsUser0) = dble(iuser-iasp2)
      w(ia+IsFBCD0) = dble(ifbcd-iasp2)
      w(ia+IsCCij0) = dble(iccij-iasp2)
C--   Store address of first spline
      iaR = iws_iaRoot()
      iat = iws_IaFirstTag(w,iaR)
      if(int(w(iat+IaSfirst0)).eq.0) w(iat+IaSfirst0) = dble(iasp2)
C--   Return address
      ispS2make     = iasp2

      return
      end

C     =====================================
      subroutine sspS2fill(w, iasp2, fvals)
C     =====================================

C--   Fill and construct 2-dimensional spline

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*), fvals(maxn0,*)
      dimension k4(0:4), imin(1), imax(1)

C--   Inline pointer function
      kk4(i,j,n,m) = k4(0)+k4(1)*i+k4(2)*j+k4(3)*n+k4(4)*m+iaFF
C--   Get addresses
      call sspGetIaTwoD(w,iasp2,iat,iau,nus,iav,nvs,iaFF,iaCC)
C--   Get pointer coefficients 4-dim
      ia0   = iws_IaKARRAY(w,iaFF)
      k4(0) = int(w(ia0))
      k4(1) = int(w(ia0+1))
      k4(2) = int(w(ia0+2))
      k4(3) = int(w(ia0+3))
      k4(4) = int(w(ia0+4))
C--   Base address of range tables (in 2nd dim of u and v tables)
      ianu  = iau+nus-1
      ianv  = iav+nvs-1
C--   Create temporary set
      itset = iws_newset(w)
      imin  = 1
      imax  = max(nus,nvs)
      itabf = iws_wtable(w,imin,imax,1)
      itabb = iws_wtable(w,imin,imax,1)
      itabc = iws_wtable(w,imin,imax,1)
      itabd = iws_wtable(w,imin,imax,1)
      iaf   = iws_BeginTbody(w,itabf)
      iab   = iws_BeginTbody(w,itabb)
      iac   = iws_BeginTbody(w,itabc)
      iad   = iws_BeginTbody(w,itabd)
C--   For each iu do spline in v
      do iu = 1,nus
C--     Fill f
        nv  = int(w(ianu+iu))
        if(nv.ne.0) then
          iafm1  = iaf-1
          do iv = 1,nv
            w(iafm1+iv) = fvals(iu,iv)
          enddo
C--       Make spline coefficients
          call smb_spline(nv,w(iav),w(iaf),w(iab),w(iac),w(iad))
C--       Store spline coefficients
          icf = kk4(iu,1,0,0)
          icb = kk4(iu,1,0,1)
          icc = kk4(iu,1,0,2)
          icd = kk4(iu,1,0,3)
          ivv = 0
          do iv = 1,nv
            ivm1       = iv-1
            w(icf+ivv) = w(iaf+ivm1)
            w(icb+ivv) = w(iab+ivm1)
            w(icc+ivv) = w(iac+ivm1)
            w(icd+ivv) = w(iad+ivm1)
            ivv        = ivv+k4(2)
          enddo
        endif
      enddo
C--   Spline tables in u for each iv
      do iv = 1,nvs
        nu  = int(w(ianv+iv))
        if(nu.ne.0) then
C--       Loop over j=f,b,c,d
          do j = 0,3
C--         Spline in u
            iaf = kk4(1,iv,0,j)
            iab = kk4(1,iv,1,j)
            iac = kk4(1,iv,2,j)
            iad = kk4(1,iv,3,j)
            call smb_spline(nu,w(iau),w(iaf),w(iab),w(iac),w(iad))
          enddo
        endif
      enddo
C--   Get rid of temporary tables
      call sws_wswipe(w,itset)

      return
      end

C     ==============================================================
      subroutine
     +      sspGetIaTwoD(w, ia, iat, iau, nus, iav, nvs, iaFF, iaCC)
C     ==============================================================

C--   Get table adresses

C--   w          (in): workspace
C--   ia         (in): address of 2-dim spline
C--   iat       (out): 1st word of tag field
C--   iau,nus   (out): first word of unode table and number of words
C--   iav,nuv   (out): first word of vnode table and number of words
C--   iaFF      (out): table addresses FBCD table
C--   iaCC      (out): table addresses Cij table

      implicit double precision(a-h,o-z)

      include 'splint.inc'

      dimension w(*)

      iat      = iws_IaFirstTag(w,ia)
      iatabu   = int(w(iat+IsUtab0)) + ia
      iatabv   = int(w(iat+IsVtab0)) + ia
      nus      = int(w(iat+NwUtab0))
      nvs      = int(w(iat+NwVtab0))
      iau      = iws_BeginTbody(w,iatabu)
      iav      = iws_BeginTbody(w,iatabv)
      iaFF     = int(w(iat+IsFBCD0)) + ia
      iaCC     = int(w(iat+IsCCij0)) + ia

      return
      end

C     ===================================
      subroutine sspRangeYT(w, ia, rscut)
C     ===================================

C--   Fill range tables nt and ny of 2-dim spline
C--   The range tables are stored in dim-2 of the u and v node tables

C--   w     (in): workspace
C--   ia    (in): address of 2-dim spline
C--   rscut (in): ln(rs*rs)

C--   For fixed y we have t < tmax(y) with tmax(y) = rscut-y
C--   For fixed t we have y < ymax(t) with ymax(t) = rscut-t

C--   The coordinates are called (y,t) or equivalently (u,v)

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

C--   Get addresses
      call sspGetIaTwoD(w,ia,iat,iau,nu,iav,nv,iaFF,iaCC)
C--   Base addresses: note that u and v tables are 2dim
      iayy = iau-1        !base address table y(iy)   node points
      iatm = iau+nu-1     !base address table itm(iy) to be filled
      iatt = iav-1        !base address table t(it)   node points
      iaym = iav+nv-1     !base address table iym(it) to be filled
C--   Go...
      if(rscut.eq.0.D0) then
C--     No cut set t-range to nv and y-range to nu
        do iy = 1,nu
          w(iatm+iy) = nv
        enddo
        do it = 1,nv
          w(iaym+it) = nu
        enddo
        w(iat+Nodesa0) = dble(nu*nv)
      else
C--     Initialise
        call smb_vfill(w(iaym+1),nv,0.D0)
C--     Set t-range for all y-nodes
        do iy = 1,nu-1
          yy = w(iayy+iy)
          tt = rscut-yy
          it = ispGetBin(tt,w(iav),nv)
          if(it.eq.-1) then
            itm = 0
          elseif(it.eq.-nv) then
            itm = nv
          else
            itm = min(it+1,nv)
          endif
          w(iatm+iy+1) = dble(itm)
          if(itm.ne.0) w(iaym+itm)  = dble(iy+1)
        enddo
        w(iatm+1) = w(iatm+2)
C--     Set y-range for all t-nodes
        iylast = int(w(iaym+nv))
        nodesa = iylast
        do it = nv-1,1,-1
          iy         = int(w(iaym+it))
          iynew      = max(iy,iylast)
          w(iaym+it) = dble(iynew)
          iylast     = iynew
          nodesa     = nodesa+iynew
        enddo
        w(iat+Nodesa0) = nodesa
      endif

      return
      end

C     ===============================================
      double precision function dspS2fun(w, ia, u, v)
C     ===============================================

C--   Compute 2-dim spline function

      implicit double  precision (a-h,o-z)
      logical lmb_gt

      include 'splint.inc'

      save iarem, iat, iau, nus, iav, nvs, iaFF, iaCC, k4
      dimension w(*), k4(0:4), cc(0:3,0:3), aa(0:3)
      data iarem/0/

C--   Inline pointer function
      kk4(i,j,n,m) = k4(0)+k4(1)*i+k4(2)*j+k4(3)*n+k4(4)*m+iaFF

      dspS2fun = 0.D0

C--   Get addresses and pointer coefficients
      if(ia.ne.iarem) then
        call sspGetIaTwoD(w,ia,iat,iau,nus,iav,nvs,iaFF,iaCC)
        ia0   = iws_IaKARRAY(w,iaFF)
        k4(0) = int(w(ia0))
        k4(1) = int(w(ia0+1))
        k4(2) = int(w(ia0+2))
        k4(3) = int(w(ia0+3))
        k4(4) = int(w(ia0+4))
        iarem = ia
      endif
C--   Get node bins
      iu = ispGetBin(u,w(iau),nus)
      iv = ispGetBin(v,w(iav),nvs)
C--   Check kinematic limit
      if(iu.gt.0 .and. iv.gt.0) then
        ianu = iau+nus
        ianv = iav+nvs
        jv   = min(iv+1,nvs)
        iuma = int(w(ianv+jv-1))
        umax = w(iau+iuma-1)
        ju   = min(iu,nus)
        ivma = int(w(ianu+ju-1))
        vmax = w(iav+ivma-1)
        if(lmb_gt(u,umax,-deps0) .or. lmb_gt(v,vmax,-deps0)) return
      endif
C--   Copy spline coefficients
      ian = kk4(iu,iv,0,0)
      do n = 0,3
        iam = ian
        do m = 0,3
          cc(n,m) = w(iam)
          iam     = iam+k4(4)
        enddo
        ian = ian+k4(3)
      enddo
C--   Set polynomial degree
      ndegu = 3
      ndegv = 3
      if(iu.lt.0) ndegu = int(w(iat+NedegU0))  !u out of range
      if(iv.lt.0) ndegv = int(w(iat+NedegV0))  !v out of range
C--   Go ...
      du = u - w(iau+abs(iu)-1)
      dv = v - w(iav+abs(iv)-1)
C--   Get v-coefficients from spline in u
      aa(0) = dspPol3(du,cc(0,0),ndegu)
      aa(1) = dspPol3(du,cc(0,1),ndegu)
      aa(2) = dspPol3(du,cc(0,2),ndegu)
      aa(3) = dspPol3(du,cc(0,3),ndegu)
C--   S2fun(u,v)
      dspS2fun = dspPol3(dv,aa,ndegv)
      
      return
      end

C=======================================================================
C===  2-dim integration from bin origin ================================
C=======================================================================

C     =========================================================
      double precision function dspBintYij(w, iasp, y1, y2, dt)
C     =========================================================

C--   Integration over y at fixed t:  Int_yi^y exp(-u) Pol3(u,t) du

C--   w     (in): workspace
C--   iasp  (in): address of 2-dim spline
C--   y1    (in): lower bin coordinate
C--   y2    (in): integration limit
C--   dt    (in): fixed value of local t coordinate

      implicit double precision (a-h,o-z)
      logical lmb_le

      include 'splint.inc'

      save iarem, iat, iau, nus, iav, nvs, iaFF, iaCC, k2
      dimension w(*), k2(0:2), eminu(0:3)
      data iarem/0/

C--   Inline pointer function
      kk2(n,m) = k2(0)+k2(1)*n+k2(2)*m+iaCC

      dy = y2-y1
C--   Empty domain
      if(lmb_le(dy,0.D0,-deps0)) then
        dspBintYij = 0.D0
        return
      endif
C--   Get addresses and pointer coefficients of Cnm array
      if(iasp.ne.iarem) then
        call sspGetIaTwoD(w,iasp,iat,iau,nus,iav,nvs,iaFF,iaCC)
        ia0   = iws_IaKARRAY(w,iaCC)
        k2(0) = int(w(ia0))
        k2(1) = int(w(ia0+1))
        k2(2) = int(w(ia0+2))
        iarem = iasp
      endif
C--   Get E-integrals
      call sspEminu(dy,eminu)
C--   Compute integral
      val = 0.D0
      dtm = 1.D0
      iam = kk2(0,0)
      
      do m = 0,3
        ian  = iam
        valn = 0.D0
        do n = 0,3
          valn = valn + w(ian)*eminu(n)
          ian  = ian + k2(1)
        enddo
        val = val + valn*dtm
        dtm = dtm*dt
        iam = iam + k2(2)
      enddo
C--   Done (almost)
      dspBintYij = exp(-y1)*val

      return
      end

C     =========================================================
      double precision function dspBintTij(w, iasp, t1, t2, dy)
C     =========================================================

C--   Integration over t at fixed y: Int_ti^t exp(v) Pol3(y,v) dv

C--   w     (in): workspace
C--   iasp  (in): address of 2-dim spline
C--   t1    (in): lower bin coordinate
C--   t2    (in): integration limit
C--   dy     (in): fixed value of local y coordinate

      implicit double precision (a-h,o-z)
      logical lmb_le

      include 'splint.inc'

      save iarem, iat, iau, nus, iav, nvs, iaFF, iaCC, k2
      dimension w(*), k2(0:2), eplus(0:3)
      data iarem/0/

C--   Inline pointer function
      kk2(n,m) = k2(0)+k2(1)*n+k2(2)*m+iaCC

      dt = t2-t1
C--   Empty domain
      if(lmb_le(dt,0.D0,-deps0)) then
        dspBintTij = 0.D0
        return
      endif
C--   Get addresses and pointer coefficients of Cnm array
      if(iasp.ne.iarem) then
        call sspGetIaTwoD(w,iasp,iat,iau,nus,iav,nvs,iaFF,iaCC)
        ia0   = iws_IaKARRAY(w,iaCC)
        k2(0) = int(w(ia0))
        k2(1) = int(w(ia0+1))
        k2(2) = int(w(ia0+2))
        iarem = iasp
      endif
C--   Get E-integrals
      call sspEplus(dy,eplus)
C--   Compute integral
      val = 0.D0
      dyn = 1.D0
      ian = kk2(0,0)
      
      do n = 0,3
        iam  = ian
        valm = 0.D0
        do m = 0,3
          valm = valm + w(iam)*eplus(m)
          iam  = iam + k2(2)
        enddo
        val = val + valm*dyn
        dyn = dyn*dy
        ian = ian + k2(1)
      enddo
C--   Done (almost)
      dspBintTij = exp(t1)*val

      return
      end

C     ==============================================================
      double precision function dspBintYTij(w, iasp, y1, y2, t1, t2)
C     ==============================================================

C--   Integral over y,t: Int_yi^y Int_ti^t exp(-u) exp(v) Pol3(u,v) dudv

C--   w     (in): workspace
C--   iasp  (in): address of 2-dim spline
C--   y1    (in): lower bin coordinate
C--   y2    (in): integration limit
C--   t1    (in): lower bin coordinate
C--   t2    (in): integration limit

      implicit double precision (a-h,o-z)
      logical lmb_le

      include 'splint.inc'

      save iarem, iat, iau, nus, iav, nvs, iaFF, iaCC, k2
      dimension w(*), k2(0:2), eminu(0:3), eplus(0:3)
      data iarem/0/

C--   Inline pointer function
      kk2(n,m) = k2(0)+k2(1)*n+k2(2)*m+iaCC

      dy = y2-y1
      dt = t2-t1
C--   Empty domain
      if(lmb_le(dy,0.D0,-deps0) .or. lmb_le(dt,0.D0,-deps0)) then
        dspBintYTij = 0.D0
        return
      endif
C--   Get addresses and pointer coefficients of Cnm array
      if(iasp.ne.iarem) then
        call sspGetIaTwoD(w,iasp,iat,iau,nus,iav,nvs,iaFF,iaCC)
        ia0   = iws_IaKARRAY(w,iaCC)
        k2(0) = int(w(ia0))
        k2(1) = int(w(ia0+1))
        k2(2) = int(w(ia0+2))
        iarem = iasp
      endif
C--   Get E-integrals
      call sspEminu(dy,eminu)
      call sspEplus(dt,eplus)
C--   Compute integral
      val = 0.D0
      ian = kk2(0,0)
      do n = 0,3
        valm = 0.D0
        iam  = ian
        do m = 0,3
          valm = valm + w(iam)*eplus(m)
          iam  = iam + k2(2)
        enddo
        val = val + valm*eminu(n)
        ian = ian + k2(1)
      enddo
C--   Done (almost)
      dspBintYTij = exp(-y1)*exp(t1)*val

      return
      end

C=======================================================================
C===  2-dim integration not from bin origin (with cut) =================
C=======================================================================

C     ==========================================================
      double precision function
     +      dspBintYYTT(w, iasp, iy, it, y1, y2, t1, t2, rs, np)
C     ==========================================================

C--   Compute integral over y-t rectangle within the bin

C--   w     (in): workspace
C--   iasp  (in): address of 2-dim spline
C--   iy    (in): node-point index in y < n_y
C--   it    (in): node-point index in t < n_t
C--   y1,2  (in): y-range will be brought in-range
C--   t1,2  (in): t-range will be brought in-range

      implicit double precision (a-h,o-z)
      logical lmb_le, lmb_ne, lmb_gt

      include 'splint.inc'

      dimension w(*)

      save iarem, iat, iau, nus, iav, nvs, iaFF, iaCC
      data iarem/0/

      common /pgaus/ y1gaus, t1gaus, scgaus, iagaus
      external dspGausFun

      dspBintYYTT = 0.D0
C--   Get addresses
      if(iasp.ne.iarem) then
        call sspGetIaTwoD(w,iasp,iat,iau,nus,iav,nvs,iaFF,iaCC)
        iagaus = iasp
        iarem  = iasp
      endif
C--   Lower and upper bin limits
      call sspBinLims(w,iasp,iy,it,ymi,yma,tmi,tma)
C--   Bring integration limits in range
      yy1 = max(y1,ymi)
      yy2 = min(y2,yma)
      tt1 = max(t1,tmi)
      tt2 = min(t2,tma)
C--   Check empty integration domain
      if(lmb_le(yy2,yy1,-deps0)) return
      if(lmb_le(tt2,tt1,-deps0)) return
C--   Check roots
      if(lmb_le(rs,0.D0,-deps0)) then
        sc = 0.D0
      else
        sc = log(rs*rs)
      endif

C--   See if sc limit crosses bin: 0 = below, 1 = cross, 2 = above
      icross = ispCrossSc(ymi, yma, tmi, tma, sc)

      if(icross.eq.2) then
C--     Bin above cut: set integral to zero
        dspBintYYTT = 0.D0

      elseif(icross.eq.0 .or. np.lt.2) then
C--     Bin below cut or ignore cut: do not re-parameterise spline
        call sspGetCoefs(w,iasp,iy,it,ymi,tmi)
        dspBintYYTT = dspBintYTij(w,iasp,ymi,yy2,tmi,tt2)
        if(lmb_ne(yy1,ymi,-deps0) .or. lmb_ne(tt1,tmi,-deps0)) then
C--       More to do because lower integration limit is not bin origin
          dspBintYYTT =
     +    dspBintYYTT + dspBintYTij(w,iasp,ymi,yy1,tmi,tt1)
     +                - dspBintYTij(w,iasp,ymi,yy1,tmi,tt2)
     +                - dspBintYTij(w,iasp,ymi,yy2,tmi,tt1)
        endif

      elseif(icross.eq.1 .and. np.ge.2) then
C--     Bin crossed by cut: re-parameterise spline to (yy1,tt1)
        call sspGetCoefs(w,iasp,iy,it,yy1,tt1)
C--     Snip integration domain
        call sspSnipSnip(yy1,yy2,tt1,tt2,sc,tta,ttb)
C--     Rectangle integration from tt1 to tta
        dspBintYYTT = dspBintYTij(w,iasp,yy1,yy2,tt1,tta)
C--     Gauss integration from tta to ttb
        if(lmb_gt(ttb,tta,-deps0)) then
          y1gaus      = yy1
          t1gaus      = tt1
          scgaus      = sc
          gint        = 0.D0
          if(np.eq.2) then
            gint = dmb_gaus2(dspGausFun,tta,ttb)
          elseif(np.eq.3) then
            gint = dmb_gaus3(dspGausFun,tta,ttb)
          elseif(np.eq.4) then
            gint = dmb_gaus4(dspGausFun,tta,ttb)
          else
            gint = dmb_gauss(dspGausFun,tta,ttb,1.D-7)
          endif
          dspBintYYTT = dspBintYYTT + gint
        endif

      else

        stop 'dspBINTYYTT: cant decide rectangle or Gauss integration'
        
      endif

      return
      end

C     =======================================
      double precision function dspGausFun(t)
C     =======================================

C--   Gauss quadrature input function

      implicit double precision (a-h,o-z)

      common /pgaus/ y1gaus, t1gaus, scgaus, iagaus

      include 'splint.inc'
      include 'spliws.inc'

      y2         = scgaus-t
      dt         = t-t1gaus
      dspGausFun = exp(t)*dspBintYij(w, iagaus, y1gaus, y2, dt)

      return
      end

C     =============================================================
      double precision function
     +            dspSpIntYT(w, iasp, y1, y2, t1, t2, rs, np, ierr)
C     =============================================================

C--   Compute integral over a collection of bins

C--   w     (in): workspace
C--   iasp  (in): address of 2-dim spline
C--   y1,2  (in): y-range
C--   t1,2  (in): t-range
C--   ierr (out): rs incompatible with rscut on iasp

      implicit double precision (a-h,o-z)

      include 'splint.inc'
      logical lmb_gt, lmb_le

      dimension w(*)

      save iarem, iat, iau, nus, iav, nvs, iaFF, iaCC
      data iarem/0/

      ierr       = 0
      dspSpIntYT = 0.D0
C--   Get addresses
      if(iasp.ne.iarem) then
        call sspGetIaTwoD(w,iasp,iat,iau,nus,iav,nvs,iaFF,iaCC)
        iarem = iasp
      endif
C--   Check rs
      rsc = w(iat+IrsCut0)
      if(lmb_gt(rsc,0.D0,-deps0)) then
        if(lmb_le(rs,0.D0,-deps0)) ierr = 1
        if(lmb_gt(rs,rsc ,-deps0)) ierr = 1
      endif
      if(ierr.eq.1) return

C--   Find node index limits
      iymi  = ispGetBin(y1,w(iau),nus)
      iyma  = ispGetBin(y2,w(iau),nus)
      itmi  = ispGetBin(t1,w(iav),nvs)
      itma  = ispGetBin(t2,w(iav),nvs)
C--   Catch out of range
      if(iymi.le.0 .or. iyma.le.0) stop
     +   ' SPLINT::dspSpIntYT: problem with y-limits out of range'
      if(itmi.le.0 .or. itma.le.0) stop
     +   ' SPLINT::dspSpIntYT: problem with t-limits out of range'
C--   Reset bin-number when ymax or tmax = upper spline limit
      if(iyma.eq.nus) iyma = iyma-1
      if(itma.eq.nvs) itma = itma-1
C--   Compute integral
      val = 0.D0
      do iy = iymi,iyma
        do it = itmi,itma
          val = val + dspBintYYTT(w,iasp,iy,it,y1,y2,t1,t2,rs,np)
        enddo
      enddo

      dspSpIntYT = val

      return
      end

C     ==================================================
      subroutine sspSnipSnip(y1, y2, t1, t2, sc, ta, tb)
C     ==================================================

C--   Find the cut intercepts ta (tb) with y2 (y1)
C--   The region t1-ta is for rectangle integration
C--   The region ta-tb is for Gauss quadrature
C--
C--   y1,2   (in): integration limits in y
C--   t1,2   (in): integration limits in t
C--   sc     (in): log(rs*rs)
C--   ta    (out): intercept with y2 (brought in range t1,t2)
C--   tb    (out): intercept with y1 (brought in range t1,t2)

      implicit double precision (a-h,o-z)
      logical lmb_le

      include 'splint.inc'

      if(lmb_le(sc,0.D0,-deps0)) then
C--     No cut
        ta = t2
        tb = t2
      else
C--   Intercepts with cut
        ta = max(sc-y2,t1)
        ta = min(ta,t2)
        tb = max(sc-y1,t1)
        tb = min(tb,t2)
      endif

*      if(ta.eq.t1) then
*        write(6,*) ' SnipSnip : Gauss only'
*      elseif(ta.eq.t2) then
*        write(6,*) ' SnipSnip : Rectangle only'
*      else
*        write(6,*) ' SnipSnip : Rectangle and Gauss'
*      endif

      return
      end

C     ===============================================
      subroutine sspGetCoefs(w, iasp, iy, it, y1, t1)
C     ===============================================

C--   Copy spline coefficients of bin (iy,it) to c(0:3,0:3) array in w
C--   If (y1,t1) are not lower bin edges then first transform to (y1,t1)

      implicit double precision (a-h,o-z)
      logical lmb_eq

      include 'splint.inc'

      dimension w(*), k2(0:2), k4(0:4), nfac(0:3)
      data nfac/1, 1, 2, 6/

      save iarem, iat, iau, nus, iav, nvs, iaFF, iaCC, k2, k4
      data iarem/0/
*      save ncopy, nrpar
*      data ncopy/0/
*      data nrpar/0/

C--   Inline pointer functions
      kk2(n,m)     = k2(0)+k2(1)*n+k2(2)*m+iaCC
      kk4(i,j,n,m) = k4(0)+k4(1)*i+k4(2)*j+k4(3)*n+k4(4)*m+iaFF

C--   Get addresses
      if(iasp.ne.iarem) then
        call sspGetIaTwoD(w,iasp,iat,iau,nus,iav,nvs,iaFF,iaCC)
        ia0   = iws_IaKARRAY(w,iaCC)
        k2(0) = int(w(ia0))
        k2(1) = int(w(ia0+1))
        k2(2) = int(w(ia0+2))
        ia0   = iws_IaKARRAY(w,iaFF)
        k4(0) = int(w(ia0))
        k4(1) = int(w(ia0+1))
        k4(2) = int(w(ia0+2))
        k4(3) = int(w(ia0+3))
        k4(4) = int(w(ia0+4))
        iarem = iasp
      endif
C--   Lower bin limits
      yi = w(iau+iy-1)
      ti = w(iav+it-1)
C--   Check if y1 and t1 are bin limits
      if(lmb_eq(y1,yi,-deps0) .and. lmb_eq(t1,ti,-deps0)) then
*        ncopy = ncopy+1
*        write(6,*),' just copy ', ncopy
C--     Copy coefficients
        iac = kk2(0,0)-1
        iam = kk4(iy,it,0,0)
        do m = 0,3
          ian = iam
          do n = 0,3
            iac = iac+1
            w(iac) = w(ian)
            ian = ian+k4(3)
           enddo
          iam = iam+k4(4)
        enddo
      else
C--     Re-parameterise to reference point (y1,t1)
*        nrpar = nrpar+1
*        write(6,*),' reparam   ', nrpar
        iac = kk2(0,0)-1
        do m = 0,3
          do n = 0,3
            iac = iac+1
            w(iac) = dspDerSp2(w,iasp,n,m,iy,it,y1,t1)/(nfac(n)*nfac(m))
          enddo
        enddo
      endif

      return
      end

C     ================================================================
      double precision function dspDerSp2(w, iasp, n, m, iy, it, y, t)
C     ================================================================

C--   Return d^(n+m)S / dy^n dt^m at y and t in bin (iy,it)
C--   w       (in): workspace
C--   iasp    (in): address of 2-dim spline
C--   n,m     (in): nm-th derivative
C--   iy,it   (in): bin index
C--   y,t     (in): interpolation point

      implicit double precision(a-h,o-z)

      include 'splint.inc'

      dimension w(*), k4(0:4), nfac(0:3), c(0:3), a(0:3)
      data nfac/1, 1, 2, 6/

      save iarem, iat, iau, nus, iav, nvs, iaFF, iaCC, k4, nfac
      data iarem/0/

C--   Inline pointer functions
      kk4(i,j,n,m) = k4(0)+k4(1)*i+k4(2)*j+k4(3)*n+k4(4)*m+iaFF

C--   Get addresses
      if(iasp.ne.iarem) then
        call sspGetIaTwoD(w,iasp,iat,iau,nus,iav,nvs,iaFF,iaCC)
        ia0   = iws_IaKARRAY(w,iaFF)
        k4(0) = int(w(ia0))
        k4(1) = int(w(ia0+1))
        k4(2) = int(w(ia0+2))
        k4(3) = int(w(ia0+3))
        k4(4) = int(w(ia0+4))
        iarem = iasp
      endif
C--   Local coordinates
      yi = w(iau+iy-1)
      dy = y-yi
      ti = w(iav+it-1)
      dt = t-ti
C--   Four interpolations in y
      iaj = kk4(iy,it,n,m)
      do j = 0,3-m
        iai = iaj
        do i = 0,3-n
          c(i) = w(iai)*nfac(i+n)/nfac(i)
          iai  = iai + k4(3)
        enddo
        a(j) = dspPol3(dy,c,3-n)*nfac(j+m)/nfac(j)
        iaj  = iaj + k4(4)
      enddo
C--   One interpolation in t
      dspDerSp2 = dspPol3(dt,a,3-m)
      
      return
      end

C=======================================================================
C===  Access to spline object  =========================================
C=======================================================================

C     =====================================
      integer function ispSplineType(w, ia)
C     =====================================

C--   Return spline dimension +-1,2  (0 = not a spline)

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

      ityp  = 0
      iobj  = iws_ObjectType(w,ia)        !2 = table set
      if(iobj.eq.2) then
        iatag = iws_IaFirstTag(w,ia)
        if(int(w(iatag+ImarkS0)).eq.MarkSp0) then
          ityp = int(w(iatag+IdimSp0))
        endif
      endif
      ispSplineType = ityp

      return
      end

C     ===================================
      integer function ispReadOnly(w, ia)
C     ===================================

C--   Return yes(1)/no(0) if spline is read-only

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

      iatag       = iws_IaFirstTag(w,ia)
      ispReadOnly = int(w(iatag+IrOnly0))

      return
      end

C     =================================================================
      subroutine sspSpLims(w, ia, nu, umi, uma, nv, vmi, vma, ndim, nb)
C     =================================================================

C--   Get node-point limits (not converted to x,q)

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

      iatag   = iws_IaFirstTag(w,ia)
      ndim    = int(w(iatag+IdimSp0))
      nb      = int(w(iatag+Nodesa0))
      nu      = int(w(iatag+NwUtab0))
      nv      = int(w(iatag+NwVtab0))

      iau     = int(w(iatag+IsUtab0)) + ia
      iu1     = iws_BeginTbody(w,iau)
      iu2     = iu1+nu-1
      umi     = w(iu1)
      uma     = w(iu2)

      if(nv.eq.0) then
        vmi   = 0.D0
        vma   = 0.D0
      else
        iav   = int(w(iatag+IsVtab0)) + ia
        iv1   = iws_BeginTbody(w,iav)
        iv2   = iv1+nv-1
        vmi   = w(iv1)
        vma   = w(iv2)
      endif

      return
      end

C=======================================================================
C==== Spline utilities ==============================================
C=======================================================================

C     =======================================
      integer function ispIaFromI(w, iasp, i)
C     =======================================

C--   Get address of spline userstore(i);  0 if i out of range

      implicit double precision (a-h,o-z)

      include 'splint.inc'

      dimension w(*)

      iatag = iws_IaFirstTag(w,iasp)

      if(i.lt.1 .or. i.gt.nusr0) then
        ispIaFromI = 0
      else
        iau        = int(w(iatag+IsUser0))+iasp
        ispIaFromI = iau + i - 1
      endif

      return
      end

C     ========================================
      integer function  ispIyFromY(w, ias2, y)
C     ========================================

C--   Get node bin index below y
C--
C--   w       (in): workspace
C--   ias2    (in): address of 2-dim spline
C--   y       (in): value of y
C--
C--   returns -1 if y below range and -ny if y above range

      implicit double precision(a-h,o-z)

      dimension w(*)

      save iarem, iau, nus
      data iarem/0/

C--   Get addresses
      if(ias2.ne.iarem) then
        call sspGetIaTwoD(w,ias2,iat,iau,nus,iav,nvs,iaFF,iaCC)
        iarem = ias2
      endif
      ispIyFromY  = ispGetBin(y,w(iau),nus)

      return
      end

C     ========================================
      integer function  ispItFromT(w, ias2, t)
C     ========================================

C--   Get node bin index below y
C--
C--   w       (in): workspace
C--   ias2    (in): address of 2-dim spline
C--   t       (in): value of t
C--
C--   returns -1 if t below range and -nt if t above range

      implicit double precision(a-h,o-z)

      dimension w(*)

      save iarem, iav, nvs
      data iarem/0/

C--   Get addresses
      if(ias2.ne.iarem) then
        call sspGetIaTwoD(w,ias2,iat,iau,nus,iav,nvs,iaFF,iaCC)
        iarem = ias2
      endif
      ispItFromT  = ispGetBin(t,w(iav),nvs)

      return
      end

C     ============================================
      logical function lspIsaFbin(w, ias2, iy, it)
C     ============================================

C--   True if (iy,it) is a filled bin

      implicit double precision (a-h,o-z)
      logical ify, ift

      dimension w(*)

      save iarem, iau, nus, iav, nvs
      data iarem/0/

C--   Get addresses
      if(ias2.ne.iarem) then
        call sspGetIaTwoD(w,ias2,iat,iau,nus,iav,nvs,iaFF,iaCC)
        iarem = ias2
      endif

      inu = iau+nus
      inv = iav+nvs
      ify = (iy .lt. int(w(inu+it-1)))
      ift = (it .lt. int(w(inv+iy-1)))
      if(ify.neqv.ift) stop 'lspISAFBIN: assignement problem'

      lspIsaFbin = ify

      return
      end
      
C     ======================================================
      subroutine sspBinLims(w, ias2, iy, it, y1, y2, t1, t2)
C     ======================================================

C--   Get bin boundaries

      implicit double precision (a-h,o-z)

      dimension w(*)

      save iarem, iau, nus, iav, nvs
      data iarem/0/

C--   Get addresses
      if(ias2.ne.iarem) then
        call sspGetIaTwoD(w,ias2,iat,iau,nus,iav,nvs,iaFF,iaCC)
        iarem = ias2
      endif

      if(iy.lt.1 .or. iy.ge.nus) stop 'sspBINLIMS: iy out of range'
      if(it.lt.1 .or. it.ge.nvs) stop 'sspBINLIMS: it out of range'

      y1 = w(iau+iy-1)
      y2 = w(iau+iy)
      t1 = w(iav+it-1)
      t2 = w(iav+it)

      return
      end

C     ===============================================
      integer function ispCrossSc(y1, y2, t1, t2, sc)
C     ===============================================

C--   0 = below cut, 1 = crossed by cut, 2 above cut
C--
C--   y1,2    (in): bin limits in y
C--   t1,2    (in): bin limits in t
C--   sc      (in): log roots^2; 0=no cut
C--
C--   Returns 0 if no cut

      implicit double precision(a-h,o-z)
      logical lmb_ge, lmb_le

      include 'splint.inc'

      if(lmb_le(sc,0.D0,-deps0)) then
        ispCrossSc = 0
      else
        if(lmb_le(t2+y2,sc,-deps0)) then
          ispCrossSc = 0
        elseif(lmb_ge(t1+y1,sc,-deps0)) then
          ispCrossSc = 2
        else
          ispCrossSc = 1
        endif
      endif
        
      return
      end

C     ===============================================
      double precision function dspRsMax(w, ias2, sc)
C     ===============================================

C--   Return root-s max cut, for given rs cut
C--   w       (in): workspace
C--   ias2    (in): address of 2-dim spline
C--   sc      (in): log roots^2; 0=no cut

      implicit double precision(a-h,o-z)
      logical lmb_le

      include 'splint.inc'

      dimension w(*)

      save iarem, iau, nus, iav, nvs
      data iarem/0/

C--   Get addresses
      if(ias2.ne.iarem) then
        call sspGetIaTwoD(w,ias2,iat,iau,nus,iav,nvs,iaFF,iaCC)
        iarem = ias2
      endif

      if(lmb_le(sc,0.D0,-deps0)) then
        dspRsMax = 0.D0
      else
        rsm = 0.D0
        do it = 1,nvs-1
          do iy = 1,nus-1
            call sspBinLims(w, ias2, iy, it, y1, y2, t1, t2)
            icross = ispCrossSc(y1, y2, t1, t2, sc)
            if(icross.eq.1) then
              rsm = max(rsm,y2+t2)
            endif
          enddo
        enddo
        rsm      = sqrt(exp(rsm))
        dspRsMax = dble(int(rsm)+1)
      endif

      return
      end

C=======================================================================
C==== Nonspline utilities ==============================================
C=======================================================================

C     =============================================
      double precision function dspPol3(x, coef, n)
C     =============================================

C--   Horner to compute f = c0 + c1*x + c2*x^2 + c3*x^3

      implicit double precision (a-h,o-z)

      dimension coef(0:3)

      sum = coef(n)
      do i = n-1,0,-1
        sum = coef(i) + x*sum
      enddo

      dspPol3 = sum

      return
      end

C     ===================================
      integer function ispGetBin(u, x, n)
C     ===================================

C--   Find bin in array x for argument u
C--   Returns -1 for u < xmin and -n for u > xmax
C--   Binary search if u is not in the same bin as the previous call.
C--   Code adapted from smb_seval in MBUTIL
C--   Soft comparisons make search stable against rouding errors

      implicit double precision (a-h,o-z)
      logical lmb_lt, lmb_ge, lmb_le

      include 'splint.inc'

      dimension x(*)

      save i
      data i/1/

C--   Quick check
      if ( i .ge. n )      i = 1
      if ( lmb_lt(u,x(i)  ,-deps0) )   go to 10
      if ( lmb_lt(u,x(i+1),-deps0) )   go to 30
C--   Binary search
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( lmb_lt(u,x(k),-deps0) ) j = k
      if ( lmb_ge(u,x(k),-deps0) ) i = k
      if ( j .gt. i+1 ) go to 20
C--   Bin found
   30 continue
      if( lmb_ge(u,x(1),-deps0) .and. lmb_le(u,x(n),-deps0)) then
        ispGetBin =  i
      else
        ispGetBin = -i
      endif

      return
      end

C     =============================
      subroutine sspEplus(x, Eplus)
C     =============================

C--   Eplus(n) = Int_0^x z^n exp(z) dz for n = 0,..,3

      implicit double precision (a-h,o-z)

      dimension Eplus(0:3)

C--   Check x.gt.0
      if(x.lt.0.D0) stop ' SPLINT::sspEplus: x < 0'

      Eplus(0) = exp(x) - 1.D0
      do i = 1,3
        Eplus(i) = x**i * exp(x) - i*Eplus(i-1)
      enddo

      return
      end

C     =============================
      subroutine sspEminu(x, Eminu)
C     =============================

C--   Eminu(n) = Int_0^x z^n exp(-z) dz for n = 0,..,3

      implicit double precision (a-h,o-z)

      dimension Eminu(0:3)

C--   Check x.gt.0
      if(x.lt.0.D0) stop ' SPLINT::sspEminu: x < 0'

      Eminu(0) = 1.D0 - exp(-x)
      do i = 1,3
        Eminu(i) = i*Eminu(i-1) - x**i * exp(-x)
      enddo

      return
      end
