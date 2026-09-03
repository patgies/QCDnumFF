
C--   This is the file srcEvolveOld.f containing the old evolution routines

C--   subroutine sqcPdfMat
C--
C--   subroutine sqcPdIdef(tmatpq,ierr)
C--   subroutine sqcAllInp(idg,func)
C--   subroutine sqcEfrmP(pval,eval)
C--   double precision function dqcEifrmP(i,pval)

C--   subroutine sqcEvolFG_old(ityp,jset,func,def,iq0,epsi,nfheavy,ierr)
C--   subroutine sqcEvolve(itype,jset,iord,it0,it1,it2,epsm,ierr)

C--   subroutine sqcNStart(itype,idout,idin,ids,idg,iyg,iord,dlam,it0)
C--   subroutine sqcGridns(itype,
C--              ipdf,ids,idg,ityp,iyg,iord,it0,it1,it2,eps,ierr)
C--   subroutine sqcNSevnf(itype,ipdf,ityp,iyg,iord,nf,iw1,iw2)
C--   subroutine sqcNSder(itype,ipdf,ityp,iyg,iord,nf,iwa,iwader)
C--   subroutine sqcNSjup(itype,idq,ids,idg,iyg,iord,dlam,ny,iwin,iwout)
C--   subroutine sqcNSjdn(itype,idq,ids,idg,iyg,iord,dlam,ny,iwin,iwout)
C--   subroutine sqcNSStoreStart(itype,ipdf,iy1,iy2,it0)
C--   subroutine sqcNSNewStart(itype,ipdf,iy1,iy2,it0,epsi)
C--   subroutine sqcNSRestoreStart(itype,ipdf,iy1,iy2,it0)
C--
C--   subroutine sqcGridsg(itype,
C--                        idf,idg,iyg,iord,it0,it1,it2,eps,ierr)
C--   subroutine sqcSGevnf(itype,idf,idg,iyg,iord,nf,iw1,iw2)
C--   subroutine sqcSGder(itype,idf,idg,iyg,iord,nf,iwa,iwader)
C--   subroutine sqcSGjup(itype,ids,idg,iyg,iord,ny,iwin,iwout)
C--   subroutine sqcSGjdn(itype,ids,idg,iyg,iord,ny,iwin,iwout)
C--   subroutine sqcSGStoreStart(itype,ids,idg,iy1,iy2,it0)
C--   subroutine sqcSGNewStart(itype,ids,idg,iy1,iy2,it0,epsi)
C--   subroutine sqcSGRestoreStart(itype,ids,idg,iy1,iy2,it0)
C--   double precision function dqcGetEps(itype,id,ny,it)
C--   subroutine sqcEvLims(it0,it1,it2,
C--              iwu1,iwu2,nflu,nup,iwd1,iwd2,nfld,ndn,ibl,ibu)

C     ====================
      subroutine sqcPdfMat
C     ====================

C--   Setup the matrices which relate the flavour basis |q> to the
C--   to the evolution basis |e>
C--
C--   |ei> = tmateq7(i,j) |qj>
C--   |qj> = tmatqe7(i,j) |ei>    NB: |q> = |q +- qbar>
C--
C--            1  2  3  4  5  6  7  8  9 10 11 12
C--   |ei>     s  2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-
C--   |qi>     d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--   Called by sqc_qcinit

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension imateq(6,6)

      data imateq /
C--      d  u  s  c  b  t
C--      1  2  3  4  5  6
     +   1, 1, 1, 1, 1, 1, !si   =  1
     +   1,-1, 0, 0, 0, 0, !ns1  =  2
     +   1, 1,-2, 0, 0, 0, !ns2  =  3
     +   1, 1, 1,-3, 0, 0, !ns3  =  4
     +   1, 1, 1, 1,-4, 0, !ns4  =  5
     +   1, 1, 1, 1, 1,-5/ !ns5  =  6

C--   Transform imateq to math notation (swap indices) -> umateq7
      do nf = 3,6
        do i = 1,6
          do j = 1,6
            umateq7(i,j,nf) = 0.D0           !initialise
          enddo
          umateq7(i,i,nf) = 1.D0             !preset diagonal
        enddo
        do i = 1,nf
          do j = 1,nf
            umateq7(i,j,nf) = imateq(j,i)    !copy nf*nf submatrix
          enddo
        enddo
      enddo
C--   Invert umateq7
      do nf = 3,6
        do i = 1,6
          do j = 1,6
            umatqe7(i,j,nf) = 0.D0           !initialise
          enddo
        enddo
C--     Invert the nf*nf submatrix
        call sqcOrtInv(umateq7(1,1,nf),umatqe7(1,1,nf),6,nf)
C--     Set diagonal
        do i = nf+1,6
            umatqe7(i,i,nf) = 1.D0
        enddo
      enddo
C--   Umat depends on nf but for tmat we take nf = 6
C--   Fill matrix tmateq7 and tmatqe7
      do i = 1,6
        do j = 1,6
          tmateq7(i,j)     = umateq7(i,j,6)
          tmateq7(i,j+6)   = 0.D0
          tmateq7(i+6,j)   = 0.D0
          tmateq7(i+6,j+6) = umateq7(i,j,6)
          tmatqe7(i,j)     = umatqe7(i,j,6)
          tmatqe7(i,j+6)   = 0.D0
          tmatqe7(i+6,j)   = 0.D0
          tmatqe7(i+6,j+6) = umatqe7(i,j,6)
        enddo
      enddo

      return
      end

C     =================================
      subroutine sqcPdIdef(tmatpq,ierr)
C     =================================

C--   The set of input distns defined (by the user) on the q+-
C--   basis is stored in the matrix tmatpq(i,j), where i = 1,...,12 is
C--   the user's input pdf index and j = 1,...12 is the q+- index:
C--            1  2  3  4  5  6  7  8  9 10 11 12
C--           d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--   Here are the various transformations performed in QCDNUM:
C--
C--   |pi> = tmatpq(i,j) |qj>  pdfs written on the q+- basis (input)
C--   |qi> = tmatqp(i,j) |pj>  inverse of the above
C--   |qi> = tmatqe7(i,j) |ej>  q->e transformation matrix (from sqcPdfMat)
C--   |ei> = tmateq7(i,j) |qj>  inverse of above (from sqcPdfMat)
C--   |pi> = tmatpe7(i,j) |ej>  pdfs written on the si/ns basis
C--   |ei> = tmatep7(i,j) |pj>  si/ns basis as lin comb of input pdfs
C--
C--   Given the matrices tmatpq, tmatqe7 and tmateq7 the routine calculates
C--   the matrices tmatpe7 and tmatep7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension tmatpq(12,12),tmatqp(12,12)
      dimension iwork(12)

C--   Invert tmatpq
      do i = 1,12
        do j = 1,12
          tmatqp(i,j)   = tmatpq(i,j)
        enddo
      enddo
      call smb_dminv(12,tmatqp,12,iwork,ierr)
C--   Oh, lala .... that doesnt look good
      if(ierr.ne.0) return

C--   tmatpe7 = tmatpq*tmatqe7
      do i = 1,12
        do j = 1,12
          sum = 0.D0
          do k = 1,12
            sum = sum + tmatpq(i,k)*tmatqe7(k,j)
          enddo
          tmatpe7(i,j) = sum
        enddo
      enddo
C--   tmatep7 = tmateq7*tmatqp
      do i = 1,12
        do j = 1,12
          sum = 0.D0
          do k = 1,12
            sum = sum + tmateq7(i,k)*tmatqp(k,j)
          enddo
          tmatep7(i,j) = sum
        enddo
      enddo

      return
      end

C     ==============================
      subroutine sqcAllInp(idg,func)
C     ==============================

C--   Calculate 2nf+1 pdfs provided by the user in func(j,x) and
C--   store these in      |p> (j=0,...,12)
C--   Then transform to   |e> (i=0,...,12)
C--   Put the transformed pdfs in stor7(iy,it=0,i), i = 0,...,12
C--   NB: pdfs stored in bin it = 0 serve as evolution start values
C--
C--   idg   (in)  gluon table identifier in global format

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension pval(0:12), eval(0:12)

C--   Initialize
      nfmax = max(abs(nfix6),3)
      do i = 0,12
        pval(i) = 0.D0
        eval(i) = 0.D0
      enddo

C--   Loop over y gridpoints
      do iy = 1,nyy2(0)
        y  = ygrid2(iy)
        x  = exp(-y)
        do j = 0,2*nfmax
          pval(j) = func(j,x)
        enddo
C--     Transform
        call sqcEfrmP(pval,eval)
C--     Store
        do id = 0,12
          iadr = iqcG5ijk(stor7,iy,0,idg+id)
          stor7(iadr) = eval(id)
        enddo
      enddo

      return
      end

C     ==============================
      subroutine sqcEfrmP(pval,eval)
C     ==============================

C--   Transform 13 parton values |p> (j=0,...,12)
C--   to the evolution basis     |e> (i=0,...,12)
C--
C--   Input:  pval(0:12)  input values  |p>
C--   Output: eval(0:12)  output values |e>

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension pval(0:12), eval(0:12)

      do i = 0,12
        eval(i) = dqcEifrmP(i,pval)
      enddo

      return
      end

C     ===========================================
      double precision function dqcEifrmP(i,pval)
C     ===========================================

C--   Transform 13 parton values            |pj> (j=0,...,12)
C--   to one element of the evolution basis |ei> (i=0,...,12)
C--
C--   Input:  i          = index of |ei> basis element (i = 0 = gluon)
C--           pval(0:12) = input vector |p>  (j = 0 = gluon)
C--
C--   Output: dqcEifrmP  = value of evolution basis element |ei>

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension pval(0:12)

      if(i.eq.0) then
C--     Gluon
        dqcEifrmP = pval(0)
      else
C--     |ei>  = tmatep7(i,j) |pj>
        sum = 0.D0
        do j = 1,12
          sum = sum + tmatep7(i,j)*pval(j)
        enddo
        dqcEifrmP = sum
      endif

      return
      end


C     ==================================================================
      subroutine sqcEvolFG_old(ityp,jset,func,def,iq0,epsi,nfheavy,ierr)
C     ==================================================================

C--   Evolve all pdfs: bits and pieces of old user interface.

C--   ierr = 1  Start point not inside grid or cuts
C--          2  At least one evolution limit below alphas cut
C--          3  Input pdfs not linearly independent

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension def(-6:6,12), pdef(12,12)

      ierr = 0

C--   Set cuts in iz make sure to include iq0
      izmic2 = izfit5(-min(itmic2,iq0))
      izmac2 = izfit5( max(itmac2,iq0))

C--   Cut range
      iq1 = itfiz5(izmic2)
      iq2 = itfiz5(izmac2)
      if(iq0.lt.iq1 .or. iq0.gt.iq2) then
        ierr = 1
        return
      endif
      
C--   Initialize pdef
      do i = 1,12
        do j = 1,12
          pdef(i,j) = 0.D0
        enddo
      enddo
C--   Predefine heavy quark qplus  in case nfmax < 4
      pdef( 7, 4) = 1.D0      !cplus
      pdef( 8,10) = 1.D0      !cmin
      pdef( 9, 5) = 1.D0      !bplus
      pdef(10,11) = 1.D0      !bmin
      pdef(11, 6) = 1.D0      !tplus
      pdef(12,12) = 1.D0      !tmin
C--   Build q+- matrix (pdef) from user input (def)
C--   Note index swap def(iflavor,ipdf) --> pdef(ipdf,iflavor)
      nfmax = max(abs(nfix6),3)
      do ipdf = 1,2*nfmax
        do j = 1,6
          pdef(ipdf,j)   = 0.5*(def(j,ipdf)+def(-j,ipdf))
          pdef(ipdf,j+6) = 0.5*(def(j,ipdf)-def(-j,ipdf))
        enddo
      enddo 

C--   Now tell qcdnum about these definitions
      call sqcPdIdef(pdef,jerr)
C--   Well that did not go OK...
      if(jerr.ne.0) then
        ierr = 3
        return
      endif

C--   Global identifier of gluon
      idg = iqcIdPdfLtoG(jset,0)
C--   Enter input distributions
      call sqcAllInp(idg,func)

C--   Off we go...
      call sqcEvolve(ityp,jset,iord6,iq0,iq1,iq2,epsi,jerr)
      if(jerr.eq.1) ierr = 2

C--   Afterburner, set all heavy quark basis pdfs to zero below treshold
      do iz = 1,nzz2
        nf  = itfiz5(-iz)
        do i = nf+1,6
          id = iqcIdPdfLtoG(jset,i)
          call sqcPsetjj(id,iz,0.D0)
          id = iqcIdPdfLtoG(jset,i+6)
          call sqcPsetjj(id,iz,0.D0)
        enddo
      enddo

      nfheavy = 6

C--   Re-set cuts in iz (which may have been stretched to include iq0)
      izmic2 = izfit5(-itmic2)
      izmac2 = izfit5( itmac2)

C--   Validate the pdfs
      do id = 0,12
        idglobal = iqcIdPdfLtoG(jset,id)
        call sqcValidate(stor7,idglobal)
      enddo

      return
      end

C     ================================================================
C     Steering routine to evolve all pdfs in the FFNS and VFNS
C     ================================================================

C     ===========================================================
      subroutine sqcEvolve(itype,jset,iord,it0,it1,it2,epsm,ierr)
C     ===========================================================

C--   Steering routine to evolve all pdfs
C--   This routine does not set the starting values
C--   itype  (in) : 1=unpol, 2=pol, 3=timelike
C--   jset   (in) : output pdf set identifier
C--   epsm  (out) : max deviation quad - lin at midpoints
C--   ierr  (out) : 1 = no alphas available at it1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension dlam(4:6)
      dimension lambda(4:6,12) !heavy quark weights in NNLO jumps
C--          nf =   4  5  6
      data lambda / 1, 1, 1,   ! 1 = e1+ singlet
     +              0, 0, 0,   ! 2 = e2+ updown ns+
     +              0, 0, 0,   ! 3 = e3+ strange ns+
     +             -3, 0, 0,   ! 4 = e4+ charm ns+
     +              0,-4, 0,   ! 5 = e5+ bottom ns+
     +              0, 0,-5,   ! 6 = e6+ top ns+
     +              0, 0, 0,   ! 7 = e1- valence
     +              0, 0, 0,   ! 8 = e2- updown ns-
     +              0, 0, 0,   ! 9 = e3- strange ns- 
     +              0, 0, 0,   !10 = e4- charm ns- 
     +              0, 0, 0,   !11 = e5- bottom ns- 
     +              0, 0, 0 /  !12 = e6- top ns- 

C--   VFNS checks
      if(abs(nfix6).eq.0) then
C--     To calculate ns jumps in NNLO we need the singlet at nf-1. This 
C--     implies that for charm a 3-flavor singlet is used. A 3-flavor
C--     singlet is only guaranteed to be available when the charm threshold
C--     is at the second grid point or larger.
        if(itchm2.lt.2 .and. iord.ge.2) stop
     +     'sqcEvolve: itchm2 .lt. 2 not allowed at NLO, NNLO ---> STOP'
C--     Furthermore it0 must be below the charm threshold
        if(it0.ge.itchm2) stop 
     +     'sqcEvolve: it0 at or above itchm2 ---> STOP'
      endif

C--   Pdf indices in main storage
      idg    =  0
      idf    =  1
      iudpl  =  2
      isspl  =  3
      ichpl  =  4
      ibopl  =  5
      itopl  =  6
      idv    =  7
      iudmi  =  8
      issmi  =  9
      ichmi  = 10
      ibomi  = 11
      itomi  = 12
C--   Global pdf identifiers
      idgg    =  iqcIdPdfLtoG(jset, 0)
      idfg    =  iqcIdPdfLtoG(jset, 1)
      iudplg  =  iqcIdPdfLtoG(jset, 2)
      issplg  =  iqcIdPdfLtoG(jset, 3)
      ichplg  =  iqcIdPdfLtoG(jset, 4)
      iboplg  =  iqcIdPdfLtoG(jset, 5)
      itoplg  =  iqcIdPdfLtoG(jset, 6)
      idvg    =  iqcIdPdfLtoG(jset, 7)
      iudmig  =  iqcIdPdfLtoG(jset, 8)
      issmig  =  iqcIdPdfLtoG(jset, 9)
      ichmig  =  iqcIdPdfLtoG(jset,10)
      ibomig  =  iqcIdPdfLtoG(jset,11)
      itomig  =  iqcIdPdfLtoG(jset,12)

      epsm = 0.D0

C--   iz range
      iz1 = izfit5(it1)
      iz2 = izfit5(it2)

C--   Loop over subgrids
      do ig = 1,nyg2
      
C--     Upper y index in subgrid      
        nyg = iqcIyMaxG(iymac2,ig)

C--  A. Singlet/gluon evolution
C--     -----------------------
        iftmp  = -4  !singlet subgrid identifier
        igtmp  = -3  !gluon   subgrid identifier
        iftmpg = iqcIdPdfLtoG(jset,iftmp)
        igtmpg = iqcIdPdfLtoG(jset,igtmp)
C--     Copy starting values at it=0 from G0 to Gi
        call sqcG0toGi(idfg,iftmpg,ig,nyg,0)
        call sqcG0toGi(idgg,igtmpg,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(iftmpg,iftmpg,nyg,0,0)
        call sqcGiFtoA(igtmpg,igtmpg,nyg,0,0)
C--     Evolve
        call sqcGridsg(itype,
     +                 iftmp,igtmp,ig,iord,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)

C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(iftmpg,ig,nyg,iz1,iz2,idfg)
        call sqcAitoF0(igtmpg,ig,nyg,iz1,iz2,idgg)

C--  B. NS+ evolution
C--     -------------
        ityp   = 1      !NS+

C--  B1 ud+ evolution
C--     -------------
        ipdf    = iudpl  !ud+
        ipdfg   = iudplg
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg  = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,2)
        dlam(5) = lambda(5,2)
        dlam(6) = lambda(6,2)
C--     Copy startvalues at it=0 from G0 to Gi
        call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)

C--  B2 s+ evolution
C--     ------------
        ipdf    = isspl  !s+
        ipdfg   = issplg
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,3)
        dlam(5) = lambda(5,3)
        dlam(6) = lambda(6,3)
C--     Copy startvalues at it=0 from G0 to Gi
        call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)

C--  B3 c+ evolution
C--     ------------
        itc     = max(it1,itchm2)
        ipdf    = ichpl  !c+
        ipdfg   = ichplg
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg  = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,4)
        dlam(5) = lambda(5,4)
        dlam(6) = lambda(6,4)
        if(abs(nfix6).eq.0) then
C--       VFNS: Copy singlet to Gi
          call sqcNStart(itype,
     +                   itemp,iftmp,iftmp,igtmp,ig,iord,dlam(4),itc)
C--       VFNS: Evolve
          if(itc.lt.it2) then
            call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,itc,itc,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(abs(nfix6).ge.4) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy singlet to Gi
          call sqcPdfCop(iftmpg,itempg)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)

C--  B4 b+ evolution
C--     ------------
        itb     = max(it1,itbot2)
        ipdf    = ibopl  !b+
        ipdfg   = iboplg
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg  = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,5)
        dlam(5) = lambda(5,5)
        dlam(6) = lambda(6,5)
        if(abs(nfix6).eq.0) then
C--       VFNS: Copy singlet to Gi
          call sqcNStart(itype,
     +                   itemp,iftmp,iftmp,igtmp,ig,iord,dlam(5),itb)
C--       VFNS: Evolve
          if(itb.lt.it2) then
            call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,itb,itb,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(abs(nfix6).ge.5) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          itempg = iqcIdPdfLtoG(jset,itemp)
          call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy singlet to Gi
          call sqcPdfCop(iftmpg,itempg)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)

C--  B5 t+ evolution
C--     ------------
        itt     = max(it1,ittop2)
        ipdf    = itopl  !t+
        ipdfg   = itoplg
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg  = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,6)
        dlam(5) = lambda(5,6)
        dlam(6) = lambda(6,6)
        if(abs(nfix6).eq.0) then
C--       VFNS: Copy singlet to Gi
          call sqcNStart(itype,
     +                   itemp,iftmp,iftmp,igtmp,ig,iord,dlam(6),itt)
C--       VFNS: Evolve
          if(itt.lt.it2) then
            call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,itt,itt,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(abs(nfix6).ge.6) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          itempg = iqcIdPdfLtoG(jset,itemp)
          call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy singlet to Gi
          call sqcPdfCop(iftmpg,itempg)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)

C--  C. NSV evolution
C--     -------------
        ityp    = 3      !NSV
        ipdf    = idv    !valence
        ipdfg   = idvg
        ivtmp   = -3     !subgrid identifier (temporary storage)
        ivtmpg = iqcIdPdfLtoG(jset,ivtmp)
        dlam(4) = lambda(4,7)
        dlam(5) = lambda(5,7)
        dlam(6) = lambda(6,7)
C--     Copy startvalues at it=0 from G0 to Gi
        call sqcG0toGi(ipdfg,ivtmpg,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(ivtmpg,ivtmpg,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +        ivtmp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(ivtmpg,ig,nyg,iz1,iz2,ipdfg)

C--  D. NS- evolution
C--     -------------
        ityp   = 2      !NS-

C--  D1 ud- evolution
C--     -------------
        ipdf    = iudmi  !ud-
        ipdfg   = iudmig
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg  = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,8)
        dlam(5) = lambda(5,8)
        dlam(6) = lambda(6,8)
C--     Copy startvalues at it=0 from G0 to Gi
        call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +        itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)

C--  D2 s- evolution
C--     ------------
        ipdf    = issmi  !s-
        ipdfg   = issmig
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,9)
        dlam(5) = lambda(5,9)
        dlam(6) = lambda(6,9)
C--     Copy startvalues from G0 to Gi
        call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +        itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)

C--  D3 c- evolution (if VFNS then only in NLLO)
C--     ----------------------------------------
        itc     = max(it1,itchm2)
        ipdf    = ichmi  !c-
        ipdfg   = ichmig
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,10)
        dlam(5) = lambda(5,10)
        dlam(6) = lambda(6,10)
        if(abs(nfix6).eq.0) then
C--       VFNS: Copy valence to Gi
          call sqcNStart(itype,itemp,ivtmp,0,0,ig,iord,dlam(4),itc)
C--       VFNS: Evolve
          if(itc.lt.it2 .and. iord.eq.3) then
            call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,itc,itc,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(abs(nfix6).ge.4) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          itempg = iqcIdPdfLtoG(jset,itemp)
          call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy valence to Gi
          call sqcPdfCop(ivtmpg,itempg)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)

C--  D4 b- evolution (if VFNS then only in NLLO)
C--     ----------------------------------------
        itb     = max(it1,itbot2)
        ipdf    = ibomi  !b-
        ipdfg   = ibomig
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg  = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,11)
        dlam(5) = lambda(5,11)
        dlam(6) = lambda(6,11)
        if(abs(nfix6).eq.0) then
C--       VFNS: Copy valence to Gi
          call sqcNStart(itype,itemp,ivtmp,0,0,ig,iord,dlam(5),itb)
C--       VFNS: Evolve
          if(itb.lt.it2 .and. iord.eq.3) then
            call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,itb,itb,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(abs(nfix6).ge.5) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy valence to Gi
          call sqcPdfCop(ivtmpg,itempg)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)
        
C--  D5 t- evolution (if VFNS then only in NLLO)
C--     ----------------------------------------
        itt     = max(it1,ittop2)
        ipdf    = itomi  !t-
        ipdfg   = itomig
        itemp   = -2     !subgrid identifier (temporary storage)
        itempg  = iqcIdPdfLtoG(jset,itemp)
        dlam(4) = lambda(4,12)
        dlam(5) = lambda(5,12)
        dlam(6) = lambda(6,12)
        if(abs(nfix6).eq.0) then
C--       VFNS: Copy valence to Gi
          call sqcNStart(itype,itemp,ivtmp,0,0,ig,iord,dlam(6),itt)
C--       VFNS: Evolve
          if(itt.lt.it2 .and. iord.eq.3) then
            call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,itt,itt,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(abs(nfix6).ge.6) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(ipdfg,itempg,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itempg,itempg,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy valence to Gi
          call sqcPdfCop(ivtmpg,itempg)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itempg,ig,nyg,iz1,iz2,ipdfg)
        
C--   End of loop over subgrids
      enddo

      return
      end

C     ================================================================
C     Nonsinglet evolution routines
C     ================================================================

C     ================================================================
      subroutine sqcNStart(itype,idout,idin,ids,idg,iyg,iord,dlam,it0)
C     ================================================================

C--   Copy singlet or valence and set proper startvalue for VFNS heavy quarks

C--   itype    (in)    1=unpol, 2=pol, 3=timelike
C--   idout    (in)    output nonsinglet table
C--   idin     (in)    input singlet or valence table
C--   ids      (in)    singlet table
C--   idg      (in)    gluon table
C--   iyg      (in)    subgrid index
C--   iord     (in)    1 = LO, 2 = NLO, 3 = NNLO
C--   dlam     (in)    heavy flavor weight
C--   it0      (in)    starting point

C--   Remark: at this point everything is in terms of A values

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

C--   Global identifiers
      igin  = iqcIdPdfLtoG(itype,idin)
      igout = iqcIdPdfLtoG(itype,idout)
C--   Copy input to output table (input is singlet or valence)
      call sqcPdfCop(igin,igout)
C--   Find out in which region it0 falls
      ibin0 = 0
      do i = 1,ntsubg6
        if(it0.eq.it1sub6(i)) ibin0 = i
      enddo
C--   Bin not found, dont set startvalue
      if(ibin0.eq.0) return
C--   Calculate iz0; izmin-itmin gives the offset between iz and it
      iz0 = it0 + iz1sub6(ibin0)-it1sub6(ibin0)
C--   Apply jumps
      if((itype.eq.1 .and. iord.eq.3) .or.
     +   (itype.eq.3 .and. iord.eq.2)) then
C--     Start scale must be at least in the second region
        if(ibin0.eq.1) stop 'sqcNStart: NNLO ibin0 .eq. 1 ---> STOP'
C--     Where to get the startvalue from
        iz1 = iz2sub6(ibin0-1) !Thats why ibin0 must be larger than 1
C--     Add discontinuity and store in iz0
        call sqcNSjup
     +          (itype,idout,ids,idg,iyg,iord,dlam,nyy2(iyg),iz1,iz0)
      endif
C--   Copy starting value to bin 0 (thats where sqcGridns looks for it)
      call sqcPCopjj(igout,iz0,igout,0)

      return
      end

C     ===============================================================
      subroutine sqcGridns(itype,
     +          ipdf,ids,idg,ityp,iyg,iord,dlam,it0,it1,it2,eps,ierr)
C     ===============================================================

C--   Steering routine for nonsinglet evolution on an equidistant y-subgrid
C--   This routine handles crossing of the flavor thresholds
C--   This routine does not set the starting values
C--
C--   1. it1 must be at lower boundary or at one of the thresholds
C--   2. it2 must be at upper boundary
C--   3. it0 must be >= it1 and <= it2
C--
C--   itype     (in) 1=unpol, 2=pol, 3=timelike
C--   ipdf      (in) nonsinglet pdf index
C--   ids       (in) singlet index
C--   idg       (in) gluon index 
C--   ityp      (in) 1=NS+, 2=NS-, 3=NSV
C--   iyg       (in) subgrid index 1,...,nyg2
C--   iord      (in) 1=LO , 2=NLO, 3=NNLO
C--   dlam(4:6) (in) heavy quark weights in NNLO nonsinglet jumps
C--   it0       (in) starting point          
C--   it1       (in) lower limit of evolution
C--   it2       (in) upper limit of evolution
C--   eps      (out) max deviation quad - lin interpolation
C--   ierr     (out) 1 = no alphas available at it1
C--
C--   Remark: ids, idg and dlam are used to calculate the NNLO jumps
C--           If, on input, ids = idg then  dont calculate the heavy
C--           flavor contribution to the jumps
C--
C--   Remark: at this point everything is in terms of A values

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension dlam(4:6)
      dimension izu1(4),izu2(4),nflu(4),izd1(4),izd2(4),nfld(4)
      
C--   Global identifier
      ipdfg  = iqcIdPdfLtoG(itype,ipdf)

C--   Check if alphas available at lowest t-bin
      if(it1.lt.itmin6) then
        ierr = 1
        return
      endif
      ierr = 0
      eps  = 0.D0

C--   Get evolution limits (kinda joblist of how to proceed)
      call sqcEvLims(it0,it1,it2,
     +               izu1,izu2,nflu,nup,izd1,izd2,nfld,ndn,idl,idu)
     
C--   Get upper y-index
      iymax = iqcIyMaxG(iymac2,iyg)

C--   Upward evolutions
      do i = 1,nup
        if(i.eq.1) then
C--       Copy starting value from bin 0 to bin iz
          call sqcPCopjj(ipdfg,0,ipdfg,izu1(i))
        else
C--       Copy starting value from previous evolution and add discontinuity
          dlm = dlam(nflu(i))
          call sqcNSjup(itype,
     +         ipdf,ids,idg,iyg,iord,dlm,iymax,izu2(i-1),izu1(i))
        endif
C--     Evolve upward with fixed nf
        call sqcNSevnf(itype,
     +                 ipdf,ityp,iyg,iord,nflu(i),izu1(i),izu2(i))
      enddo

C--   Downward evolutions (always w/linear interpolation)
      do i = 1,ndn
      
        if(i.eq.1) then
C--       Copy starting value from bin 0 to bin iz
          call sqcPCopjj(ipdfg,0,ipdfg,izd1(i))
        else
C--       Copy starting value from previous evolution and add discontinuity
          dlm = dlam(nfld(i))
          call sqcNSjdn(itype,
     +         ipdf,ids,idg,iyg,iord,dlm,iymax,izd2(i-1),izd1(i))
        endif
        
        if(ioy2.eq.2) then
        
C--       Evolve downward with current ioy2 = lin
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
     
        elseif(ioy2.eq.3 .and. niter6.lt.0) then
        
C--       Evolve downward with current ioy2 = quad
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
          
        elseif(ioy2.eq.3 .and. niter6.eq.0) then
        
C--       Convert starting A values from quad to linear and set ioy2 = 2
          call sqcGiQtoL(ipdfg,ipdfg,iymax,izd1(i),izd1(i))
C--       Evolve downward with linear interpolation
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
C--       Convert all A values from lin to quad and set ioy2 = 3
          call sqcGiLtoQ(ipdfg,ipdfg,iymax,izd2(i),izd1(i))

        elseif(ioy2.eq.3 .and. niter6.gt.0) then
        
C--       Switch to linear interpolation
          call sqcGiQtoL(ipdfg,ipdfg,iymax,izd1(i),izd1(i))
C--       Evolve downward with lin interpolation
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
C--       Remember startvalue
          call sqcNSStoreStart(itype,ipdf,1,iymax,izd1(i))
C--       Switch to quad interpolation
          call sqcGiLtoQ(ipdfg,ipdfg,iymax,izd2(i),izd2(i))
C--       Evolve upward with quad interpolation
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd2(i),izd1(i))
C--       Iteration loop
          do iter = 0,niter6
C--         Switch to linear interpolation
            call sqcGiQtoL(ipdfg,ipdfg,iymax,izd1(i),izd1(i))
C--         New starting value
            call sqcNSNewStart(itype,ipdf,1,iymax,izd1(i),eps)
C--         Finished?
            if(iter.eq.niter6) then
C--           Restore startvalue
              call sqcNSRestoreStart(itype,ipdf,1,iymax,izd1(i))
C--           Switch to quad interpolation
              call sqcGiLtoQ(ipdfg,ipdfg,iymax,izd1(i),izd1(i))
            else
C--           Evolve downward with lin interpolation
              call sqcNSevnf(itype,
     +                   ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
C--           Switch to quad interpolation
              call sqcGiLtoQ(ipdfg,ipdfg,iymax,izd2(i),izd2(i))
C--           Evolve upward with quad interpolation
              call sqcNSevnf(itype,
     +                   ipdf,ityp,iyg,iord,nfld(i),izd2(i),izd1(i))
            endif 
          enddo
        endif
        
C--   End of loop over downward evolutions        
      enddo
      
      eps0 = dqcGetEps(itype,ipdf,iymax,it0)
      eps1 = dqcGetEps(itype,ipdf,iymax,it1)
      eps2 = dqcGetEps(itype,ipdf,iymax,it2)
      eps  = max(eps0,eps1,eps2)

      return
      end

C     =========================================================
      subroutine sqcNSevnf(itype,ipdf,ityp,iyg,iord,nf,iz1,iz2)
C     =========================================================

C--   Non-singlet evolution from iz1 to iz2 at fixed nf
C--
C--   Input:  itype = 1=unpol, 2=pol, 3=timelike
C--           ipdf  = index of pdf to be evolved
C--           ityp  = 1=NS+, 2=NS-, 3=NSV
C--           iyg   = subgrid index 1,...,nyg2
C--           iord  = 1=LO,  2=NLO, 3=NNLO
C--           nf    = number of flavours
C--           iz1   = start of evolution
C--           iz2   = end of evolution
C--
C--   Output:  Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   NB: the routines sqcNSmult and sqcNSeqs can be found in qcdutil.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension ialf(mord0)

      dimension idwt(3)
C--               NS+ NS- NSV
      data idwt /  5,  6,  7  /

      dimension sbar(mxx0),ssum(mxx0),vmat(mxx0),bvec(mxx0),hvec(mxx0)
      
*mb
*mb      logical first,lprint
*mb      save first,lprint,nflast
*mb      data first /.true./
*mb      
*mb      if(first) then
*mb        first  = .false.
*mb        lprint = .true.
*mb        nflast = nf 
*mb      endif
*mb      
*mb      if(nf.ne.nflast .and. lprint) then
*mb        write(6,*) 'EVOL: change number of flavours'
*mb        lprint = .false.
*mb      endif
*mb                                
      
C--   Initialization
C--   --------------
C--   Base address alfas tables
      iset = 0                                 !take alfas from base set
      do i = 1,mord0
        ialf(i) = iqcIaAtab(1,i,iset)-1
      enddo
C--   Interpolation index
      idk     = ioy2-1
C--   Direction of evolution (isign) and first point after it1 (next)
      isign = 1
      next  = iz1+1
      if(iz2.lt.iz1) then 
        isign = -1
        next  = iz1-1
      endif
C--   Setup sbar
      do i = 1,nyy2(iyg)
        sbar(i) = 0.D0
        ssum(i) = 0.D0
      enddo
C--   Index limits of region iyg on subgrid iyg
      iy1 = 1
      iy2 = iqcIyMaxG(iymac2,iyg)

C--   Calculate vector b at input scale iz1
C--   -------------------------------------
C--   Grid spacing delta = z(next)-z (divided by 2)
      delt = 0.5*abs(zgrid2(next)-zgrid2(iz1))
C--   Transformation matrix divided by delta (Sbar)        
      do i = 1,nmaty2(ioy2)
        sbar(i) = smaty2(i,ioy2)/delt
      enddo
C--   Find t-index
      it1 = itfiz5(iz1)
C--   Weight matrix and  V = Sbar+W at t1
      do iy = 1,iy2
        vmat(iy) = sbar(iy)
      enddo
      do k = 1,iord
        id = idPij7(idwt(ityp),k,itype)
        as = stor7(ialf(k)+iz1)
        ia = iqcGaddr(stor7,1,it1,nf,iyg,id)-1
        do iy = 1,iy2
          ia = ia+1
          vmat(iy) = vmat(iy) + isign*as*stor7(ia)
        enddo  
      enddo      
C--   Address of pdf(iy=1,iz1,ipdf) in the store
      iadr = iqcPdfIjkl(1,iz1,ipdf,itype)
C--   Calculate V a = b
      call sqcNSmult(vmat,iy2,stor7(iadr),bvec,iy2)

C--   Evolution loop over t
C--   ---------------------
      do iz = next,iz2,isign
C--     Find t-index
        it = itfiz5(iz)
C--     Weight matrix and V matrix at t
        do iy = 1,iy2
          vmat(iy) = sbar(iy)
        enddo
        do k = 1,iord
          id = idPij7(idwt(ityp),k,itype)
          as = stor7(ialf(k)+iz)
          ia = iqcGaddr(stor7,1,it,nf,iyg,id)-1
          do iy = 1,iy2
            ia = ia+1
            vmat(iy) = vmat(iy) - isign*as*stor7(ia)
          enddo  
        enddo
C--     Address of pdf(1,iz,ipdf) in the store
        iadr = iqcPdfIjkl(1,iz,ipdf,itype)
C--     Solve V a = b
        call sqcNSeqs(vmat,iy2,stor7(iadr),bvec,iy2)
C--     Update b for the next iteration: not at last iteration of the loop
        if(iz.ne.iz2) then
C--       Grid spacing delta = z(next)-z (divided by 2)
          delt = 0.5*abs(zgrid2(iz+isign)-zgrid2(iz))
C--       Sum of current and next sbar; store next sbar       
          do i = 1,nmaty2(ioy2)
            sbnext  = smaty2(i,ioy2)/delt
            ssum(i) = sbar(i)+sbnext
            sbar(i) = sbnext
          enddo
C--       Calculate Ssum a = h
          call sqcNSmult(ssum,nmaty2(ioy2),stor7(iadr),hvec,iy2)
C--       Update b
          do iy = 1,iy2
            bvec(iy) = hvec(iy)-bvec(iy)
          enddo
        endif
      enddo

      return
      end

C     ===============================================================
      subroutine 
     +        sqcNSjup(itype,idq,ids,idg,iyg,iord,dlam,ny,izin,izout)
C     ===============================================================

C--   Calculate jump in nonsinglet for upward evolution
C--
C--   Input:  itype  = 1=unpol, 2=pol, 3=timelike
C--           idq    = index of nonsinglet pdf
C--           ids    = index of singlet pdf
C--           idg    = index of gluon pdf
C--           iyg    = subgrid index 1,...,nyg2         
C--           iord   = 1=LO,  2=NLO, 3=NNLO
C--           dlam   = weight of heavy quark discontinuity
C--           ny     = upper index limit yloop
C--           izin   = z-grid index where A(nf) is to be found 
C--           izout  = z-grid index where A(nf+1) will be stored
C--
C--   Output: Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   Remark: (1) Handle itype=1 iord=3 or itype=3 iord=2 otherwise copy
C--           (2) When ids = idg the heavy flavor contribution is set to zero

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      data idqq, idhq, idhg / 3, 4, 5 /

      dimension wmatqq(mxx0), wmathq(mxx0), wmathg(mxx0)
      dimension qjump(mxx0) , sjump(mxx0) , gjump(mxx0)

C--   Global identifiers
      igq  = iqcIdPdfLtoG(itype,idq)

C--   Alphas to use
      iset = 0
      as   = stor7(iqcIaAtab(izout,0,iset))
      assq = as*as

C--   Find t-index
      itin = itfiz5(izin)

C--   For spacelike evolution at NNLO
      if(itype.eq.1 .and. iord.eq.3) then
C--     Calculate quark weight Wqq
        nf = 3  !dummy variable since the A-weights do not depend on nf
        ia = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(2,2,iord,itype))-1
        do iy = 1,ny
          wmatqq(iy) = assq*stor7(ia+iy)
        enddo
C--     Add transformation matrix to Wqq
        do i = 1,nmaty2(ioy2)
          wmatqq(i) = wmatqq(i)+smaty2(i,ioy2)
        enddo
C--     Calculate qjump = W*a and store in buffer
        iaq   = iqcPdfIjkl(1,izin,idq,itype) !address of input nonsinglet
        call sqcNSmult(wmatqq,ny,stor7(iaq),qjump,ny)
C--     Add heavy flavor but only if ids .ne. idg
        if(ids.ne.idg) then
C--       Weight matrix
          iaq = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(3,2,iord,itype))-1
          iag = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(3,1,iord,itype))-1
          do iy = 1,ny
            wmathq(iy) = assq*stor7(iaq+iy)
            wmathg(iy) = assq*stor7(iag+iy)
          enddo
          ias = iqcPdfIjkl(1,izin,ids,itype) !address of input singlet
          iag = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--       Calculate jump = W*a and store in buffers sjump and gjump
          call sqcNSmult(wmathq,ny,stor7(ias),sjump,ny)
          call sqcNSmult(wmathg,ny,stor7(iag),gjump,ny)
C--       Add weighted heavy quark jump
          do i = 1,ny
            qjump(i) = qjump(i) + dlam * (sjump(i)+gjump(i))
          enddo
        endif
C--     Calculate anew by solving S*anew = qnew
        iaq  = iqcPdfIjkl(1,izout,idq,itype) !address of output nonsinglet
        call sqcNSeqs(smaty2(1,ioy2),nmaty2(ioy2),stor7(iaq),qjump,ny)

C--   For timelike evolution at NLO
      elseif(itype.eq.3 .and. iord.eq.2 .and. itlmc6.ne.0) then
C--     Add transformation matrix to Wqq = 0
        do i = 1,nmaty2(ioy2)
          wmatqq(i) = smaty2(i,ioy2)
        enddo
C--     Calculate qjump = W*a and store in buffer
        iaq   = iqcPdfIjkl(1,izin,idq,itype) !address of input nonsinglet
        call sqcNSmult(wmatqq,nmaty2(ioy2),stor7(iaq),qjump,ny)
C--     Add heavy flavor but only if ids .ne. idg
        if(ids.ne.idg) then
C--       Weight matrix
          iag = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(3,1,iord,itype))-1
          do iy = 1,ny
            wmathg(iy) = as*stor7(iag+iy)
          enddo
          iag = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--       Calculate jump = W*a and store in buffer gjump
          call sqcNSmult(wmathg,ny,stor7(iag),gjump,ny)
C--       Add weighted heavy quark jump
          do i = 1,ny
            qjump(i) = qjump(i) + dlam * gjump(i)
          enddo
        endif
C--     Calculate anew by solving S*anew = qnew
        iaq  = iqcPdfIjkl(1,izout,idq,itype) !address of output nonsinglet
        call sqcNSeqs(smaty2(1,ioy2),nmaty2(ioy2),stor7(iaq),qjump,ny)

C--   No jumps, just copy izin to izout
      else
        call sqcPCopjj(igq,izin,igq,izout)
      endif

      return
      end

C     ===============================================================
      subroutine 
     +        sqcNSjdn(itype,idq,ids,idg,iyg,iord,dlam,ny,izin,izout)
C     ===============================================================

C--   Calculate NNLO jump in nonsinglet for downward evolution
C--
C--   Input:  itype  = 1=unpol, 2=pol, 3=timelike
C--           idq    = index of nonsinglet pdf
C--           ids    = index of singlet pdf
C--           idg    = index of gluon pdf
C--           iyg    = subgrid index 1,...,nyg2         
C--           iord   = 1=LO,  2=NLO, 3=NNLO
C--           dlam   = weight of heavy quark discontinuity
C--           ny     = upper index limit yloop
C--           izin   = z-grid index where A(nf+1) is to be found 
C--           izout  = z-grid index where A(nf) will be stored
C--
C--   Output: Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   Remark: (1) Handle itype=1 iord=3 or itype=3 iord=2 otherwise copy
C--           (2) When ids = idg the heavy flavor contribution is set to zero

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      data idqq, idhq, idhg / 3, 4, 5 /

      dimension wmatqq(mxx0), wmathq(mxx0), wmathg(mxx0)
      dimension qstor(mxx0) , sstor(mxx0) , gstor(mxx0)

C--   Global identifiers
      igq  = iqcIdPdfLtoG(itype,idq)

C--   Alphas to use
      iset = 0
      as   = stor7(iqcIaAtab(izin,0,iset))
      assq = as*as

C--   Find t-index
      itin = itfiz5(izin)

C--   For spacelike evolution at NNLO
      if(itype.eq.1 .and. iord.eq.3) then
C--     Calculate quark weight Wqq
        nf = 3  !dummy variable since the A-weights do not depend on nf
        ia = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(2,2,iord,itype))-1
        do iy = 1,ny
          wmatqq(iy) = assq*stor7(ia+iy)
        enddo
C--     Add transformation matrix to Wqq
        do i = 1,nmaty2(ioy2)
          wmatqq(i) = wmatqq(i)+smaty2(i,ioy2)
        enddo
C--     Calculate qold = S*aold and store in buffer
        iaq   = iqcPdfIjkl(1,izin,idq,itype) !address of input nonsinglet
        call sqcNSmult(smaty2(1,ioy2),nmaty2(ioy2),stor7(iaq),qstor,ny)
C--     Subtract heavy flavor but only if ids .ne. idg
        if(ids.ne.idg) then
C--       Weight matrix
          iaq = iqcGaddr(stor7,iy,itin,nf,iyg,idAijk7(3,2,iord,itype))-1
          iag = iqcGaddr(stor7,iy,itin,nf,iyg,idAijk7(3,1,iord,itype))-1
          do iy = 1,ny
            wmathq(iy) = assq*stor7(iaq+iy)
            wmathg(iy) = assq*stor7(iag+iy)
          enddo
          ias = iqcPdfIjkl(1,izin,ids,itype) !address of input singlet
          iag = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--       Calculate jump = W*a and store in buffers sstor and gstor
          call sqcNSmult(wmathq,ny,stor7(ias),sstor,ny)
          call sqcNSmult(wmathg,ny,stor7(iag),gstor,ny)
C--       Subtract weighted heavy quark jump
          do i = 1,ny
            qstor(i) = qstor(i) - dlam * (sstor(i)+gstor(i))
          enddo
        endif
C--     Calculate anew by solving Wqq*anew = qstor
        iaq  = iqcPdfIjkl(1,izout,idq,itype) !address of output nonsinglet
        call sqcNSeqs(wmatqq,ny,stor7(iaq),qstor,ny)

C--   For timelike evolution at NLO
      elseif(itype.eq.3 .and. iord.eq.2 .and. itlmc6.ne.0) then
C--     Add transformation matrix to Wqq = 0
        do i = 1,nmaty2(ioy2)
          wmatqq(i) = smaty2(i,ioy2)
        enddo
C--     Calculate qold = S*aold and store in buffer
        iaq   = iqcPdfIjkl(1,izin,idq,itype) !address of input nonsinglet
        call sqcNSmult(smaty2(1,ioy2),nmaty2(ioy2),stor7(iaq),qstor,ny)
C--     Subtract heavy flavor but only if ids .ne. idg
        if(ids.ne.idg) then
C--       Weight matrix
          iag = iqcGaddr(stor7,iy,itin,nf,iyg,idAijk7(3,1,iord,itype))-1
          do iy = 1,ny
            wmathg(iy) = as*stor7(iag+iy)
          enddo
          iag = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--       Calculate jump = W*a and store in buffer gstor
          call sqcNSmult(wmathg,ny,stor7(iag),gstor,ny)
C--       Subtract weighted heavy quark jump
          do i = 1,ny
            qstor(i) = qstor(i) - dlam * gstor(i)
          enddo
        endif
C--     Calculate anew by solving Wqq*anew = qstor
        iaq  = iqcPdfIjkl(1,izout,idq,itype) !address of output nonsinglet
        call sqcNSeqs(wmatqq,nmaty2(ioy2),stor7(iaq),qstor,ny)

C--   No jumps, just copy izin to izout
      else
        call sqcPCopjj(igq,izin,igq,izout)
      endif

      return
      end

C     ==================================================
      subroutine sqcNSStoreStart(itype,ipdf,iy1,iy2,iz0)
C     ==================================================

C--   Store  startvalue (A values) into common block /stbuf/
C--
C--   itype  (in) 1=unpol, 2=pol, 3=timelike
C--   ipdf   (in) Nonsinglet table
C--   iy1    (in) First yloop index
C--   iy2    (in) Last yloop index
C--   iz0    (in) t-bin containing startvalue

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base address
      ia  = iqcPdfIjkl(iy1,iz0,ipdf,itype)-1
C--   Store startvalue
      do j = iy1,iy2
        ia = ia+1
        qtarg(j) = stor7(ia)
        qlast(j) = stor7(ia)
      enddo

      return
      end

C     ====================================================
      subroutine sqcNSNewStart(itype,ipdf,iy1,iy2,iz0,eps)
C     ====================================================

C--   Set new startvalue (A values) using common block /stbuf/
C--
C--   itype  (in) 1=unpol, 2=pol, 3=timelike
C--   ipdf   (in) Nonsinglet table
C--   iy1    (in) First yloop index
C--   iy2    (in) Last yloop index
C--   iz0    (in) z-bin containing startvalue
C--   eps    (out) Max |new-original|  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base address
      ia = iqcPdfIjkl(iy1,iz0,ipdf,itype)-1
C--   New startvalue
      eps = -999.D0
      do j = iy1,iy2
        ia = ia+1
        dif = stor7(ia)-qtarg(j)
        eps = max(eps, abs(dif))
        stor7(ia) = qlast(j)-dif
        qlast(j)  = stor7(ia)
      enddo

      return
      end

C     ====================================================
      subroutine sqcNSRestoreStart(itype,ipdf,iy1,iy2,iz0)
C     ====================================================

C--   Restore  original startvalue (A values) from common block /stbuf/
C--
C--   itype  1=unpol, 2=pol, 3=timelike
C--   ipdf   Nonsinglet table
C--   iy1    First yloop index
C--   iy2    Last yloop index
C--   iz0    z-bin containing startvalue  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base address
      ia = iqcPdfIjkl(iy1,iz0,ipdf,itype)-1
C--   Restore startvalue
      do j = iy1,iy2
        ia = ia+1
        stor7(ia) = qtarg(j)
      enddo

      return
      end

C     ================================================================
C     Singlet/Gluon evolution routines
C     ================================================================

C     ===========================================================
      subroutine sqcGridsg(itype,
     +                     idf,idg,iyg,iord,it0,it1,it2,eps,ierr)
C     ===========================================================

C--   Steering routine for singlet/gluon evolution on an equidistant subgrid
C--   This routine handles crossing of the flavor thresholds
C--   This routine does not set the starting values
C--
C--   1. it1 must be at lower boundary or at one of the thresholds
C--   2. it2 must be at upper boundary
C--   3. it0 must be >= it1 and <= it2
C--
C--   itype  (in) 1=unpol, 2=pol, 3=timelike
C--   idf    (in) pdf index for singlet
C--   idg    (in) pdf index for gluon
C--   iyg    (in) subgrid index 1,...,nyg2
C--   iord   (in) 1-LO , 2=NLO, 3=NNLO
C--   it0    (in) starting point           ) 
C--   it1    (in) lower limit of evolution ) it1 <= it0 <= it2
C--   it2    (in) upper limit of evolution )
C--   eps   (out) max deviation quad - lin interpolation
C--   ierr  (out) 1 = no alphas available at it1

C--   Remark: at this point everything is in terms of A values

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension izu1(4),izu2(4),nflu(4),izd1(4),izd2(4),nfld(4)

C--   Global identifiers
      idfg  = iqcIdPdfLtoG(itype,idf)
      idgg  = iqcIdPdfLtoG(itype,idg)
      
C--   Check if alphas available at lowest t-bin
      if(it1.lt.itmin6) then
        ierr = 1
        return
      endif
      ierr = 0
      eps  = 0.D0

C--   Get evolution limits (kinda joblist of how to proceed)
      call sqcEvLims(it0,it1,it2,
     +               izu1,izu2,nflu,nup,izd1,izd2,nfld,ndn,ibl,ibu)
     
C--   Get upper y-index
      iymax = iqcIyMaxG(iymac2,iyg)

C--   Upward evolutions
      do i = 1,nup
        if(i.eq.1) then
C--       Copy starting value from bin 0 to bin iz
          call sqcPCopjj(idfg,0,idfg,izu1(i))
          call sqcPCopjj(idgg,0,idgg,izu1(i))
        else
C--       Copy starting value from previous evolution and add discontinuity
          call sqcSGjup
     +           (itype,idf,idg,iyg,iord,iymax,izu2(i-1),izu1(i))
        endif
C--     Evolve upward with fixed nf
        call sqcSGevnf(itype,idf,idg,iyg,iord,nflu(i),izu1(i),izu2(i))
      enddo

C--   Downward evolutions (always w/linear interpolation) 
      do i = 1,ndn
      
        if(i.eq.1) then
C--       Copy starting value from bin 0 to bin iz
          call sqcPCopjj(idfg,0,idfg,izd1(i))
          call sqcPCopjj(idgg,0,idgg,izd1(i))
        else
C--       Copy starting value from previous evolution and add discontinuity
          call sqcSGjdn
     +           (itype,idf,idg,iyg,iord,iymax,izd2(i-1),izd1(i))
        endif
        
        if(ioy2.eq.2) then
        
C--       Evolve downward with current ioy2 = lin
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i)) 
     
        elseif(ioy2.eq.3 .and. niter6.lt.0) then   
                 
C--       Evolve downward with current ioy2 = quad 
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i))
     
        elseif(ioy2.eq.3 .and. niter6.eq.0) then
        
C--       Convert starting A values from quad to linear and set ioy2 = 2
          call sqcGiQtoL(idfg,idfg,iymax,izd1(i),izd1(i))
          call sqcGiQtoL(idgg,idgg,iymax,izd1(i),izd1(i))
C--       Evolve downward with linear interpolation
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i))
C--       Convert all A values from lin to quad and set ioy2 = 3
          call sqcGiLtoQ(idfg,idfg,iymax,izd2(i),izd1(i))
          call sqcGiLtoQ(idgg,idgg,iymax,izd2(i),izd1(i))

       elseif(ioy2.eq.3 .and. niter6.gt.0) then
       
C--       Switch to linear interpolation       
          call sqcGiQtoL(idfg,idfg,iymax,izd1(i),izd1(i))
          call sqcGiQtoL(idgg,idgg,iymax,izd1(i),izd1(i))
C--       Evolve downward with lin interpolation          
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i))
C--       Remember startvalue     
          call sqcSGStoreStart(itype,idf,idg,1,iymax,izd1(i))
C--       Switch to quad interpolation
          call sqcGiLtoQ(idfg,idfg,iymax,izd2(i),izd2(i))
          call sqcGiLtoQ(idgg,idgg,iymax,izd2(i),izd2(i))
C--       Evolve upward with quad interpolation     
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd2(i),izd1(i))
C--       Iteration loop
          do iter = 0,niter6     
C--         Switch to linear interpolation     
            call sqcGiQtoL(idfg,idfg,iymax,izd1(i),izd1(i))
            call sqcGiQtoL(idgg,idgg,iymax,izd1(i),izd1(i))
C--         New starting value          
            call sqcSGNewStart
     +              (itype,idf,idg,1,iymax,izd1(i),eps)
C--         Finished?
            if(iter.eq.niter6) then
C--           Restore startvalue          
              call sqcSGRestoreStart
     +                      (itype,idf,idg,1,iymax,izd1(i))
C--           Switch to quad interpolation              
              call sqcGiLtoQ(idfg,idfg,iymax,izd1(i),izd1(i))
              call sqcGiLtoQ(idgg,idgg,iymax,izd1(i),izd1(i))
            else
C--           Evolve downward with lin interpolation   
              call sqcSGevnf(itype,
     +                   idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i))
C--           Switch to quad interpolation     
              call sqcGiLtoQ(idfg,idfg,iymax,izd2(i),izd2(i))
              call sqcGiLtoQ(idgg,idgg,iymax,izd2(i),izd2(i))
C--           Evolve upward with quad interpolation            
              call sqcSGevnf(itype,
     +                   idf,idg,iyg,iord,nfld(i),izd2(i),izd1(i))
            endif 
          enddo
        endif

C--   End of loop over downward evolutions        
      enddo
      
      epf0 = dqcGetEps(itype,idf,iymax,it0)
      epf1 = dqcGetEps(itype,idf,iymax,it1)
      epf2 = dqcGetEps(itype,idf,iymax,it2)
      epg0 = dqcGetEps(itype,idg,iymax,it0)
      epg1 = dqcGetEps(itype,idg,iymax,it1)
      epg2 = dqcGetEps(itype,idg,iymax,it2)
      eps  = max(epf0,epf1,epf2,epg0,epg1,epg2)

      return
      end

C     =======================================================
      subroutine sqcSGevnf(itype,idf,idg,iyg,iord,nf,iz1,iz2)
C     =======================================================

C--   Singlet-gluon evolution from iz1 to iz2 at fixed nf
C--
C--   Input:  itype = 1=unpol, 2=pol, 3=timelike
C--           idf   = index of singlet pdf to be evolved
C--           idg   = index of gluon pdf
C--           iyg   = subgrid index 1,...,nyg2
C--           iord  = 1=LO,  2=NLO, 3=NNLO
C--           nf    = number of flavours
C--           iz1   = start of evolution
C--           iz2   = end of evolution
C--
C--   Output:  Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   NB: the routines sqcSGmult and sqcSGeqs can be found in qcdutil.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension ialf(mord0)

      dimension idwt(4)
C--                 QQ  QG  GQ  GG
      data idwt /    1,  2,  3,  4  /

      dimension sbar(mxx0,4),ssum(mxx0,4),vmat(mxx0,4)
      dimension bf(mxx0),bg(mxx0),hf(mxx0),hg(mxx0)

C--   Initialization
C--   --------------
C--   Base address alfas tables
      iset = 0
      do i = 1,mord0
        ialf(i) = iqcIaAtab(1,i,iset)-1
      enddo
C--   Interpolation index
      idk     = ioy2-1
C--   Direction of evolution (isign) and first point after it1 (next)
      isign = 1
      next  = iz1+1
      if(iz2.lt.iz1) then 
        isign = -1
        next  = iz1-1
      endif
C--   Setup sbar
      do j = 1,4
        do i = 1,nyy2(iyg)
          sbar(i,j) = 0.D0
          ssum(i,j) = 0.D0
        enddo
      enddo
      
C--   Set here yloop index range      
      iy1 = 1
      iy2 = iqcIyMaxG(iymac2,iyg)

C--   Calculate vector b at input scale iz1
C--   -------------------------------------
C--   Grid spacing delta = z(next)-z (divided by 2)
      delt = 0.5*abs(zgrid2(next)-zgrid2(iz1))
C--   Transformation matrix divided by delta (Sbar)        
      do i = 1,nmaty2(ioy2)
        sbar(i,1) = smaty2(i,ioy2)/delt
        sbar(i,4) = smaty2(i,ioy2)/delt
      enddo
C--   Weight matrix and  V = Sbar+W at t1
      do ityp = 1,4
        do iy = 1,iy2
          vmat(iy,ityp) = sbar(iy,ityp)
        enddo
C--     Find t-index
        it1 = itfiz5(iz1)
        do k = 1,iord
          id = idPij7(idwt(ityp),k,itype)
          as = stor7(ialf(k)+iz1)
          ia = iqcGaddr(stor7,1,it1,nf,iyg,id)-1
          do iy = 1,iy2
            ia = ia+1
            vmat(iy,ityp) = vmat(iy,ityp) + isign*as*stor7(ia)
          enddo
        enddo  
      enddo
C--   Address of pdf(iy=1,iz1,ipdf) in the store
      iaf = iqcPdfIjkl(1,iz1,idf,itype)
      iag = iqcPdfIjkl(1,iz1,idg,itype)
C--   Calculate V a = b
      call sqcSGmult(vmat(1,1),vmat(1,2),vmat(1,3),vmat(1,4),iy2,
     +               stor7(iaf),stor7(iag),bf,bg,iy2)

C--   Evolution loop over t
C--   ---------------------
      do iz = next,iz2,isign      
C--     Find t-index
        it = itfiz5(iz)
C--     Weight matrix and V matrix at t
        do ityp = 1,4
          do iy = 1,iy2
            vmat(iy,ityp) = sbar(iy,ityp)
          enddo
          do k = 1,iord
            id = idPij7(idwt(ityp),k,itype)
            as = stor7(ialf(k)+iz)
            ia = iqcGaddr(stor7,1,it,nf,iyg,id)-1
            do iy = 1,iy2
              ia = ia+1
              vmat(iy,ityp) = vmat(iy,ityp) - isign*as*stor7(ia)
            enddo  
          enddo
        enddo
C--     Address of pdf(1,iz,ipdf) in the store
        iaf = iqcPdfIjkl(1,iz,idf,itype)
        iag = iqcPdfIjkl(1,iz,idg,itype)
C--     Solve V a = b
        call sqcSGeqs(vmat(1,1),vmat(1,2),vmat(1,3),vmat(1,4),
     +                stor7(iaf),stor7(iag),bf,bg,iy2)

C--     Update b for the next iteration: not at last iteration of the loop
        if(iz.ne.iz2) then
C--       Grid spacing delta = z(next)-z (divided by 2)
          delt = 0.5*abs(zgrid2(iz+isign)-zgrid2(iz))
C--       Sum of current and next sbar; store next sbar       
          do i = 1,nmaty2(ioy2)
            sbnext    = smaty2(i,ioy2)/delt
            ssum(i,1) = sbar(i,1)+sbnext
            ssum(i,4) = sbar(i,4)+sbnext
            sbar(i,1) = sbnext
            sbar(i,4) = sbnext
          enddo
C--       Calculate Ssum a = h
          call sqcSGmult(ssum(1,1),ssum(1,2),ssum(1,3),ssum(1,4),
     +                   nmaty2(ioy2),stor7(iaf),stor7(iag),hf,hg,iy2)
C--       Update b
          do iy = 1,iy2
            bf(iy) = hf(iy)-bf(iy)
            bg(iy) = hg(iy)-bg(iy)
          enddo
        endif
      enddo

C--   Thats it...

      return
      end

C     =========================================================
      subroutine sqcSGjup(itype,ids,idg,iyg,iord,ny,izin,izout)
C     =========================================================

C--   Calculate NNLO singlet and gluon jump for upward evolution 
C--
C--   Input:  itype  = 1=unpol, 2=pol, 3=timelike
C--           ids    = index of singlet pdf
C--           idg    = index of gluon pdf
C--           iyg    = subgrid index 1,...,nyg2         
C--           iord   = 1=LO,  2=NLO, 3=NNLO
C--           ny     = upper index yloop
C--           izin   = z-grid index where A(nf) is to be found 
C--           izout  = z-grid index where A(nf+1) will be stored
C--
C--   Output: Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   Remark: (1) Handle itype=1 iord=3 or itype=3 iord=2 otherwise copy
C--           (2) Acts as a do-nothing when ids = idg

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      data idgq, idgg, idqq, idhq, idhg / 1, 2, 3, 4, 5 /

      dimension wmatqq(mxx0), wmatqg(mxx0), wmatgq(mxx0), wmatgg(mxx0)
      dimension fjump(mxx0) , gjump(mxx0)

C--   Global identifiers
      igs  = iqcIdPdfLtoG(itype,ids)
      igg  = iqcIdPdfLtoG(itype,idg)

C--   Alphas to use
      iset = 0
      as   = stor7(iqcIaAtab(izout,0,iset))
      assq = as*as

C--   Find t-index
      itin = itfiz5(izin)

      nf = 3  !dummy variable since the A-weights do not depend on nf

C--   For spacelike evolution at NNLO
      if(itype.eq.1 .and. iord.eq.3 .and. ids.ne.idg) then
C--     Weight matrix
        iaqq = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(2,2,iord,itype))-1
        iahq = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(3,2,iord,itype))-1
        iahg = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(3,1,iord,itype))-1
        iagq = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(1,2,iord,itype))-1
        iagg = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(1,1,iord,itype))-1
        do iy = 1,ny
          wmatqq(iy) = assq * (stor7(iaqq+iy) + stor7(iahq+iy))
          wmatqg(iy) = assq *  stor7(iahg+iy)
          wmatgq(iy) = assq *  stor7(iagq+iy)
          wmatgg(iy) = assq *  stor7(iagg+iy)
        enddo
C--     Add transformation matrix to Wqq and Wgg
        do i = 1,nmaty2(ioy2)
          wmatqq(i) = wmatqq(i)+smaty2(i,ioy2)
          wmatgg(i) = wmatgg(i)+smaty2(i,ioy2)
        enddo
C--     Address of a(iy=1,izin,ipdf) in the store
        ias   = iqcPdfIjkl(1,izin,ids,itype) !address of input quark
        iag   = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--     Calculate fnew = W*a and store in buffers
        call sqcSGmult(wmatqq,wmatqg,wmatgq,wmatgg,ny,
     &                 stor7(ias),stor7(iag),fjump,gjump,ny)
C--     Calculate anew by solving S*anew = fnew
        ias   = iqcPdfIjkl(1,izout,ids,itype) !address of output singlet
        iag   = iqcPdfIjkl(1,izout,idg,itype) !address of output gluon
        call sqcNSeqs(smaty2(1,ioy2),nmaty2(ioy2),stor7(ias),fjump,ny)
        call sqcNSeqs(smaty2(1,ioy2),nmaty2(ioy2),stor7(iag),gjump,ny)

C--   For timelike evolution at NLO
      elseif(itype.eq.3 .and. iord.eq.2 .and. ids.ne.idg .and.
     &       itlmc6.ne.0) then
      call sqcPCopjj(igg,izin,igg,izout)
      call sqcNSjup(itype,ids,ids,idg,iyg,iord,1.D0,ny,izin,izout)

C--   No jumps, just copy izin to izout
      else
        call sqcPCopjj(igs,izin,igs,izout)
        call sqcPCopjj(igg,izin,igg,izout)
      endif

      return
      end

C     =========================================================
      subroutine sqcSGjdn(itype,ids,idg,iyg,iord,ny,izin,izout)
C     =========================================================

C--   Calculate NNLO singlet and gluon jump for downward evolution 
C--
C--   Input:  itype  = 1=unpol, 2=pol, 3=timelike
C--           ids    = index of singlet pdf
C--           idg    = index of gluon pdf
C--           iyg    = subgrid index 1,...,nyg2         
C--           iord   = 1=LO,  2=NLO, 3=NNLO
C--           ny     = upper index yloop
C--           izin   = z-grid index where A(nf+1) is to be found 
C--           izout  = z-grid index where A(nf) will be stored
C--
C--   Output: Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   Remark: (1) Handle itype=1 iord=3 or itype=3 iord=2 otherwise copy
C--           (2) Acts as a do-nothing when ids = idg

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      data idgq, idgg, idqq, idhq, idhg / 1, 2, 3, 4, 5 /

      dimension wmatqq(mxx0), wmatqg(mxx0), wmatgq(mxx0), wmatgg(mxx0)
      dimension fstor(mxx0) , gstor(mxx0)

C--   Global identifiers
      igs  = iqcIdPdfLtoG(itype,ids)
      igg  = iqcIdPdfLtoG(itype,idg)

C--   Alphas to use
      iset = 0
      as   = stor7(iqcIaAtab(izin,0,iset))
      assq = as*as
      
C--   Find t-index
      itin = itfiz5(izin)

      nf = 3  !dummy variable since the A-weights do not depend on nf

C--   For spacelike evolution at NNLO
      if(itype.eq.1 .and. iord.eq.3 .and. ids.ne.idg) then
C--     Weight matrix
        iaqq = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(2,2,iord,itype))-1
        iahq = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(3,2,iord,itype))-1
        iahg = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(3,1,iord,itype))-1
        iagq = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(1,2,iord,itype))-1
        iagg = iqcGaddr(stor7,1,itin,nf,iyg,idAijk7(1,1,iord,itype))-1
        do iy = 1,ny
          wmatqq(iy) = assq * (stor7(iaqq+iy) + stor7(iahq+iy))
          wmatqg(iy) = assq *  stor7(iahg+iy)
          wmatgq(iy) = assq *  stor7(iagq+iy)
          wmatgg(iy) = assq *  stor7(iagg+iy)
        enddo
C--     Add transformation matrix to Wqq and Wgg
        do i = 1,nmaty2(ioy2)
          wmatqq(i) = wmatqq(i)+smaty2(i,ioy2)
          wmatgg(i) = wmatgg(i)+smaty2(i,ioy2)
        enddo
C--     Address of a(iy=1,izin,ipdf) in the store
        ias   = iqcPdfIjkl(1,izin,ids,itype) !address of input quark
        iag   = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--     Calculate fold = S*a and store in buffers
        call sqcNSmult(smaty2(1,ioy2),nmaty2(ioy2),stor7(ias),fstor,ny)
        call sqcNSmult(smaty2(1,ioy2),nmaty2(ioy2),stor7(iag),gstor,ny)
C--     Calculate anew by solving W*anew = fold
        ias   = iqcPdfIjkl(1,izout,ids,itype) !address of output singlet
        iag   = iqcPdfIjkl(1,izout,idg,itype) !address of output gluon
        call sqcSGeqs(wmatqq,wmatqg,wmatgq,wmatgg,
     &                stor7(ias),stor7(iag),fstor,gstor,ny)

C--   For timelike evolution at NLO
      elseif(itype.eq.3 .and. iord.eq.2 .and. ids.ne.idg .and.
     &       itlmc6.ne.0) then
      call sqcPCopjj(igg,izin,igg,izout)
      call sqcNSjdn(itype,ids,ids,idg,iyg,iord,1.D0,ny,izin,izout)

C--   No jumps, just copy izin to izout
      else
        call sqcPCopjj(igs,izin,igs,izout)
        call sqcPCopjj(igg,izin,igg,izout)
      endif

      return
      end

C     =====================================================
      subroutine sqcSGStoreStart(itype,ids,idg,iy1,iy2,iz0)
C     =====================================================

C--   Store  startvalue (A values) into common block /stbuf/
C--
C--   itype  1=unpol, 2=pol, 3=timelike
C--   ids    Singlet table
C--   idg    Gluon table
C--   iy1    First yloop index
C--   iy2    Last  yloop index
C--   iz0    z-bin containing startvalue  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base addresses
      ias  = iqcPdfIjkl(iy1,iz0,ids,itype)-1
      iag  = iqcPdfIjkl(iy1,iz0,idg,itype)-1
C--   Store startvalues
      do j = iy1,iy2
        ias = ias+1
        iag = iag+1
        qtarg(j) = stor7(ias)
        gtarg(j) = stor7(iag)
        qlast(j) = stor7(ias)
        glast(j) = stor7(iag)
      enddo

      return
      end

C     =======================================================
      subroutine sqcSGNewStart(itype,ids,idg,iy1,iy2,iz0,eps)
C     =======================================================

C--   Set new startvalue (A values) using common block /stbuf/
C--
C--   itype  1=unpol, 2=pol, 3=timelike
C--   ids    Singlet table
C--   idg    Gluon table
C--   iy1    First yloop index
C--   iy2    Last yloop index
C--   iz0    z-bin containing startvalue
C--   eps    (out) Max |new-original|  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base addresses
      ias = iqcPdfIjkl(iy1,iz0,ids,itype)-1
      iag = iqcPdfIjkl(iy1,iz0,idg,itype)-1
C--   New startvalues
      eps = -999.D0
      do j = iy1,iy2
        ias  = ias+1
        iag  = iag+1
        difs = stor7(ias)-qtarg(j)
        difg = stor7(iag)-gtarg(j)
        eps  = max(eps, abs(difs))
        eps  = max(eps, abs(difg))
        stor7(ias) = qlast(j)-difs
        stor7(iag) = glast(j)-difg
        qlast(j)   = stor7(ias)
        glast(j)   = stor7(iag)
      enddo

      return
      end

C     =======================================================
      subroutine sqcSGRestoreStart(itype,ids,idg,iy1,iy2,iz0)
C     =======================================================

C--   Restore  original startvalue (A values) from common block /stbuf/
C--
C--   itype  1=unpol, 2=pol, 3=timelike
C--   ids    Singlet table
C--   idg    Gluon table
C--   iy1    First yloop index
C--   iy2    Last yloop index
C--   iz0    z-bin containing startvalue  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base addresses
      ias = iqcPdfIjkl(iy1,iz0,ids,itype)-1
      iag = iqcPdfIjkl(iy1,iz0,idg,itype)-1
C--   Restore startvalues
      do j = iy1,iy2
        ias = ias+1
        iag = iag+1
        stor7(ias) = qtarg(j)
        stor7(iag) = gtarg(j)
      enddo

      return
      end

C     ===================================================
      double precision function dqcGetEps(itype,id,ny,it)
C     ===================================================

C--   Get max deviation quad-lin at midpoints
C--   This works on a subgrid with A coefficients
C--
C--   itype   (in)  : 1=unpol, 2=pol, 3=timelike
C--   id      (in)  : Pdf identifier
C--   ny      (in)  : upper loop index in y
C--   it      (in)  : t-grid point

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      
      dimension epsi(mxx0)
      
      dqcGetEps = 0.D0
      if(ioy2.ne.3) return

C--   Base address
      iz = izfit5(it)
      ia = iqcPdfIjkl(1,iz,id,itype)
C--   Now get vector of deviations        
      call sqcDHalf(ioy2,stor7(ia),epsi,ny)
C--   Max deviation
      do iy = 1,ny
        dqcGetEps = max(dqcGetEps,abs(epsi(iy)))
      enddo
      
      return
      end

C     =============================================================
      subroutine sqcEvLims(it0,it1,it2,
     +               izu1,izu2,nflu,nup,izd1,izd2,nfld,ndn,ibl,ibu)
C     =============================================================

C--   Setup ranges for up and downward evolution 
C--
C--   it0         (in)   tgrid start point
C--   it1         (in)   lower tgrid limit of the evolution
C--   it2         (in)   upper tgrid limit of the evolution
C--   izu1(4)     (out)  lower zgrid limits of upward evolutions
C--   izu2(4)     (out)  upper zgrid limits of upward evolutions
C--   nflu(4)     (out)  number of flavors of each upward evolution
C--   nup         (out)  number of upward evolutions
C--   izd1(4)     (out)  upper zgrid limits of downward evolutions
C--   izd2(4)     (out)  lower zgrid limits of downward evolutions
C--   nfld(4)     (out)  number of flavors of each downward evolution
C--   ndn         (out)  number of downward evolutions
C--   ibl         (out)  lowest nf bin accessed
C--   ibu         (out)  highest nf bin accessed

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'

      dimension izu1(4),izu2(4),nflu(4),izd1(4),izd2(4),nfld(4)

C--   Initialize
      nup = 0
      ndn = 0
      idl = 0
      idu = 0
      do i = 1,4
        izu1(i) = 0
        izu2(i) = 0
        nflu(i) = 0
        izd1(i) = 0
        izd2(i) = 0
        nfld(i) = 0
      enddo
C--   Adjust upper and lower evolution limits to grid limits
      jt1    = max(it1,1)
      jt2    = min(it2,ntt2)
C--   Check it0 is in between grid limits
      if(it0.lt.jt1 .or. it0.gt.jt2) return !it0 out of range 

C--   Find out in which region it0, jt1 and jt2 falls
      ibin0 = 0
      ibin1 = 0
      ibin2 = 0
      do i = 1,ntsubg6
        if(it1sub6(i).le.it0 .and. it0.le.it2sub6(i)) ibin0 = i
        if(it1sub6(i).le.jt1 .and. jt1.le.it2sub6(i)) ibin1 = i
        if(it1sub6(i).le.jt2 .and. jt2.le.it2sub6(i)) ibin2 = i
      enddo
      if(ibin0.eq.0) return   !it0 out of range
      if(ibin1.eq.0) return   !it1 out of range
      if(ibin2.eq.0) return   !it2 out of range
C--   Calculate iz0; izmin-itmin gives the offset between iz and it
      iz0 = it0 + iz1sub6(ibin0)-it1sub6(ibin0)
      izd = it1 + iz1sub6(ibin1)-it1sub6(ibin1)
      izu = it2 + iz1sub6(ibin2)-it1sub6(ibin2)
C--   Upward from iz0, evolve flavor by flavor
      do i = ibin0,ibin2
        iz1 = max(iz0,iz1sub6(i))
        iz2 = min(iz2sub6(i),izu)
        nf  = nfsubg6(i)
        if(iz1.lt.iz2) then
          nup       = nup+1
          izu1(nup) = iz1
          izu2(nup) = iz2
          nflu(nup) = nf
        endif
      enddo
C--   Downward from iz0, evolve flavor by flavor
      do i = ibin0,ibin1,-1
        iz1 = min(iz0,iz2sub6(i))
        iz2 = max(iz1sub6(i),izd)
        nf  = nfsubg6(i)
        if(iz2.lt.iz1) then
          ndn       = ndn+1
          izd1(ndn) = iz1
          izd2(ndn) = iz2
          nfld(ndn) = nf
        endif
      enddo
C--   Range of regions accessed
      ibl = ibin1
      ibu = ibin2

      return
      end
      
      
      
      
