

C--   This is the file srcMatch.f with matching routines
C--
C--   subroutine sqcDoJumps(ityp,mset,nf0,iz0,idz,pdiff)
C--
C--   subroutine sqcMatchUnpIntrins(itype,mset,ig,iord,iz,ny)
C--   subroutine sqcMatchUnpDynamic(itype,mset,ig,iord,iz,ny)
C--   subroutine sqcMatchPol(itype,mset,ig,iord,iz,ny)
C--   subroutine sqcMatchTml(itype,mset,ig,iord,iz,ny)
C--
C--   subroutine sqcPMDelta(ww,iw,ic,idi,ido,fff,iord,ig,ny,iz)
C--   subroutine sqcJump11(ww,iw,cc,iai,iao,fff,iord,ig,ny,iz,iter)
C--   subroutine sqcFMatch11(vv,ff,fp,ny)
C--   subroutine sqcBMatch11(vv,ff,fp,ny,iter)
C--   subroutine sqcJump22(ww,iw,cc,iai,iao,fff,iord,ig,ny,iz,iter)
C--   subroutine sqcFMatch22(vv1,vv2,vv3,vv4,ff1,ff2,fp1,fp2,ny)
C--   subroutine sqcBMatch22(vv1,vv2,vv3,vv4,ff1,ff2,fp1,fp2,ny,iter)
C--   subroutine sqcJumpNN(ww,iw,cc,idim,iai,iao,fff,iord,ig,ny,iz,n)
C--   subroutine sqcGetGSH(ide,ig,ny,iz,iaf,fff,iepm)
C--   subroutine sqcPutGSH(ide,ig,ny,iz,iaf,fff,iepm)
C--   subroutine sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C     ==================================================
      subroutine sqcDoJumps(ityp,mset,nf0,iz0,idz,pdiff)
C     ==================================================

C--   Calculate jumps at the thresholds
C--
C--   ityp         (in) : unpol/pol/timelike [1,3]
C--   mset         (in) : pdf set
C--   nf0          (in) : nf at start of next evolution
C--   iz0          (in) : iz at start of next evolution
C--   idz          (in) : +1 forward, -1 backward evolution
C--   pdiff       (out) : difference array, set to zero here
C--
C--   The startvalues ff0 are stored in bin iz = -ig for all subgrids ig
C--   This routine updates ff0 by applying the jumps

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension pdiff(2,mxx0,*)

C--   Jump is not transferred via pdiff so set it to zero
      do id = 1,12
        do iy = 1,mxx0
          pdiff(1,iy,id) = 0.D0
          pdiff(2,iy,id) = 0.D0
        enddo
      enddo

C--   nf above threshold is nf0 for idz = +1 and nf0+1 for idz = -1
      nfabove = nf0 + (1-idz)/2
      nfbelow = nfabove-1
C--   iz above threshold is iz0 for idz = +1 and iz0+1 for idz = -1
      izabove = iz0 + (1-idz)/2
      izbelow = izabove-1
C--   iz1,2
      if(idz .eq. 1) then
        iz1 = izbelow
        iz2 = izabove
      else
        iz1 = izabove
        iz2 = izbelow
      endif

      if(idbug6.ge.1) then
        nf1 = itfiz5(-iz1)
        nf2 = itfiz5(-iz2)
        iq1 = itfiz5(iz1)
        iq2 = itfiz5(iz2)
        write(lunerr1,'('' JUMPNF iq1,2 = '',2I5,''   nf = '',2I3)')
     +        iq1,iq2,nf1,nf2
      endif

C--   Check
      nfa = itfiz5(-izabove)
      nfb = itfiz5(-izbelow)
      if(nfabove.ne.nfa) stop ' sqcDoJumps nfup problem'
      if(nfbelow.ne.nfb) stop ' sqcDoJumps nfdn problem'

      if(ityp.eq.1) then
C--     Unpolarised
        if(abs(nfix6).eq.0 .and. idz.eq.1) then
C--       Loop over subgrids (from coarse to dense)
          do ig = nyg2,1,-1
C--         Number of subgrid points after cuts
            ny = iqcIyMaxG(iymac2,ig)
            call sqcMatchUnpDynamic(ityp,mset,ig,iord6,iz1,ny)
          enddo
        else
C--       Loop over subgrids (from coarse to dense)
          nfixrem = nfix6
          nfix6   = 1
          do ig = nyg2,1,-1
C--         Number of subgrid points after cuts
            ny = iqcIyMaxG(iymac2,ig)
            call sqcMatchUnpIntrins(ityp,mset,ig,iord6,iz1,ny)
          enddo
          nfix6 = nfixrem
        endif

      elseif(ityp.eq.2) then
C--     Polarised
C--     Loop over subgrids (from coarse to dense)
        do ig = nyg2,1,-1
C--       Number of subgrid points after cuts
          ny = iqcIyMaxG(iymac2,ig)
          call sqcMatchPol(ityp,mset,ig,iord6,iz1,ny)
        enddo

      elseif(ityp.eq.3) then
C--     Time-like
C--     Loop over subgrids (from coarse to dense)
        do ig = nyg2,1,-1
C--       Number of subgrid points after cuts
          ny =  iqcIyMaxG(iymac2,ig)
          call sqcMatchTml(ityp,mset,ig,iord6,iz1,ny)
        enddo
      endif

      return
      end


C=======================================================================
C==   New code  ========================================================
C=======================================================================

C     =======================================================
      subroutine sqcMatchUnpIntrins(itype,mset,ig,iord,iz,ny)
C     =======================================================

C--   Forward and reverse matching of intrinsic unpolarised pdfs
C--
C--   itype (in): 1=unpol, 2=pol, 3=timelike
C--   mset  (in): pdf set
C--   ig    (in): subgrid index 1,...,nyg2
C--   iord  (in): 1=LO,  2=NLO, 3=NNLO
C--   iz    (in): iz point before matching
C--   ny    (in): upper index yloop

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension iw1(3), cc1(3), iw2(2,2,3), cc2(2,2,3), iai(2), iao(2)

      dimension fff(14*mxx0)
      dimension ide(13),iaf(14)

      Logical Lbelow

*mb      write(6,*) 'Match intrinsic'

C--   Only for unpolarised dynamic heavy flavours
      if(itype.ne.1 .or. abs(nfix6).ne.1)
     +                   stop 'sqcMatchUnpIntrins wrong type'

C--   Initialise
      do k = 1,3
        iw1(k) = 0
        cc1(k) = 0.D0
        do j = 1,2
          do i = 1,2
            iw2(i,j,k) = 0
            cc2(i,j,k) = 0.D0
          enddo
        enddo
      enddo

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C--   Get global ids
      do i = 0,12
        ide(i+1) = iqcIdPdfLtoG(mset,i)
      enddo
C--   Transform EPM from memory to GSH and put into fff
      call sqcGetGSH(ide,ig,ny,iz1,iaf,fff,0)

C--   Alphas to use
      iset = 0
      as   = stor7(iqcIaAtab(iza,0,iset))
      assq = as*as

C--   Identifiers and addresses
      idg   = 0
      ids   = 1
      idv   = 7
      idhp  = nfa
      idhm  = nfa+6
      iag   = idg*ny+1
      ias   = ids*ny+1
      iav   = idv*ny+1
      iahp  = idhp*ny+1
      iahm  = idhm*ny+1

C--   Forward GSH matching NNLO ----------------------------------------
      if(idz.eq.+1 .and. iord.le.3) then
C--     Get ~G' and ~H+'
        iw2(1,1,1) = idAijk7(1,1,1,itype)                            !A0
        cc2(1,1,1) = 1.D0
        iw2(1,1,3) = idAijk7(1,1,3,itype)                       !as2 Agg
        cc2(1,1,3) = assq
        iw2(1,2,2) = idAijk7(1,3,2,itype)                        !as Agh
        cc2(1,2,2) = as
        iw2(2,1,3) = idAijk7(3,1,3,itype)                       !as2 Ahg
        cc2(2,1,3) = assq
        iw2(2,2,1) = idAijk7(3,3,1,itype)                            !A0
        cc2(2,2,1) = 1.D0
        iw2(2,2,2) = idAijk7(3,3,2,itype)                        !as Ahh
        cc2(2,2,2) = as
        iai(1)     = iag
        iai(2)     = iahp
        iao(1)     = iag
        iao(2)     = iahp
        call sqcJump22(stor7,iw2,cc2,iai,iao,fff,iord,ig,ny,iz,0)
C--     Transform ~G' --> G'
        iw1(1) = 0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(1,2,3,itype)                           !as2 Agq
        cc1(3) = assq
        call sqcPMDelta(stor7,iw1,cc1,ids,idg,fff,iord,ig,ny,iz)
C--     Transform ~H+' --> H+'
        iw1(1) = 0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(3,2,3,itype)                           !as2 Ahq
        cc1(3) = assq
        call sqcPMDelta(stor7,iw1,cc1,ids,idhp,fff,iord,ig,ny,iz)
C--     Get S' and E2+',...,En+'
        iw1(1) = idAijk7(2,2,1,itype)                                !A0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(2,2,3,itype)                           !as2 Aqq
        cc1(3) = assq
        do i = 1,nfb
          ia = i*ny+1
          call sqcJump11(stor7,iw1,cc1,ia,ia,fff,iord,ig,ny,iz,0)
        enddo
C--     Get V' and E2-',...,En-'
        iw1(1) = idAijk7(2,2,1,itype)                                !A0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(2,2,3,itype)                           !as2 Aqq
        cc1(3) = assq
        do i = 7,nfb+6
          ia = i*ny+1
          call sqcJump11(stor7,iw1,cc1,ia,ia,fff,iord,ig,ny,iz,0)
        enddo
C--     Get H-'
        iw1(1) = idAijk7(2,2,1,itype)                                !A0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(3,3,2,itype)                            !as Ahh
        cc1(3) = as
        call sqcJump11(stor7,iw1,cc1,iahm,iahm,fff,iord,ig,ny,iz,0)

C--   Reverse GSH matching NNLO ----------------------------------------
      elseif(idz.eq.-1 .and. iord.le.3) then
*mb        write(6,*) 'Hurray reverse match'
C--     Get H-
        iw1(1) = idAijk7(2,2,1,itype)                                !A0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(3,3,2,itype)                            !as Ahh
        cc1(3) = as
        call sqcJump11(stor7,iw1,cc1,iahm,iahm,fff,iord,ig,ny,iz,10)
C--     Get V and E2-,...,En-
        iw1(1) = idAijk7(2,2,1,itype)                                !A0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(2,2,3,itype)                           !as2 Aqq
        cc1(3) = assq
        do i = 7,nfb+6
          ia = i*ny+1
          call sqcJump11(stor7,iw1,cc1,ia,ia,fff,iord,ig,ny,iz,0)
        enddo
C--     Get S and E2+,...,En+
        iw1(1) = idAijk7(2,2,1,itype)                                !A0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(2,2,3,itype)                           !as2 Aqq
        cc1(3) = assq
        do i = 1,nfb
          ia = i*ny+1
          call sqcJump11(stor7,iw1,cc1,ia,ia,fff,iord,ig,ny,iz,0)
        enddo
C--     Transform H+' --> ~H+'
        iw1(1) = 0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(3,2,3,itype)                          !-as2 Ahq
        cc1(3) = -assq
        call sqcPMDelta(stor7,iw1,cc1,ids,idhp,fff,iord,ig,ny,iz)
C--     Transform G' --> ~G'
        iw1(1) = 0
        cc1(1) = 1.D0
        iw1(3) = idAijk7(1,2,3,itype)                          !-as2 Agq
        cc1(3) = -assq
        call sqcPMDelta(stor7,iw1,cc1,ids,idg,fff,iord,ig,ny,iz)
C--     Get G' and H+'
        iw2(1,1,1) = idAijk7(1,1,1,itype)                            !A0
        cc2(1,1,1) = 1.D0
        iw2(1,1,3) = idAijk7(1,1,3,itype)                       !as2 Agg
        cc2(1,1,3) = assq
        iw2(1,2,2) = idAijk7(1,3,2,itype)                        !as Agh
        cc2(1,2,2) = as
        iw2(2,1,3) = idAijk7(3,1,3,itype)                       !as2 Ahg
        cc2(2,1,3) = assq
        iw2(2,2,1) = idAijk7(3,3,1,itype)                            !A0
        cc2(2,2,1) = 1.D0
        iw2(2,2,2) = idAijk7(3,3,2,itype)                        !as Ahh
        cc2(2,2,2) = as
        iai(1)     = iag
        iai(2)     = iahp
        iao(1)     = iag
        iao(2)     = iahp
        call sqcJump22(stor7,iw2,cc2,iai,iao,fff,iord,ig,ny,iz,10)

C--   Forward GSH matching NLO -----------------------------------------
      elseif(idz.eq.+1 .and. iord.eq.2) then
C--     Get G'
        iw1(1) = 0
        cc1(1) = 1.D0
        iw1(2) = idAijk7(1,3,2,itype)                            !as Agh
        cc1(2) = as
        call sqcPMDelta(stor7,iw1,cc1,idg,idhp,fff,iord,ig,ny,iz)
C--     Get H+'
        iw1(1) = idAijk7(3,3,2,itype)                                !A0
        cc1(1) = 1.D0
        iw1(2) = idAijk7(3,3,2,itype)                            !as Ahh
        cc1(2) = as
        call sqcJump11(stor7,iw1,cc1,iahp,iahp,fff,iord,ig,ny,iz,0)
C--     Get H-'
        call sqcJump11(stor7,iw1,cc1,iahm,iahm,fff,iord,ig,ny,iz,0)

C--   Reverse GSH matching NLO -----------------------------------------
      elseif(idz.eq.-1 .and. iord.eq.2) then
C--     Get H-
        iw1(1) = idAijk7(3,3,2,itype)                                !A0
        cc1(1) = 1.D0
        iw1(2) = idAijk7(3,3,2,itype)                            !as Ahh
        cc1(2) = as
        call sqcJump11(stor7,iw1,cc1,iahm,iahm,fff,iord,ig,ny,iz,0)
C--     Get H+
        call sqcJump11(stor7,iw1,cc1,iahp,iahp,fff,iord,ig,ny,iz,0)
C--     Get G'
        iw1(1) = 0
        cc1(1) = 1.D0
        iw1(2) = idAijk7(1,3,2,itype)                           !-as Agh
        cc1(2) = -as
        call sqcPMDelta(stor7,iw1,cc1,idg,idhp,fff,iord,ig,ny,iz)

      endif
C--   End of GSH matching ----------------------------------------------

C--   Transform GSH from fff to EPM and put into memory
      call sqcPutGSH(ide,ig,ny,iz2,iaf,fff,0)

      return
      end

C     =======================================================
      subroutine sqcMatchUnpDynamic(itype,mset,ig,iord,iz,ny)
C     =======================================================

C--   Forward and reverse matching of dynamic unpolarised pdfs
C--
C--   itype (in): 1=unpol, 2=pol, 3=timelike
C--   mset  (in): pdf set
C--   ig    (in): subgrid index 1,...,nyg2
C--   iord  (in): 1=LO,  2=NLO, 3=NNLO
C--   iz    (in): iz point before matching
C--   ny    (in): upper index yloop

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension iw(3), cc(3)

      dimension fff(14*mxx0)
      dimension ide(13),iaf(14)

      Logical Lbelow

*mb      write(6,*) 'Match dynamic'

C--   Only for unpolarised dynamic heavy flavours
      if(itype.ne.1 .or. abs(nfix6).ne.0)
     +                   stop 'sqcMatchUnpDynamic wrong type'

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C--   Get global ids
      do i = 0,12
        ide(i+1) = iqcIdPdfLtoG(mset,i)
      enddo
C--   Transform EPM from memory to GSH and put into fff
      call sqcGetGSH(ide,ig,ny,iz1,iaf,fff,0)

C--   Alphas to use
      iset = 0
      as   = stor7(iqcIaAtab(iza,0,iset))
      assq = as*as

C--   Identifiers and addresses
      idg   = 0
      ids   = 1
      idv   = 7
      idhp  = nfa
      idhm  = nfa+6
      iag   = idg*ny+1
      ias   = ids*ny+1
      iav   = idv*ny+1
      iahp  = idhp*ny+1
      iahm  = idhm*ny+1

C--   No NLO, thank you
      iw(2) = 0
      cc(2) = 1.D0

C--   Forward GSH matching NNLO ----------------------------------------
      if(idz.eq.+1 .and. iord.eq.3) then

C--     Get ~H+'
        iw(1) = 0
        cc(1) = 1.D0
        iw(3) = idAijk7(3,1,3,itype)                            !as2 Ahg
        cc(3) = assq
        call sqcJump11(stor7,iw,cc,iag,iahp,fff,iord,ig,ny,iz,0)
C--     Transform ~H+' --> H+'
        iw(1) = 0
        cc(1) = 1.D0
        iw(3) = idAijk7(3,2,3,itype)                            !as2 Ahq
        cc(3) = assq
        call sqcPMDelta(stor7,iw,cc,ids,idhp,fff,iord,ig,ny,iz)

C--     Get ~G'
        iw(1) = idAijk7(1,1,1,itype)                                 !A0
        cc(1) = 1.D0
        iw(3) = idAijk7(1,1,3,itype)                            !as2 Agg
        cc(3) = assq
        call sqcJump11(stor7,iw,cc,iag,iag,fff,iord,ig,ny,iz,0)
C--     Transform ~G' --> G'
        iw(1) = 0
        cc(1) = 1.D0
        iw(3) = idAijk7(1,2,3,itype)                            !as2 Agq
        cc(3) = assq
        call sqcPMDelta(stor7,iw,cc,ids,idg,fff,iord,ig,ny,iz)

C--     Get S' and E2+',...,En+'
        iw(1) = idAijk7(2,2,1,itype)                                 !A0
        cc(1) = 1.D0
        iw(3) = idAijk7(2,2,3,itype)                            !as2 Aqq
        cc(3) = assq
        do i = 1,nfb
          ia = i*ny+1
          call sqcJump11(stor7,iw,cc,ia,ia,fff,iord,ig,ny,iz,0)
        enddo

C--     Get V' and E2-',...,En-'
        iw(1) = idAijk7(2,2,1,itype)                                 !A0
        cc(1) = 1.D0
        iw(3) = idAijk7(2,2,3,itype)                            !as2 Aqq
        cc(3) = assq
        do i = 7,nfb+6
          ia = i*ny+1
          call sqcJump11(stor7,iw,cc,ia,ia,fff,iord,ig,ny,iz,0)
        enddo

C--   Reverse GSH matching NNLO ----------------------------------------
      elseif(idz.eq.-1 .and. iord.eq.3) then

C--     Get V and E2-,...,En-
        iw(1) = idAijk7(2,2,1,itype)                                 !A0
        cc(1) = 1.D0
        iw(3) = idAijk7(2,2,3,itype)                            !as2 Aqq
        cc(3) = assq
        do i = 7,nfb+6
          ia = i*ny+1
          call sqcJump11(stor7,iw,cc,ia,ia,fff,iord,ig,ny,iz,0)
        enddo

C--     Set H- to zero
        do jy = 0,ny-1
          fff(iahm+jy) = 0.D0
        enddo

C--     Get S and E2+,...,En+
        iw(1) = idAijk7(2,2,1,itype)                                 !A0
        cc(1) = 1.D0
        iw(3) = idAijk7(2,2,3,itype)                            !as2 Aqq
        cc(3) = assq
        do i = 1,nfb
          ia = i*ny+1
          call sqcJump11(stor7,iw,cc,ia,ia,fff,iord,ig,ny,iz,0)
        enddo

C--     Transform G' --> ~G'
        iw(1) = 0
        cc(1) = 1.D0
        iw(3) = idAijk7(1,2,3,itype)                           !-as2 Agq
        cc(3) = -assq
        call sqcPMDelta(stor7,iw,cc,ids,idg,fff,iord,ig,ny,iz)
C--     Get G
        iw(1) = idAijk7(1,1,1,itype)                                 !A0
        cc(1) = 1.D0
        iw(3) = idAijk7(1,1,3,itype)                            !as2 Agg
        cc(3) = assq
        call sqcJump11(stor7,iw,cc,iag,iag,fff,iord,ig,ny,iz,0)

C--     Set H+ to zero
        do jy = 0,ny-1
          fff(iahp+jy) = 0.D0
        enddo

      endif
C--   End of GSH matching ----------------------------------------------

C--   Transform GSH from fff to EPM and put into memory
      call sqcPutGSH(ide,ig,ny,iz2,iaf,fff,0)

      return
      end

C     ================================================
      subroutine sqcMatchPol(itype,mset,ig,iord,iz,ny)
C     ================================================

C--   Forward and reverse matching of polarised pdfs
C--
C--   itype (in): 1=unpol, 2=pol, 3=timelike
C--   mset  (in): pdf set
C--   ig    (in): subgrid index 1,...,nyg2
C--   iord  (in): 1=LO,  2=NLO, 3=NNLO
C--   iz    (in): iz point before matching
C--   ny    (in): upper index yloop

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension fff(14*mxx0)
      dimension ide(13),iaf(14)

      Logical Lbelow

C--   Only for unpolarised dynamic heavy flavours
      if(itype.ne.2) stop 'sqcMatchPol wrong type'

      jord = iord                                !avoid compiler warning

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C--   Get global ids
      do i = 0,12
        ide(i+1) = iqcIdPdfLtoG(mset,i)
      enddo
C--   Transform EPM from memory to GSH and put into fff
      call sqcGetGSH(ide,ig,ny,iz1,iaf,fff,0)

C--   Transform GSH from fff to EPM and put into memory
      call sqcPutGSH(ide,ig,ny,iz2,iaf,fff,0)

      return
      end

C     ================================================
      subroutine sqcMatchTml(itype,mset,ig,iord,iz,ny)
C     ================================================

C--   Forward and reverse matching of time-like pdfs
C--
C--   itype (in): 1=unpol, 2=pol, 3=timelike
C--   mset  (in): pdf set
C--   ig    (in): subgrid index 1,...,nyg2
C--   iord  (in): 1=LO,  2=NLO, 3=NNLO
C--   iz    (in): iz point before matching
C--   ny    (in): upper index yloop

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension iw(3), cc(3)

      dimension fff(14*mxx0)
      dimension ide(13),iaf(14)

      Logical Lbelow

C--   Only for unpolarised dynamic heavy flavours
      if(itype.ne.3) stop 'sqcMatchTml wrong type'

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C--   Get global ids
      do i = 0,12
        ide(i+1) = iqcIdPdfLtoG(mset,i)
      enddo
C--   Transform EPM from memory to GSH and put into fff
      call sqcGetGSH(ide,ig,ny,iz1,iaf,fff,0)

C--   Alphas to use
      iset = 0
      as   = stor7(iqcIaAtab(iza,0,iset))
      assq = as*as

C--   Identifiers and addresses
      idg   = 0
      idh   = nfa

C--   Only NLO, thank you
      iw(1) = 0
      cc(1) = 1.D0
      iw(3) = 0
      cc(3) = 1.D0

      if(idz.eq.+1 .and. iord.ge.2) then
C--     Forward matching
C--     Get H+'
        iw(2) = idAijk7(3,1,2,itype)                             !as Ahg
        cc(2) = as
        call sqcPMDelta(stor7,iw,cc,idg,idh,fff,iord,ig,ny,iz)
      elseif(idz.eq.-1 .and. iord.ge.2) then
C--     Reverse matching
C--     Get H+
        iw(2) = idAijk7(3,1,2,itype)                             !as Ahg
        cc(2) = -as
        call sqcPMDelta(stor7,iw,cc,idg,idh,fff,iord,ig,ny,iz)
      endif

C--   Transform GSH from fff to EPM and put into memory
      call sqcPutGSH(ide,ig,ny,iz2,iaf,fff,0)

      return
      end

C     =========================================================
      subroutine sqcPMDelta(ww,iw,cc,idi,ido,fff,iord,ig,ny,iz)
C     =========================================================

C--   Computes Fout = Fout + c Aij * Fin
C--
C--   ww    (in) : weight store
C--   iw(3) (in) : weight table identifiers Aij(k) for each order k
C--   cc(3) (in) : coefficient of Aij(k) for each order k
C--   idi   (in) : identifier of Fin
C--   ido   (in) : identifier of Fout
C--   fff   (in) : pdf store
C--   iord  (in) : order
C--   ig    (in) : subgrid index
C--   ny    (in) : number of y points
C--   iz    (in) : threshold point
C--

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension ww(*), iw(*), cc(*), fff(*)

      dimension jw(1,1,3), cj(1,1,3), jai(1), jao(1), jad(1)

      logical Lbelow

      jw(1,1,1) = iw(1)
      jw(1,1,2) = iw(2)
      jw(1,1,3) = iw(3)
      cj(1,1,1) = cc(1)
      cj(1,1,2) = cc(2)
      cj(1,1,3) = cc(3)

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

      idd   = 13                                          !scratch space
      iai   = idi*ny+1
      iao   = ido*ny+1
      iad   = idd*ny+1
      jai(1) = iai
      jao(1) = iao
      jad(1) = iad
      call sqcJump11(ww,iw,cc,iai,iad,fff,iord,ig,ny,izb,0)
      do iy = 1,ny
        fff(iao) = fff(iao) + fff(iad)
        iai      = iai+1
        iao      = iao+1
        iad      = iad+1
      enddo

      return
      end

C     =============================================================
      subroutine sqcJump11(ww,iw,cc,iai,iao,fff,iord,ig,ny,iz,iter)
C     =============================================================

C--   1x1 matching: does forward or reverse depending on iz

C--   ww         (in)  weight store
C--   iw(3)      (in)  weight table identifiers Aij(k) for each order k
C--   cc(3)      (in)  coefficient of Aij(k) for each order k
C--   iai        (in)  input  pdf address in fff
C--   iao        (in)  output pdf address in fff
C--   fff        (in)  pdf store
C--   iord       (in)  perturbative order
C--   ig         (in)  subgrid index
C--   ny         (in)  upper index y-subgrid after cuts
C--   iz         (in)  threshold point
C--   n          (in)  number of coupled equations
C--
C--   Remark: iai and iao can be the same

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension ww(*), iw(*), cc(*), fff(*)

C--   Working arrays in this routine
C--
C--   vvv       linear array with weight Aij
C--   aaa       linear array with spline coefficients

      dimension vvv(mxx0),aaa(mxx0)

      logical Lbelow

*mb      write(6,'('' Jump11 ig '', I5)') ig

C--   Initialize
      do i = 1,mxx0
        vvv(i) = 0.D0
        aaa(i) = 0.D0
      enddo

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C--   t-index irrelevant for weights type 1 or 2
      it1 = itfiz5(izb)

C--   V = weight at it1
      do k = 1,iord
        if(iw(k).ne.0) then
          as = cc(k)
          ia = iqcGaddr(ww,1,it1,nfb,ig,iw(k))-1
          do iy = 1,ny
            ia = ia+1
            vvv(iy) = vvv(iy) + as*ww(ia)
          enddo
        endif
      enddo

C--   Now do the jumps
      if(Lbelow) then
C--     Upward jump convert Fin to A and multiply Va = Fout
        call sqcFMatch11(vvv,fff(iai),fff(iao),ny)
*        call sqcNNFjtoAj(fff(iai),aaa,ny)
*        call sqcNSmult(vvv,ny,aaa,fff(iao),ny)
      else
C--     Downward jump solve Va = Fin and convert A to Fout
        call sqcBMatch11(vvv,fff(iao),fff(iai),ny,iter)
*        call sqcNSeqs(vvv,ny,aaa,fff(iai),ny)
*        call sqcNNAjtoFj(aaa,fff(iao),ny)
      endif

      return
      end

C     ===================================
      subroutine sqcFMatch11(vv,ff,fp,ny)
C     ===================================

C--   1x1 forward matching
C--
C--   ff      (in) : f  below threshold
C--   fp     (out) : f' above threshold

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension vv(*), ff(*), fp(*)
      dimension aa(mxx0)

      call sqcNNFjtoAj(ff,aa,ny)
      call sqcNSmult(vv,ny,aa,fp,ny)

      return
      end

C     ========================================
      subroutine sqcBMatch11(vv,ff,fp,ny,iter)
C     ========================================

C--   1x1 reverse matching of F' to F
C--
C--   ff  (out) : F  below threshold
C--   fp   (in) : F' above threshold
C--   iter (in) : 0 = no iteration, otherwise max number of steps
C--
C--   NB: ff(out) can be the same array as fp(in) --> in-place matching

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dimension vv(*)
      dimension ff(*),fp(*)
      dimension aa(mxx0), dd(mxx0)
      dimension dp(mxx0), gg(mxx0)

      if(iter.le.0) then
C--     Try to solve the 2x2 matching equation directly
        call sqcNSeqs(vv,ny,aa,fp,ny)
        call sqcNNAjtoFj(aa,ff,ny)
       else
C--     Iteration will not work when ff and fp are the same array
C--     Thus we use a local output array gg and copy later to ff
C--     Substract A0: V --> W
        do i = 1,2
          vv(i) = vv(i) - smaty2(i,ioy2)
         enddo
C--     First guess F = F'
        call smb_Vcopy(fp,gg,ny)
C--     Iterate
        vnrem = 1.D11
        do i = 1,iter
C--       Get D = W * F
          call sqcFMatch11(vv,gg,dd,ny)
C--       Converged?
C--       Delta = F' - F - D
          call smb_VminV(fp,gg,dp,ny)
          call smb_VminV(dp,dd,dp,ny)
          vnorm = dmb_Vnorm(2,dp,ny)
          if(vnorm.ge.vnrem) goto 100
          vnrem = vnorm
C--       Get F = F' - D
          call smb_VminV(fp,dd,gg,ny)
        enddo
100     continue
C--     Copy local array gg to output array ff
        call smb_Vcopy(gg,ff,ny)

C--     Add A0: W --> V
        do i = 1,2
          vv(i) = vv(i) + smaty2(i,ioy2)
        enddo
      endif

      return
      end


C     =============================================================
      subroutine sqcJump22(ww,iw,cc,iai,iao,fff,iord,ig,ny,iz,iter)
C     =============================================================

C--   2x2 matching: does forward or reverse depending on iz

C--   ww         (in)  weight store
C--   iw(i,j,k)  (in)  weight table identifiers of Aij for each order k
C--   cc(i,j,k)  (in)  coefficient of Aij for each order k
C--   iai(i)     (in)  input  pdf addresses in fff, for i = 1,2
C--   iao(i)     (in)  output pdf addresses in fff, for i = 1,2
C--   fff        (in)  pdf store
C--   iord       (in)  perturbative order
C--   ig         (in)  subgrid index
C--   ny         (in)  upper index y-subgrid after cuts
C--   iz         (in)  threshold point
C--   iter       (in)  0=no iterative reverse matching, otherwise #steps
C--
C--   Remark: iai and iao can be the same

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension ww(*), fff(*), iai(*), iao(*)
      dimension iw(2,2,*), cc(2,2,*)

C--   Working arrays in this routine
C--   vvv       linear array with weights Aij
C--   ivv(i,j)  address of first word of Aij in vvv

      dimension vvv(4*mxx0)
      dimension ivv(2,2)

      logical Lbelow

*mb      write(6,'('' Jump22 ig '', I5)') ig

C--   Initialize
      do i = 1,4*mxx0
        vvv(i) = 0.D0
      enddo
      do j = 1,2
        do i = 1,2
          ivv(i,j) = 0
        enddo
      enddo

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C--   t-index irrelevant for weights type 1 or 2
      it1 = itfiz5(izb)

C--   V = weight at it1
      iw0  = 0
      do i = 1,2
        do j = 1,2
          do k = 1,iord
            if(iw(i,j,k).ne.0) then
              as = cc(i,j,k)
              ia = iqcGaddr(ww,1,it1,nfb,ig,iw(i,j,k))-1
              do iy = 1,ny
                ia = ia+1
                vvv(iw0+iy) = vvv(iw0+iy) + as*ww(ia)
              enddo
            endif
          enddo
          ivv(i,j) = iw0+1
          iw0      = iw0+ny
        enddo
      enddo

C--   Now do the jumps
      iv1 = ivv(1,1)
      iv2 = ivv(1,2)
      iv3 = ivv(2,1)
      iv4 = ivv(2,2)
      ii1 = iai(1)
      ii2 = iai(2)
      io1 = iao(1)
      io2 = iao(2)
      if(Lbelow) then
C--     Upward jump
        call sqcFMatch22(vvv(iv1),vvv(iv2),vvv(iv3),vvv(iv4),
     +                   fff(ii1),fff(ii2),fff(io1),fff(io2),ny)
      else
C--     Reverse jump
        call sqcBMatch22(vvv(iv1),vvv(iv2),vvv(iv3),vvv(iv4),
     +                   fff(io1),fff(io2),fff(ii1),fff(ii2),ny,iter)
      endif

      return
      end

C     ==========================================================
      subroutine sqcFMatch22(vv1,vv2,vv3,vv4,ff1,ff2,fp1,fp2,ny)
C     ==========================================================

C--   2x2 forward matching
C--
C--   ff1,2   (in) : f  below threshold
C--   fp1,2  (out) : f' above threshold

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension vv1(*), vv2(*), vv3(*), vv4(*)
      dimension ff1(*), ff2(*), fp1(*), fp2(*)
      dimension aa1(mxx0), aa2(mxx0)

      call sqcNNFjtoAj(ff1,aa1,ny)
      call sqcNNFjtoAj(ff2,aa2,ny)

      call sqcSGmult(vv1,vv2,vv3,vv4,ny,aa1,aa2,fp1,fp2,ny)

      return
      end

C     ===============================================================
      subroutine sqcBMatch22(vv1,vv2,vv3,vv4,ff1,ff2,fp1,fp2,ny,iter)
C     ===============================================================

C--   2x2 reverse matching of F' to F
C--
C--   ff1,2  (out) : F  below threshold
C--   fp1,2   (in) : F' above threshold
C--   iter    (in) : 0 = no iteration, otherwise max number of steps
C--
C--   NB: ff(out) can be the same array as fp(in) --> in-place matching

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dimension vv1(*), vv2(*), vv3(*), vv4(*)
      dimension ff1(*), ff2(*), fp1(*), fp2(*)
      dimension aa1(mxx0), aa2(mxx0), dd1(mxx0), dd2(mxx0)
      dimension dp1(mxx0), dp2(mxx0), gg1(mxx0), gg2(mxx0)

      if(iter.le.0) then
C--     Try to solve the 2x2 matching equation directly
        call sqcQSGeqs(vv1,vv2,vv3,vv4,aa1,aa2,fp1,fp2,ny)
        call sqcNNAjtoFj(aa1,ff1,ny)
        call sqcNNAjtoFj(aa2,ff2,ny)
      else
C--     Iteration will not work when ff and fp are the same array
C--     Thus we use a local output array gg and copy later to ff
C--     Substract A0: V --> W
        do i = 1,2
          vv1(i) = vv1(i) - smaty2(i,ioy2)
          vv4(i) = vv4(i) - smaty2(i,ioy2)
        enddo
C--     First guess F = F'
        call smb_Vcopy(fp1,gg1,ny)
        call smb_Vcopy(fp2,gg2,ny)
C--     Iterate
        vnrem = 1.D11
        do i = 1,iter
C--       Get D = W * F
          call sqcFMatch22(vv1,vv2,vv3,vv4,gg1,gg2,dd1,dd2,ny)
C--       Converged?
C--       Delta = F' - F - D
          call smb_VminV(fp1,gg1,dp1,ny)
          call smb_VminV(fp2,gg2,dp2,ny)
          call smb_VminV(dp1,dd1,dp1,ny)
          call smb_VminV(dp2,dd2,dp2,ny)
          vnorm = dmb_Vnorm(2,dp1,ny)+dmb_Vnorm(2,dp2,ny)
          if(vnorm.ge.vnrem) goto 100
          vnrem = vnorm
C--       Get F = F' - D
          call smb_VminV(fp1,dd1,gg1,ny)
          call smb_VminV(fp2,dd2,gg2,ny)
        enddo
100     continue
C--     Copy local array gg to output array ff
        call smb_Vcopy(gg1,ff1,ny)
        call smb_Vcopy(gg2,ff2,ny)

C--     Add A0: W --> V
        do i = 1,2
          vv1(i) = vv1(i) + smaty2(i,ioy2)
          vv4(i) = vv4(i) + smaty2(i,ioy2)
        enddo
      endif

      return
      end

C     ===============================================================
      subroutine sqcJumpNN(ww,iw,cc,idim,iai,iao,fff,iord,ig,ny,iz,n)
C     ===============================================================

C--   NxN matching: does forward or reverse depending on iz

C--   ww         (in)  weight store
C--   iw(i,j,k)  (in)  weight table identifiers of Aij for each order k
C--   cc(i,j,k)  (in)  coefficient of Aij for each order k
C--   idim       (in)  first 2 dims of iw and cc
C--   iai(i)     (in)  input  pdf addresses in fff, for i = 1,n
C--   iao(i)     (in)  output pdf addresses in fff, for i = 1,n
C--   fff        (in)  pdf store
C--   iord       (in)  perturbative order
C--   ig         (in)  subgrid index
C--   ny         (in)  upper index y-subgrid after cuts
C--   iz         (in)  threshold point
C--   n          (in)  number of coupled equations
C--
C--   Remark: iai and iao can be the same

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension ww(*), fff(*), iai(*), iao(*)
      dimension iw(idim,idim,*), cc(idim,idim,*)

C--   Working arrays in this routine
C--
C--   vvv       linear array with weights Aij
C--   ivv(i,j)  address of first word of Aij in vvv
C--   aaa       linear array with spline coefficients
C--   iaa(i)    address of first word of pdf(i) in aaa and bbb

      dimension vvv(mce0*mce0*mxx0),aaa(mce0*mxx0)
      dimension ivv(mce0,mce0),iaa(mce0)

      logical Lbelow

*mb      write(6,'('' JumpNN ig '', I5)') ig

C--   Initialize
      do i = 1,mce0*mce0*mxx0
        vvv(i) = 0.D0
      enddo
      do i = 1,mce0*mxx0
        aaa(i) = 0.D0
      enddo
      do j = 1,mce0
        iaa(j) = 0
        do i = 1,mce0
          ivv(i,j) = 0
        enddo
      enddo

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C--   t-index irrelevant for weights type 1 or 2
      it1 = itfiz5(izb)

C--   V = weight at it1
      iw0  = 0
      ia0  = 0
      do i = 1,n
        do j = 1,n
          do k = 1,iord
            if(iw(i,j,k).ne.0) then
              as = cc(i,j,k)
              ia = iqcGaddr(ww,1,it1,nfb,ig,iw(i,j,k))-1
              do iy = 1,ny
                ia = ia+1
                vvv(iw0+iy) = vvv(iw0+iy) + as*ww(ia)
              enddo
            endif
          enddo
          ivv(i,j) = iw0+1
          iw0      = iw0+ny
        enddo
        iaa(i) = ia0+1
        ia0    = ia0+ny
      enddo

C--   Now do the jumps
      if(Lbelow) then
C--     Upward jump convert Fin to A and multiply Va = Fout
        do i = 1,n
          call sqcNNFjtoAj(fff(iai(i)),aaa(iaa(i)),ny)
        enddo
        call sqcNNmult(vvv,ivv,aaa,iaa,fff,iao,n,ny,ny,mce0)
      else
C--     Downward jump solve Va = Fin and convert A to Fout
        call sqcNNeqs(vvv,ivv,aaa,iaa,fff,iai,n,ny,mce0,ierr)
        if(ierr.ne.0) stop 'sqcJumpNN error sqcNNeqs'
        do i = 1,n
          call sqcNNAjtoFj(aaa(iaa(i)),fff(iao(i)),ny)
        enddo
      endif

      return
      end

C     ===============================================
      subroutine sqcGetGSH(ide,ig,ny,iz,iaf,fff,iepm)
C     ===============================================

C--   Read Ei pdfs of subgrid ig from central memory, transform to the
C--   GSH basis, and store the result in a local array
C--   For GSH see writeup appendix 'Forward and Reverse Matching'

C --  ide(13)      (in) : global ids of E0 (gluon) ... E12
C--   ig           (in) : subgrid index
C--   ny           (in) : upper index y-subgrid after cuts
C--   iz           (in) : threshold point
C--   iaf(14)     (out) : addresses of GSH pdfs in fff, 14=scratch
C--   fff(14*ny)  (out) : linear array with GSH pdfs
C--   iepm         (in) : 1 = get EPM (no transformation to GSH)

C--   With n the number of flavours above threshold the SGH indexing is
C--            1  2  ...  n+1 ... 8 ... n+7 ... 13
C--            G  S  ...   H  ... V ...  H  ... 13

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      dimension ide(*), iaf(*), fff(*)

      logical Lbelow

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

C--   Copy bin iq = -ig of Ei basis to fff, for subgrid ig
      ia0 = 1
      do i = 1,13
        call sqcNNgetVj(stor7,ide(i),-ig,ig,ny,fff(ia0))
        iaf(i) = ia0
        ia0    = ia0+ny
      enddo
      iaf(14) = ia0                                         !scratch pdf

      if(Lbelow .or. iepm.eq.1) return

C--   Transform Singlet (2) and Eh+ (nf+1) to GSH
      ias = iaf(2)-1
      iah = iaf(nfa+1)-1
      do iy = 1,ny
        ias = ias+1
        iah = iah+1
        fff(iah) = (fff(ias)-fff(iah))/nfa
        fff(ias) =  fff(ias)-fff(iah)
      enddo
C--   Transform Valence (8) and Eh- (nf+7) to GSH
      iav = iaf(8)-1
      iah = iaf(nfa+7)-1
      do iy = 1,ny
        iav = iav+1
        iah = iah+1
        fff(iah) = (fff(iav)-fff(iah))/nfa
        fff(iav) =  fff(iav)-fff(iah)
      enddo

      return
      end

C     ===============================================
      subroutine sqcPutGSH(ide,ig,ny,iz,iaf,fff,iepm)
C     ===============================================

C--   Read GSH pdfs from a local array, transform to the Ei basis
C--   and store the result for subgrid ig in central memory
C--   For GSH see writeup appendix 'Forward and Reverse Matching'

C --  ide(13)     (in) : global ids of E0 (gluon) ... E12
C--   ig          (in) : subgrid index
C--   ny          (in) : upper index y-subgrid after cuts
C--   iz          (in) : threshold point
C--   iaf(14)     (in) : addresses of GSH pdfs in fff, 14=scratch
C--   fff(14*ny)  (in) : linear array with GSH pdfs
C--   iepm        (in) : 1 = put EPM (no transformation from GSH)

C--   With n the number of flavours above threshold the SGH indexing is
C--            1  2  ...  n+1 ... 8 ... n+7 ... 13
C--            G  S  ...   H  ... V ...  H  ... 13

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      dimension ide(*), iaf(*), fff(*)

      logical Lbelow

C--   All you want to know about iz
      call sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)

      if(.not.Lbelow .and. iepm.ne.1) then
C--     Transform GSH to Singlet (2) and Eh+ (nf+1)
        ias = iaf(2)-1
        iah = iaf(nfa+1)-1
        do iy = 1,ny
          ias = ias+1
          iah = iah+1
          fff(ias) =  fff(ias)+fff(iah)
          fff(iah) =  fff(ias)-nfa*fff(iah)
        enddo
C--     Transform GSH to Valence (8) and Eh- (nf+7)
        iav = iaf(8)-1
        iah = iaf(nfa+7)-1
        do iy = 1,ny
          iav = iav+1
          iah = iah+1
          fff(iav) =  fff(iav)+fff(iah)
          fff(iah) =  fff(iav)-nfa*fff(iah)
        enddo
      endif

C--   Copy fff to bin iq = -ig of Ei basis, for subgrid g
      do i = 1,13
        call sqcNNputVj(stor7,ide(i),-ig,ig,ny,fff(iaf(i)))
        call sqcNNputVj(stor7,ide(i), 0 ,ig,ny,fff(iaf(i)))
      enddo

      return
      end

C     ===========================================================
      subroutine sqcWhatIz(iz,iz1,iz2,idz,izb,iza,nfb,nfa,Lbelow)
C     ===========================================================

C--   Get threshold characteristics useful for matching routines

C--   iz       (in) : iz point (--> stop if not at threshold)
C--   iz1     (out) : start point of match
C--   iz2     (out) : end   point of match
C--   idz     (out) : +1 (-1) for forward (reverse) match
C--   izb     (out) : iz below threshold
C--   iza     (out) : iz above threshold
C--   nfb     (out) : nf below threshold
C--   nfa     (out) : nf above threshold
C--   Lbelow  (out) : .true. if iz below threshold

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      logical Lbelow

      nf    = itfiz5(-iz)
      nfzp1 = itfiz5(-(iz+1))
      nfzm1 = itfiz5(-(iz-1))

      if    (iz.ne.nzz2 .and. nfzp1.eq.nf+1) then
        Lbelow  = .true.
        iz1     = iz
        iz2     = iz+1
        idz     = +1
        izb     = iz
        iza     = iz+1
        nfb     = nf
        nfa     = nf+1
      elseif(iz.ne.1    .and. nfzm1.eq.nf-1) then
        Lbelow  = .false.
        iz1     = iz
        iz2     = iz-1
        idz     = -1
        izb     = iz-1
        iza     = iz
        nfb     = nf-1
        nfa     = nf
      else
        stop 'sqcWhatIz: iz not at threshold'
      endif

      return
      end







