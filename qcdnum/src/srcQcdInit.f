
C--   This is the file srcQcdInit.f containing the qcdnum init routines

C--   subroutine sqcIniCns
C--   subroutine sqcSetLun(lun,fname)
C--   subroutine sqcBanner(lun)
C--   subroutine sqcReftoo(lun)
C--   subroutine sqcIniStore(nwlast,ierr)

C=======================================================================
C==   Initialization routines ==========================================
C=======================================================================

C     ====================
      subroutine sqcIniCns
C     ====================

C---  Initialize constants.
C---  Called by sqc_qcinit.
 
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc' 
      include 'qconst.inc'
      include 'qvers1.inc'
      include 'qgrid2.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'pstor8.inc'
      include 'qfast9.inc'
      include 'qcard9.inc'

C--   Fixed parameters
C--   ----------------

      pi     = 3.14159265359
      proton = 0.9382796
      eutron = 0.9395731
      ucleon = (proton + eutron) / 2.

C--   Constants inhereted from Ouarou and Virchaux (original QCDNUM) 
      c1s3   = 1./3.
      c2s3   = 2./3.
      c4s3   = 4./3.
      c5s3   = 5./3.
      c8s3   = 8./3.
      c14s3  = 14./3.
      c16s3  = 16./3.
      c20s3  = 20./3.
      c28s3  = 28./3.
      c38s3  = 38./3.
      c40s3  = 40./3.
      c44s3  = 44./3.
      c52s3  = 52./3.
      c136s3 = 136./3.
      c11s6  = 11./6.
      c2s9   = 2./9.
      c4s9   = 4./9.
      c10s9  = 10./9.
      c14s9  = 14./9.
      c16s9  = 16./9.
      c40s9  = 40./9.
      c44s9  = 44./9.
      c62s9  = 62./9.
      c112s9 = 112./9.
      c182s9 = 182./9.
      c11s12 = 11./12.
      c35s18 = 35./18.
      c11s3  = 11./3.
      c22s3  = 22./3.
      c61s12 = 61./12.
      c215s1 = 215./12.
      c29s12 = 29./12.
      cpi2s3 = pi**2/3.
      cpia   = 67./18. - cpi2s3/2.
      cpib   = 4.*cpi2s3
      cpic   = 17./18. + 3.5*cpi2s3
      cpid   = 367./36. - cpi2s3
      cpie   = 5. - cpi2s3
      cpif   = cpi2s3 - 218./9.

      cca    = 3.
      ccf    = (cca*cca-1.)/(2.*cca)
      ctf    = 0.5
      catf   = cca*ctf
      cftf   = ccf*ctf

C--   Initialize storage
C--   ------------------
      do i = 1,nwf0
        stor7(i) = 0.D0
      enddo
      do i = 1,nwp0
        pars7(i) = 0.D0
      enddo
      Lqswrite6 = .true.
      do i = 1,mqs0
        qstore6(i) = 0.D0
      enddo

C--   Parameter lists and version
      ipver6 = 1
      itver6 = 1
      do j = 0,mpl0
        do i = 1,mpar0
          parlist6(i,j) = 0.D0
        enddo
      enddo
      do i = 1,mpl0
        keycount6(i) = 0
      enddo

C--   Pdf set with id = -1 is base set with scratch pdfs and alfas tables
      do i = -1,mset0
        isetf7(i)  = 0
        Lfill7(i)  = .false.
        ikeyf7(i)  = 0
      enddo
C--   Base set
      ikeyf7(0) = 1
      
C--   No grid available
      Lygrid2 = .false.
      Ltgrid2 = .false.
      igver2  =  0
C--   Weight tables not initialised
      Lwtini7 = .false.
C--   No fast structure functions
      nmax9   = 0
      nnff9   = 0
C--   No parameter lists available
      nplist6 = 0
C--   Pdf scope
      Lscopechek6 = .false.
      iscopeslot6 = 1                      !current parameters

C--   Set parameter bits
      ipbits8 = 0
      call smb_sbit1(ipbits8,infbit8)      !no nfmap
      call smb_sbit1(ipbits8,iasbit8)      !no alfas tables
      call smb_sbit1(ipbits8,izcbit8)      !no iz cuts
      call smb_sbit1(ipbits8,ipsbit8)      !no param store

C--   Blank addon package subroutine name (smb_cfill comes from mbutil)
C--   -----------------------------------------------------------------
      call smb_cfill(' ',usrnam3)

C--   Default parameters which can be changed by the user
C--   ---------------------------------------------------

C--   Default values for qpars6
      aepsi6   = 1.D-9        !tolerance floating-point comparison
      gepsi6   = 1.D-7        !requested accuracy Gauss integration
      dflim6   = 0.5D0        !max abs deviation of spline from function
      aslim6   = 10.D0        !largest allowed value of alphas
      qnull6   = 1.D11        !qcdnum null value
      qlimd6   = 0.1D0        !lowest allowed mu2 value
      qlimu6   = 1.D11        !largest allowed mu2 value
      niter6   = 1            !number of dnward evolution iterations
      itlmc6   = 1            !apply timelike matching conditions
      idbug6   = 0            !no debug printout
      inew6    = 1            !use new matching code

C--   Default FFNS with nf = 3
      nfix6     = 3
      call sqcThrFFNS(nfix6)

      q0alf6 = 8315.25D0      !MZ^2
      alfq06 = 0.118D0        !alphas(MZ^2)

C--   Default for renormalization scale r2 = aar6*f2 + bbr6
      aar6   = 1.D0
      bbr6   = 0.D0
      itmin6 = 1

C--   Default is NLO
      iord6  = 2

C--   Predefined keys
C--                 123456789012
      qkeys9( 1) = 'SETLUN  QKEY'
      qkeys9( 2) = 'SETVAL  QKEY'
      qkeys9( 3) = 'SETINT  QKEY'
      qkeys9( 4) = 'GXMAKE  QKEY'
      qkeys9( 5) = 'GQMAKE  QKEY'
      qkeys9( 6) = 'FILLWT  QKEY'
      qkeys9( 7) = 'SETORD  QKEY'
      qkeys9( 8) = 'SETALF  QKEY'
      qkeys9( 9) = 'SETCBT  QKEY'
      qkeys9(10) = 'MIXFNS  QKEY'
      qkeys9(11) = 'SETABR  QKEY'
      qkeys9(12) = 'SETCUT  QKEY'
      qkeys9(13) = 'QCSTOP  QKEY'
      do i = 14,mky0
        qkeys9(i) = '        FREE'
      enddo

C--   Evdglap initialisation
C--   nnopt6(1,2,3,...) = # perturbative terms at order 1,2,3,...
C--   inopt6 is a n-digit encoding e.g. 123 = 1,2,3 terms at order 1,2,3
      nnopt6(0) = mord0
      inopt6    = 0
      ifac      = 10**mord0
      do i = 1,mord0
        ifac      = ifac/10
        inopt6    = inopt6 + ifac*i
        nnopt6(i) = i
      enddo

      return
      end

C     ===============================
      subroutine sqcSetLun(lun,fname)
C     ===============================

C--   Redirect QCDNUM output

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*(*) fname

      lunerr1 = lun
      if(lun.ne.6) then
        open(unit=lun,file=fname,status='unknown')
      endif

      return
      end

C     =========================
      subroutine sqcBanner(lun)
C     =========================

C--   QCDNUM banner printout

      implicit double precision (a-h,o-z)

      include 'qvers1.inc'
      
      write(lun,'('' '')')
      write(lun,'(''                  ///                 '',
     &            ''                 .().                 '')')
      write(lun,'(''                 (..)                 '',
     &            ''                 (--)                 '')')
      write(lun,'(''  +----------ooO--()--Ooo-------------'',
     &            ''-------------ooO------Ooo---------+   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  |    #####      ######    ######    '',
     &            '' ##    ##   ##    ##   ##     ##  |   '')')
      write(lun,'(''  |   ##   ##    ##    ##   ##   ##   '',
     &            '' ###   ##   ##    ##   ###   ###  |   '')')
      write(lun,'(''  |  ##     ##   ##    ##   ##    ##  '',
     &            '' ####  ##   ##    ##   #### ####  |   '')')
      write(lun,'(''  |  ##     ##   ##         ##    ##  '',
     &            '' ## ## ##   ##    ##   ## ### ##  |   '')')
      write(lun,'(''  |  ##     ##   ##         ##    ##  '',
     &            '' ##  ####   ##    ##   ##  #  ##  |   '')')
      write(lun,'(''  |   ##   ##    ##    ##   ##   ##   '',
     &            '' ##   ###   ##    ##   ##     ##  |   '')')
      write(lun,'(''  |    #####      ######    ######    '',
     &            '' ##    ##    ######    ##     ##  |   '')')
      write(lun,'(''  |        ##                         '',
     &            ''                                  |   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  |    Version '',A10,''  '',A8,''   '',
     &            ''      Author m.botje@nikhef.nl    |   '')') 
     &            cvers1,cdate1
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  +-----------------------------------'',
     &            ''----------------------------------+   ''//)')
     
      return
      end
      
C     =========================
      subroutine sqcReftoo(lun)
C     =========================

C--   QCDNUM banner printout
     
     
      write(lun,'(''  +-----------------------------------'',
     &            ''----------------------------------+   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  |    If you use QCDNUM, please refer'',
     &            '' to:                              |   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  |    M. Botje, Comput. Phys. Commun.'',
     &            '' 182(2011)490, arXiV:1005.1481    |   '')')
*      write(lun,'(''  |                                   '',
*     &            ''     Erratum   arXiV:1602.08383   |   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  +-----------------------------------'',
     &            ''----------------------------------+   '')')
    

      return
      end

C=======================================================================
C===  Base set 0 (parameter list, alfas tables, scratch pdfs  ==========
C=======================================================================


C     ===================================
      subroutine sqcIniStore(nwlast,ierr)
C     ===================================

C--   Initialize store (called by grid routines after grid is defined)
C--
C--   nwlast  (out) :  # memeory words used (< 0 no space)
C--   ierr   (out)  0 = OK
C--                -1 = empty set of tables (never occurs)
C--                -2 = not enough space
C--                -3 = iset count MST0 exceeded
C--                -4 = iset exist but with smaller n or different ifrst

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension itypes(6)

C--   Initialise
      call smb_IFill(itypes,6,0)
C--   Book pdf and alfas tables of set 0; do nothing if tables exist
      iset  =  0
      ntab  =  nwrk0
      ifrst =  1
      noalf =  0
      call sqcPdfBook(iset,ntab,ifrst,noalf,nwlast,ierr)
C--   Lfill7 has no meaning for this set but better set it to true
      Lfill7(iset) = .true.
C--   Set increments for fast addressing in pdf table (type 5)
      idg    = iqcIdPdfLtoG(iset,ifrst)
      iag    = iqcG5ijk(stor7,1,1,idg)
      inciy7 = iqcG5ijk(stor7,2,1,idg)   - iag
      inciz7 = iqcG5ijk(stor7,1,2,idg)   - iag
      incid7 = iqcG5ijk(stor7,1,1,idg+1) - iag
C--   Set increments for fast addressing in pointer table (type 7)
      idg    = 1000*isetf7(iset)+701
      incpt7 = iqcG7ij(stor7,1,idg+1) - iqcG7ij(stor7,1,idg)
C--   Set alfas table adresses for evolution routines
      do k = 1,mord0
        id = iqcIdAtab(k,iset)
        do j = 1,2
          do i = 1,2
            idEijk7(i,j,k) = id
          enddo
        enddo
      enddo

      return
      end

