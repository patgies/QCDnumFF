 
C--   This is hqweits.f containing the heavy quark weight routines
C--
C--   subroutine hqfillw(istf,qmas,aq2,bq2,nwords)
C--   subroutine hqdumpwCPP(lun,file,ls)
C--   subroutine hqdumpw(lun,file)
C--   subroutine hqreadwCPP(lun,file,ls,nwords,ierr)
C--   subroutine hqreadw(lun,file,nwords,ierr)
C--   subroutine hqindex(ntabs,istf,qmass,ifirst,ntot)
C--   subroutine hqfillL(w,id0,qmas,aq,bq)
C--   subroutine hqfill2(w,id0,qmas,aq,bq)
C--   double precision function dhqAchi(qmu2)
C--    
C--   double precision function dhqC0LG(chi,qmu2,nf)
C--   double precision function dhqC1LG(chi,qmu2,nf)
C--   double precision function dhqC1BLG(chi,qmu2,nf)
C--   double precision function dhqC1LQ(chi,qmu2,nf)
C--   double precision function dhqC1BLQ(chi,qmu2,nf)
C--   double precision function dhqD1LQ(chi,qmu2,nf)
C--   double precision function dhqD1BLQ(chi,qmu2,nf)
C--    
C--   double precision function dhqC02G(chi,qmu2,nf)
C--   double precision function dhqC12G(chi,qmu2,nf)
C--   double precision function dhqC1B2G(chi,qmu2,nf)
C--   double precision function dhqC12Q(chi,qmu2,nf)
C--   double precision function dhqC1B2Q(chi,qmu2,nf)
C--   double precision function dhqD12Q(chi,qmu2,nf)
C--   double precision function dhqD1B2Q(chi,qmu2,nf)
C--

C-----------------------------------------------------------------------
CXXHDR    void hqfillw(int istf, double *qmass, double aq, double bq, int &nused);
C-----------------------------------------------------------------------
CXXHFW  #define fhqfillw FC_FUNC(hqfillw,HQFILLW)
CXXHFW    void fhqfillw(int*, double*, double*, double*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void hqfillw(int istf, double *qmass, double aq, double bq, int &nused)
CXXWRP  {
CXXWRP    fhqfillw(&istf,qmass,&aq,&bq,&nused);
CXXWRP  }
C-----------------------------------------------------------------------
 
C     ============================================      
      subroutine hqfillw(istf,qmas,aq2,bq2,nwords)
C     ============================================

C--   Input   istf     1=FL, 2=F2, 3=both
C--           qmas(3)  c,b,t mass, if < 1GeV skip heavy quark
C--           aq2,bq2  Q2 = a*mu2 + b
C--   Output  nwords   #words used for tables < 0 not enough space

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'hqstf.inc'
      include 'hqstore.inc'
      include 'hqflags.inc'
      include 'hqparms.inc'
      
      dimension qmas(3),ifirst(6),pars(11)
      dimension itypes(6)
      data itypes /0, 0, 0, 0, 0, 0/
      
      character*1 cbt(3)
      data cbt /'c','b','t'/

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('HQFILLW')
C--   Check user input
      if(istf.lt.1 .or. istf.gt.3) stop 
     +              'HQFILLW: input ISTF not in range [1-3]'
C--   Limit allowed range aq2 and bq2
      if(aq2.lt.0.1D0 .or. aq2.gt.10.D0) then
        stop 'HQFILLW: Coefficient AQ outside range [0.1,10]'
      endif
      if(abs(bq2).gt.100.D0) then
        stop 'HQFILLW: Coefficient BQ outside range [-100,100]'
      endif                
C--   Setup table ids: ifirst(6) in order FLc,b,t, F2c,b,t; 0=no table 
      ntabs = 7   !# tables per stf
      call hqindex(ntabs,istf,qmas,ifirst,ntot)
C--   No tables to be generated
      if(ntot.eq.0) stop 'HQFILLW: no tables to be generated'
C--   Set hqstf initialized
      if(ihqini.ne.12345) then
        ihqini = 12345
        ihqpdf = 1
      endif
C--   Copy to common block hqparms
      hqmas(1) = qmas(1)
      hqmas(2) = qmas(2)
      hqmas(3) = qmas(3)
      do i = 1,6
         idtab0(i) = ifirst(i)
      enddo
      hqaaaa = aq2
      hqbbbb = bq2
C--   Partition store
      itypes(3) = -ntot
      np  = 20
      new = 1     !only one table set, thank you
      call MakeTab(hqstor,nhqstor,itypes,np,new,jsetw,nwords)
C--   Store in common block hqstore.inc
      nhused = nwords             
C--   Get QCDNUM logical unit number for message
      call getint('lunq',lun)
      write(lun,'(/'' HQFILLW: start weight calculations'')') 
C--   FL tables c,b,t
      do i = 1,3
        if(ifirst(i).ne.0) then
          call hqfillL(hqstor,ifirst(i),qmas(i),hqaaaa,hqbbbb)
          write(lun,'(''          FL'',A,'' done ...'')') cbt(i) 
        endif
      enddo
C--   F2 tables c,b,t
      do i = 1,3
        if(ifirst(i+3).ne.0) then          
          call hqfill2(hqstor,ifirst(i+3),qmas(i),hqaaaa,hqbbbb)
          write(lun,'(''          F2'',A,'' done ...'')') cbt(i)
        endif
      enddo 
      write(lun,'('' HQFILLW: calculations completed'')')
C--   Store quark masses and table ids
      do i = 1,3
        pars(i) = hqmas(i)
      enddo
      do i = 1,6
        pars(i+3) = idtab0(i)
      enddo
      pars(10) = hqaaaa
      pars(11) = hqbbbb
      call setparw(hqstor,jsetw,pars,11)
C--   Clear subroutine name          
      call clrUmsg
            
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void hqdumpw(int lun, string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fhqdumpwcpp FC_FUNC(hqdumpwcpp,HQDUMPWCPP)
CXXHFW    void fhqdumpwcpp(int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void hqdumpw(int lun, string fname)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fhqdumpwcpp(&lun,cfname,&ls);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================
      subroutine hqdumpwCPP(lun,file,ls)
C     ==================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      character*(100) file

      if(ls.gt.100) stop 'hqdumpwCPP: file name size > 100 characters'

      call hqdumpw(lun,file(1:ls))

      return
      end
 
C     ============================
      subroutine hqdumpw(lun,file)
C     ============================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'hqstf.inc'
      include 'hqstore.inc'

      character*(*) file

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('HQDUMPW')
C--   Write (chvers in hqstf.inc)
      call DumpTab(hqstor,1,lun,file,chvers)
C--   Clear subroutine name
      call clrUmsg     

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void hqreadw(int lun, string fname, int &nused, int &ierr);
C-----------------------------------------------------------------------
CXXHFW  #define fhqreadwcpp FC_FUNC(hqreadwcpp,HQREADWCPP)
CXXHFW    void fhqreadwcpp(int*, char*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void hqreadw(int lun, string fname, int &nused, int &ierr)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fhqreadwcpp(&lun,cfname,&ls,&nused,&ierr);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================================
      subroutine hqreadwCPP(lun,file,ls,nwords,ierr)
C     ==============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      character*(100) file

      if(ls.gt.100) stop 'hqreadwCPP: file name size > 100 characters'

      call hqreadw(lun,file(1:ls),nwords,ierr)

      return
      end

C     ========================================
      subroutine hqreadw(lun,file,nwords,ierr)
C     ========================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'hqstf.inc'
      include 'hqstore.inc'
      include 'hqflags.inc'
      include 'hqparms.inc'

      character*(*) file

      dimension pars(11)

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('HQREADW')
C--   Read (chvers in hqstf.inc)
      new = 1  !always one table set thank you
      call ReadTab(hqstor,nhqstor,lun,file,chvers,new,isetw,nwords,ierr)
      if(ierr.ne.0) return
C--   Store in common block hqstore.inc
      nhused = nwords       
C--   Set hqstf initialized
      if(ihqini.ne.12345) then
        ihqini = 12345
        ihqpdf = 1
      endif
C--   Read quark masses, table ids and scale parameters
      call getparw(hqstor,isetw,pars,11)
      do i = 1,3
        hqmas(i) = pars(i)
      enddo
      do i = 1,6
        idtab0(i) = int(anint(pars(i+3)))
      enddo
      hqaaaa = pars(10)
      hqbbbb = pars(11)  
C--   Clear subroutine name
      call clrUmsg
C--   Store iset
      jsetw = isetw

      return
      end
      
C     ================================================
      subroutine hqindex(ntabs,istf,qmass,ifirst,ntot)
C     ================================================

C--   Assign tables to FLcbt and F2cbt
C--
C--   Input:  ntabs     # tables per stf
C--           istf      1=F2, 2=FL, 3=both
C--           qmass(3)  cbt mass deselect if mass < 1 GeV
C--   Output: ifirst(6) first id in order FLcbt, F2cbt (0=no table)
C--           ntot      total number of tables to generate

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dimension qmass(*), ifirst(*)
      
      id    = 301-ntabs   !base address
      do i = 1,6
        ifirst(i) = 0
      enddo  
      ntot = 0 
      
      do ist = 1,2
        if(ist.eq.1 .and. (istf.eq.1.or.istf.eq.3))  then
          do iq = 1,3
            if(qmass(iq).ge.1.D0) then
              id             = id+ntabs
              ifirst(iq)     = id
              ntot           = ntot+ntabs
            endif
          enddo   
        elseif(ist.eq.2 .and. (istf.eq.2.or.istf.eq.3)) then
          do iq = 1,3
            if(qmass(iq).ge.1.D0) then
              id             = id+ntabs
              ifirst(iq+3)   = id
              ntot           = ntot+ntabs
            endif
          enddo  
        endif        
      enddo
      
      return
      end      
      
C     ====================================
      subroutine hqfillL(w,id0,qmas,aq,bq)
C     ====================================      

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'hqstf.inc'
      include 'hqstore.inc'
 
      dimension w(*)
      
      common /hqpass/ qmass, ascale, bscale
      
      external dhqAchi
      external dhqC0LG
      external dhqC1LG,dhqC1BLG,dhqC1LQ,dhqC1BLQ,dhqD1LQ,dhqD1BLQ

C--   Pass to common block
      qmass  = qmas
      ascale = aq
      bscale = bq

C--   Global identifier
      idg = 1000*jsetw + id0

C--   LO:
      call makeWtA(w,idg  ,dhqC0LG ,dhqAchi)
C--   NLO:
      call makeWtA(w,idg+1,dhqC1LG ,dhqAchi)
      call makeWtA(w,idg+2,dhqC1BLG,dhqAchi)
      call makeWtA(w,idg+3,dhqC1LQ ,dhqAchi)
      call makeWtA(w,idg+4,dhqC1BLQ,dhqAchi)
      call makeWtA(w,idg+5,dhqD1LQ ,dhqAchi)
      call makeWtA(w,idg+6,dhqD1BLQ,dhqAchi)
      
      return
      end
                                                
C     ====================================
      subroutine hqfill2(w,id0,qmas,aq,bq)
C     ====================================      

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'hqstf.inc'
      include 'hqstore.inc'
 
      dimension w(*)
      
      common /hqpass/ qmass, ascale, bscale
      
      external dhqAchi
      external dhqC02G
      external dhqC12G,dhqC1B2G,dhqC12Q,dhqC1B2Q,dhqD12Q,dhqD1B2Q

C--   Pass to common block
      qmass  = qmas
      ascale = aq
      bscale = bq

C--   Global identifier
      idg = 1000*jsetw + id0

C--   LO:
      call makeWtA(w,idg  ,dhqC02G ,dhqAchi)
C--   NLO:
      call makeWtA(w,idg+1,dhqC12G ,dhqAchi)
      call makeWtA(w,idg+2,dhqC1B2G,dhqAchi)
      call makeWtA(w,idg+3,dhqC12Q ,dhqAchi)
      call makeWtA(w,idg+4,dhqC1B2Q,dhqAchi)
      call makeWtA(w,idg+5,dhqD12Q ,dhqAchi)
      call makeWtA(w,idg+6,dhqD1B2Q,dhqAchi)
      
      return
      end
      
C     =======================================
      double precision function dhqAchi(qmu2)
C     =======================================      

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
 
      common /hqpass/ qmass, ascale, bscale
      
C--   Q2 may become negative, depending on ascale and bscale
C--   We put here a Q2 limit of 0.25 GeV2 which then leads
C--   to wrong tables at lower Q2. This does not matter because
C--   in the stf routines the lower Q2 limit is set to 0.5 GeV2.
      
      qsq     = max((ascale*qmu2 + bscale),0.25D0)
      dhqAchi = 1.D0 + 4.D0*qmass*qmass/qsq
      
      return
      end
      
C--   -----------------------------------------------------------
C--   FL coefficient functions      -----------------------------   
C--   -----------------------------------------------------------      
      
C     ==============================================
      double precision function dhqC0LG(chi,qmu2,nf)
C     ==============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
 
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale

      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c0Lg   = 0.
      xi     = 1.D0/factor
      eta    = xi * (1.-x)/(4.D0*x) - 1.
      c0Lg   = C0_Lg(eta,xi)
C--   Since xi = Q2/m2 the line below takes care of the
C--   factor Q2/(2pi m2) in front of the integral
      c0Lg   = c0Lg * xi / (2.D0*pi)
      dhqC0LG = c0Lg/x
 
      return
      end
      
C     ==============================================
      double precision function dhqC1LG(chi,qmu2,nf)
C     ==============================================
      
C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale

      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c1lg   = 0.
      xi     = 1./factor
      eta    = xi * (1.-x)/(4.*x) - 1.
      bet    = sqrt(eta/(1.+eta))
      rho    = 1./(1.+eta)
C--   Eq. (9)+(10) in hep-ph/9411431; CATF in Eq. (10) should be CFTF
      c1lg =   catf * h1_alg(eta,xi) +
     +         cftf * h1_flg(eta,xi) +
     +  catf * bet  * gfun_l(eta,xi) +
     +  catf * rho  * efun_la(eta,xi) +
     +  cftf * rho  * efun_lf(eta,xi)
C--   This takes care of the factor in front of the integral
      c1lg   = c1lg*4.*pi/factor
      dhqC1LG = c1lg/x
 
      return
      end
 
C     ===============================================
      double precision function dhqC1BLG(chi,qmu2,nf)
C     ===============================================
      
C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c1blg  = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi    = 1./factor
        eta   = xi * (1.-x)/(4.*x) - 1.
        bet   = sqrt(eta/(1.+eta))
        rho   = 1./(1.+eta)
C--     Eq. (12) in hep-ph/9411431        
        c1blg = catf * h1bar_lg(eta,xi) +
     +   catf * bet  * gbar_l(eta,xi) +
     +   catf * rho  * ebar_la(eta,xi)
C--     This takes care of the factor in front of the integral     
        c1blg = c1blg*4.*pi/factor
C--     And this of the factor ln(mu2/m2)
        c1blg = c1blg*log(qmu2/(qmass*qmass))        
*mb      endif
      dhqC1BLG = c1blg/x
 
      return
      end
 
C     ==============================================
      double precision function dhqC1LQ(chi,qmu2,nf)
C     ==============================================
      
C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c1lq   = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi   = 1./factor
        eta  = xi * (1.-x)/(4.*x) - 1.
        bet  = sqrt(eta/(1.+eta))
        bet3 = bet*bet*bet
        rho  = 1./(1.+eta)
C--     Eq. (26) in hep-ph/9411431         
        c1lq = cftf * h1_hlq(eta,xi) +
     +  cftf * bet3 * gfun_l(eta,xi)
C--     This takes care of the factor in front of the integral     
        c1lq = c1lq*4.*pi/factor
*mb      endif
      dhqC1LQ = c1lq/x
 
      return
      end
 
C     ===============================================
      double precision function dhqC1BLQ(chi,qmu2,nf)
C     ===============================================
      
C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c1blq  = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi    = 1./factor
        eta   = xi * (1.-x)/(4.*x) - 1.
        bet   = sqrt(eta/(1.+eta))
        bet3  = bet*bet*bet
        rho   = 1./(1.+eta)
C--     Eq. (27) in hep-ph/9411431        
        c1blq = cftf * h1bar_hlq(eta,xi) +
     +   cftf * bet3 * gbar_l(eta,xi)
C--     This takes care of the factor in front of the integral     
        c1blq = c1blq*4.*pi/factor
C--     And this of the factor ln(mu2/m2)
        c1blq = c1blq*log(qmu2/(qmass*qmass))        
*mb      endif
      dhqC1BLQ = c1blq/x
 
      return
      end
 
C     ==============================================
      double precision function dhqD1LQ(chi,qmu2,nf)
C     ==============================================
      
C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      d1lq   = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi   = 1./factor
        eta  = xi * (1.-x)/(4.*x) - 1.
C--     Eq. (28) in hep-ph/9411431        
        if(qsq.le.1.5) then
          d1lq = cftf * h1f_llq(eta,xi)
        else
          d1lq = cftf * h1_llq(eta,xi)
        endif
C--     This takes care of the factor in front of the integral       
        d1lq = d1lq*4.*pi/factor
*mb      endif
      dhqD1LQ = d1lq/x
 
      return
      end
 
C     ===============================================
      double precision function dhqD1BLQ(chi,qmu2,nf)
C     ===============================================
      
C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac      
C--   See p.5 in hep-ph/9411431
      dhqD1BLQ = 0.
 
      return
      end
      
C--   -----------------------------------------------------------
C--   F2 coefficient functions      -----------------------------
C--   -----------------------------------------------------------      
      
C     ==============================================
      double precision function dhqC02G(chi,qmu2,nf)
C     ==============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
 
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c02g   = 0.
*mb      if(x.lt.(1.D0/factor)) then
        xi   = 1.D0/factor
        eta  = xi * (1.-x)/(4.D0*x) - 1.
        c02g = C0_Lg(eta,xi)+C0_Tg(eta,xi)
C--     Since xi = Q2/m2 the line below takes care of the
C--     factor Q2/(2pi m2) in front of the integral        
        c02g = c02g * xi / (2.D0*pi)
*mb      endif
      dhqC02G = c02g/x
 
      return
      end

C     ==============================================
      double precision function dhqC12G(chi,qmu2,nf)
C     ==============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale

      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c12g   = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi   = 1./factor
        eta  = xi * (1.-x)/(4.*x) - 1.
        bet  = sqrt(eta/(1.+eta))
        rho  = 1./(1.+eta)
C--     Eq. (9)+(10) in hep-ph/9411431; CATF in Eq. (10) should be CFTF         
        c12g = catf * (h1_alg(eta,xi)+h1_atg(eta,xi)) +
     +         cftf * (h1_flg(eta,xi)+h1_ftg(eta,xi)) +
     +  catf * bet  * (gfun_l(eta,xi)+gfun_t(eta,xi)) +
     +  catf * rho  * (efun_la(eta,xi)+efun_ta(eta,xi)) +
     +  cftf * rho  * (efun_lf(eta,xi)+efun_tf(eta,xi))
C--     This takes care of the factor in front of the integral
        c12g = c12g*4.*pi/factor
*mb      endif
      dhqC12G = c12g/x
 
      return
      end
 
C     ===============================================
      double precision function dhqC1B2G(chi,qmu2,nf)
C     ===============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'ffcons.inc'
            
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c1b2g  = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi    = 1./factor
        eta   = xi * (1.-x)/(4.*x) - 1.
        bet   = sqrt(eta/(1.+eta))
        rho   = 1./(1.+eta)
C--     Eq. (12) in hep-ph/9411431        
        c1b2g = catf * (h1bar_lg(eta,xi)+h1bar_tg(eta,xi)) +
     +   catf * bet  * (gbar_l(eta,xi)+gbar_t(eta,xi)) +
     +   catf * rho  * (ebar_la(eta,xi)+ebar_ta(eta,xi))
C--     This takes care of the factor in front of the integral
        c1b2g = c1b2g*4.*pi/factor
C--     And this of the factor ln(mu2/m2)
        c1b2g = c1b2g*log(qmu2/(qmass*qmass))        
*mb      endif
      dhqC1B2G = c1b2g/x
 
      return
      end
 
C     ==============================================
      double precision function dhqC12Q(chi,qmu2,nf)
C     ==============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c12q   = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi   = 1./factor
        eta  = xi * (1.-x)/(4.*x) - 1.
        bet  = sqrt(eta/(1.+eta))
        bet3 = bet*bet*bet
        rho  = 1./(1.+eta)
C--     Eq. (26) in hep-ph/9411431         
        c12q = cftf * (h1_hlq(eta,xi)+h1_htq(eta,xi)) +
     +  cftf * bet3 * (gfun_l(eta,xi)+gfun_t(eta,xi))
C--     This takes care of the factor in front of the integral     
        c12q = c12q*4.*pi/factor
*mb      endif
      dhqC12Q = c12q/x
 
      return
      end
 
C     ===============================================
      double precision function dhqC1B2Q(chi,qmu2,nf)
C     ===============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      c1b2q  = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi    = 1./factor
        eta   = xi * (1.-x)/(4.*x) - 1.
        bet   = sqrt(eta/(1.+eta))
        bet3  = bet*bet*bet
        rho   = 1./(1.+eta)
C--     Eq. (27) in hep-ph/9411431        
        c1b2q = cftf * (h1bar_hlq(eta,xi)+h1bar_htq(eta,xi)) +
     +   cftf * bet3 * (gbar_l(eta,xi)+gbar_t(eta,xi))
C--     This takes care of the factor in front of the integral     
        c1b2q = c1b2q*4.*pi/factor
C--     And this of the factor ln(mu2/m2)
        c1b2q = c1b2q*log(qmu2/(qmass*qmass))        
*mb      endif
      dhqC1B2Q = c1b2q/x
 
      return
      end
 
C     ==============================================
      double precision function dhqD12Q(chi,qmu2,nf)
C     ==============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale
 
      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      d12q   = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi   = 1./factor
        eta  = xi * (1.-x)/(4.*x) - 1.
C--     Eq. (28) in hep-ph/9411431        
        if(qsq.le.1.5) then
          d12q = cftf * (h1f_llq(eta,xi)+h1f_ltq(eta,xi))
        else
          d12q = cftf * (h1_llq(eta,xi)+h1_ltq(eta,xi))
        endif
C--     This takes care of the factor in front of the integral        
        d12q = d12q*4.*pi/factor
*mb      endif
      dhqD12Q = d12q/x
 
      return
      end
 
C     ===============================================
      double precision function dhqD1B2Q(chi,qmu2,nf)
C     ===============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'ffcons.inc'
      
      common /hqpass/ qmass, ascale, bscale

      jf     = nf                                 !avoid compiler warning
      qsq    = max((ascale*qmu2 + bscale),0.25D0) !force Q2 > 0.25
      factor = qmass*qmass/qsq
      afac   = 1.D0 + 4.D0*factor
      x      = chi/afac
      d1b2q  = 0.
*mb      if(x.lt.(1./(1.+4.*factor))) then
        xi   = 1./factor
        eta  = xi * (1.-x)/(4.*x) - 1.
C--     Eq. (29) in hep-ph/9411431        
        if(qsq.le.1.5) then
          d1b2q = cftf * h1bar_ltq(eta,xi)
        else
          d1b2q = 0.
        endif
C--     This takes care of the factor in front of the integral        
        d1b2q = d1b2q*4.*pi/factor
C--     And this of the factor ln(mu2/m2)
        d1b2q = d1b2q*log(qmu2/(qmass*qmass))        
*mb      endif
      dhqD1B2Q = d1b2q/x
 
      return
      end
      
