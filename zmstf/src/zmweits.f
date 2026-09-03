
C--   file zmweits.f containing the zmstf weight routines

C--   subroutine zmfillw(nwords)
C--   double precision function beta0(iq,nf)
C--   double precision function tbet0(iq,nf)
c--   double precision function beta1(iq,nf)
C--   subroutine zmdumpwCPP(lun,file,ls)
C--   subroutine zmdumpw(lun,file)
C--   subroutine zmreadwCPP(lun,file,ls,nwords,ierr)
C--   subroutine zmreadw(lun,file,nwords,ierr)
C--   subroutine zmwfileCPP(file,ls)
C--   subroutine zmwfile(file)
C--   subroutine zmwtids

C-----------------------------------------------------------------------
CXXHDR    void zmfillw(int &nused);
C-----------------------------------------------------------------------
CXXHFW  #define fzmfillw FC_FUNC(zmfillw,ZMFILLW)
CXXHFW    void fzmfillw(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zmfillw(int &nused)
CXXWRP  {
CXXWRP    fzmfillw(&nused);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==========================
      subroutine zmfillw(nwords)
C     ==========================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'
      include 'zmwidee.inc'

      external dzmConst1,dzmAchi
      external dzmC2G,dzmC2Q,dzmCLG,dzmCLQ,dzmD3Q
      external dzmC2NN2A,dzmC2NS2B,dzmC2NN2C,dzmC2NC2A,dzmC2NC2C
      external dzmC2S2A ,dzmC2G2A ,dzmC2G2C ,dzmCLNN2A,dzmCLNN2C
      external dzmCLNC2A,dzmCLNC2C,dzmCLS2A ,dzmCLG2A ,dzmC3NP2A
      external dzmC3NS2B,dzmC3NP2C,dzmC3NM2A,dzmC3NM2C
      external dzmCLG3A ,dzmCLS3A ,dzmCLNP3A,dzmCLNP3C
      
      external beta0, tbet0, beta1

      dimension itypes(6)
      data itypes /-4, -38, 0, 0, 0, 0/

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('ZMFILLW')

C--   Initialize
      if(izini.ne.12345) then
        izini  = 12345
        izpdf  = 1
        ascale = 1.D0
        bscale = 0.D0
        jscale = 0
      endif
      
C--   Get QCDNUM logical unit number for messages
      call getint('lunq',lun)      
      
C--   Check if splitting functions weight tables are available      
      if(idSpfun('PQQ',1,1).eq.-1) then
        write(lun,
     &   '(/'' ZMFILLW: no spltting function weights available''/
     &   ''          please call FILLWT or READWT before ZMFILLW'')')
       stop
      endif  
      
C--   Book tables in store zmstor of nzmstor words
      npar = 0
      new  = 1     !only one table set, thank you
      call MakeTab(zmstor,nzmstor,itypes,npar,new,jsetw,nwords)
C--   Store in common block zmstore.inc
      nzused = nwords      
C--   Setup the table identifiers in common block /zmwidee/
      call zmwtids
      
      write(lun,'(/'' ZMFILLW: start weight calculations'',4I4)') 
     +                (abs(itypes(i)),i=1,4)

C--   Global id
      ig0 = 1000*jsetw

C--   LO: 
      call MakeWtD(zmstor,ig0+idwtLO,dzmConst1,dzmAchi)

C--   NLO:
      call makeWtA(zmstor,ig0+idC2G1,dzmC2G,dzmAchi)
      call makeWtB(zmstor,ig0+idC2Q1,dzmC2Q,dzmAchi,0)

      call makeWtA(zmstor,ig0+idCLG1,dzmCLG,dzmAchi)
      call makeWtA(zmstor,ig0+idCLQ1,dzmCLQ,dzmAchi)

      call copyWgt(zmstor,ig0+idC2Q1,ig0+idC3Q1,0)
      call makeWtA(zmstor,ig0+idC3Q1,dzmD3Q,dzmAchi)

C--   NNLO
      call makeWtA(zmstor,ig0+idC2P2,dzmC2NN2A,dzmAchi)
      call makeWtB(zmstor,ig0+idC2P2,dzmC2NS2B,dzmAchi,1)
      call makeWtD(zmstor,ig0+idC2P2,dzmC2NN2C,dzmAchi)

      call makeWtA(zmstor,ig0+idC2M2,dzmC2NC2A,dzmAchi)
      call makeWtB(zmstor,ig0+idC2M2,dzmC2NS2B,dzmAchi,1)
      call makeWtD(zmstor,ig0+idC2M2,dzmC2NC2C,dzmAchi)

      call copyWgt(zmstor,ig0+idC2P2,ig0+idC2S2,0)
      call makeWtA(zmstor,ig0+idC2S2,dzmC2S2A,dzmAchi)

      call makeWtA(zmstor,ig0+idC2G2,dzmC2G2A,dzmAchi)
      call makeWtD(zmstor,ig0+idC2G2,dzmC2G2C,dzmAchi)

      call makeWtA(zmstor,ig0+idCLP2,dzmCLNN2A,dzmAchi)
      call makeWtD(zmstor,ig0+idCLP2,dzmCLNN2C,dzmAchi)

      call makeWtA(zmstor,ig0+idCLM2,dzmCLNC2A,dzmAchi)
      call makeWtD(zmstor,ig0+idCLM2,dzmCLNC2C,dzmAchi)

      call copyWgt(zmstor,ig0+idCLP2,ig0+idCLS2,0)
      call makeWtA(zmstor,ig0+idCLS2,dzmCLS2A,dzmAchi)

      call makeWtA(zmstor,ig0+idCLG2,dzmCLG2A,dzmAchi)

      call makeWtA(zmstor,ig0+idC3P2,dzmC3NP2A,dzmAchi)
      call makeWtB(zmstor,ig0+idC3P2,dzmC3NS2B,dzmAchi,1)
      call makeWtD(zmstor,ig0+idC3P2,dzmC3NP2C,dzmAchi)

      call makeWtA(zmstor,ig0+idC3M2,dzmC3NM2A,dzmAchi)
      call makeWtB(zmstor,ig0+idC3M2,dzmC3NS2B,dzmAchi,1)
      call makeWtD(zmstor,ig0+idC3M2,dzmC3NM2C,dzmAchi)
      
C--   NNLO FL'
      call makeWtA(zmstor,ig0+idCLG3,dzmCLG3A,dzmAchi)

      call makeWtA(zmstor,ig0+idCLN3,dzmCLNP3A,dzmAchi)
      call makeWtD(zmstor,ig0+idCLN3,dzmCLNP3C,dzmAchi)
      
      call copyWgt(zmstor,ig0+idCLN3,ig0+idCLS3,0)
      call makeWtA(zmstor,ig0+idCLS3,dzmCLS3A,dzmAchi)
      
C--   Now go for the scale dependent pieces, all type-2
C--   F2
C--   C_{2,x}^(1,1)
      iC2s11 = ig0+216
      call CopyWgt(zmstor,idSpfun('PQQ',1,1),iC2s11,0)
      iC2g11 = ig0+217
      call CopyWgt(zmstor,idSpfun('PQG',1,1),iC2g11,0)
      iC2p11 = ig0+216
      iC2m11 = ig0+216
C--   C_{2,s}^(2,1)      
      iC2s21 = ig0+218
      call CopyWgt(zmstor,idSpfun('PQQ',2,1),iC2s21,0)
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idC2Q1,iC2s21,1)
      call WcrossW(zmstor,idSpfun('PGQ',1,1),ig0+idC2G1,iC2s21,1)
      call WtimesF(zmstor,beta0,ig0+idC2Q1,iC2s21,-1)
C--   C_{2,g}^(2,1)
      iC2g21 = ig0+219
      call CopyWgt(zmstor,idSpfun('PQG',2,1),iC2g21,0)
      call WcrossW(zmstor,idSpfun('PQG',1,1),ig0+idC2Q1,iC2g21,1)
      call WcrossW(zmstor,idSpfun('PGG',1,1),ig0+idC2G1,iC2g21,1)
      call WtimesF(zmstor,beta0,ig0+idC2G1,iC2g21,-1)
C--   C_{2,+}^(2,1)
      iC2p21 = ig0+220
      call CopyWgt(zmstor,idSpfun('PPL',2,1),iC2p21,0)
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idC2Q1,iC2p21,1)
      call WtimesF(zmstor,beta0,ig0+idC2Q1,iC2p21,-1)
C--   C_{2,-}^(2,1)
      iC2m21 = ig0+221
      call CopyWgt(zmstor,idSpfun('PMI',2,1),iC2m21,0)
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idC2Q1,iC2m21,1)
      call WtimesF(zmstor,beta0,ig0+idC2Q1,iC2m21,-1)
C--   C_{2,s}^(2,2)
      iC2s22 = ig0+222
      call WcrossW(zmstor,idSpfun('PQQ',1,1),iC2s11,iC2s22,0)
      call WcrossW(zmstor,idSpfun('PGQ',1,1),iC2g11,iC2s22,1)
      call WtimesF(zmstor,beta0,iC2s11,iC2s22,-1)
      call ScaleWt(zmstor,0.5D0,iC2s22)
C--   C_{2,g}^(2,2)
      iC2g22 = ig0+223
      call WcrossW(zmstor,idSpfun('PQG',1,1),iC2s11,iC2g22,0)
      call WcrossW(zmstor,idSpfun('PGG',1,1),iC2g11,iC2g22,1)
      call WtimesF(zmstor,beta0,iC2g11,iC2g22,-1)
      call ScaleWt(zmstor,0.5D0,iC2g22)
C--   C_{2,+}^(2,2)
      iC2p22 = ig0+224
      call WcrossW(zmstor,idSpfun('PQQ',1,1),iC2p11,iC2p22,0)
      call WtimesF(zmstor,beta0,iC2p11,iC2p22,-1)
      call ScaleWt(zmstor,0.5D0,iC2p22)
C--   C_{2,-}^(2,2)
      iC2m22 = ig0+224
      
C--   FL
C--   C_{L,s}^(2,1)
      iCLs21 = ig0+225
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idCLQ1,iCLs21,0)
      call WcrossW(zmstor,idSpfun('PGQ',1,1),ig0+idCLG1,iCLs21,1)
      call WtimesF(zmstor,beta0,ig0+idCLQ1,iCLs21,-1)
C--   C_{L,g}^(2,1)
      iCLg21 = ig0+226
      call WcrossW(zmstor,idSpfun('PQG',1,1),ig0+idCLQ1,iCLg21,0)
      call WcrossW(zmstor,idSpfun('PGG',1,1),ig0+idCLG1,iCLg21,1)
      call WtimesF(zmstor,beta0,ig0+idCLG1,iCLg21,-1)
C--   C_{L,+}^(2,1)
      iCLp21 = ig0+227
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idCLQ1,iCLp21,0)
      call WtimesF(zmstor,beta0,ig0+idCLQ1,iCLp21,-1) 
C--   C_{L,-}^(2,1)
      iCLm21 = ig0+227

C--   xF3
C--   C_{3,x}^(1,1)
      iC3p11 = ig0+216
      iC3m11 = ig0+216
C--   C_{3,+}^(2,1)
      iC3p21 = ig0+228
      call CopyWgt(zmstor,idSpfun('PPL',2,1),iC3p21,0)
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idC3Q1,iC3p21,1)
      call WtimesF(zmstor,beta0,ig0+idC3Q1,iC3p21,-1)
C--   C_{3,-}^(2,1)
      iC3m21 = ig0+229
      call CopyWgt(zmstor,idSpfun('PMI',2,1),iC3m21,0)
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idC3Q1,iC3m21,1)
      call WtimesF(zmstor,beta0,ig0+idC3Q1,iC3m21,-1)
C--   C_{3,+}^(2,2)
      iC3p22 = ig0+230
      call WcrossW(zmstor,idSpfun('PQQ',1,1),iC3p11,iC3p22,0)
      call WtimesF(zmstor,beta0,iC3p11,iC3p22,-1)
      call ScaleWt(zmstor,0.5D0,iC3p22)
C--   C_{3,-}^(2,2)
      iC3m22 = ig0+230
      
C--   FL'
C--   C_{L,s}^(3,1)
      iCLs31 = ig0+231
      call WcrossW(zmstor,idSpfun('PQQ',2,1),ig0+idCLQ1,iCLs31,0)
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idCLS2,iCLs31,1)
      call WcrossW(zmstor,idSpfun('PGQ',2,1),ig0+idCLG1,iCLs31,1)
      call WcrossW(zmstor,idSpfun('PGQ',1,1),ig0+idCLG2,iCLs31,1)
      call WtimesF(zmstor,beta1,ig0+idCLQ1,iCLs31,-1)
      call WtimesF(zmstor,tbet0,ig0+idCLS2,iCLs31,-1)
C--   C_{L,g}^(3,1)
      iCLg31 = ig0+232
      call WcrossW(zmstor,idSpfun('PQG',2,1),ig0+idCLQ1,iCLg31,0)
      call WcrossW(zmstor,idSpfun('PQG',1,1),ig0+idCLS2,iCLg31,1)
      call WcrossW(zmstor,idSpfun('PGG',2,1),ig0+idCLG1,iCLg31,1)
      call WcrossW(zmstor,idSpfun('PGG',1,1),ig0+idCLG2,iCLg31,1)
      call WtimesF(zmstor,beta1,ig0+idCLG1,iCLg31,-1)
      call WtimesF(zmstor,tbet0,ig0+idCLG2,iCLg31,-1)
C--   C_{L,+}^(3,1)
      iCLp31 = ig0+233
      call WcrossW(zmstor,idSpfun('PPL',2,1),ig0+idCLQ1,iCLp31,0)
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idCLP2,iCLp31,1)
      call WtimesF(zmstor,beta1,ig0+idCLQ1,iCLp31,-1)
      call WtimesF(zmstor,tbet0,ig0+idCLP2,iCLp31,-1)
C--   C_{L,-}^(3,1)
      iCLm31 = ig0+234
      call WcrossW(zmstor,idSpfun('PMI',2,1),ig0+idCLQ1,iCLm31,0)
      call WcrossW(zmstor,idSpfun('PQQ',1,1),ig0+idCLM2,iCLm31,1)
      call WtimesF(zmstor,beta1,ig0+idCLQ1,iCLm31,-1)
      call WtimesF(zmstor,tbet0,ig0+idCLM2,iCLm31,-1)
C--   C_{L,s}^(3,2)
      iCLs32 = ig0+235
      call WcrossW(zmstor,idSpfun('PQQ',1,1),iCLs21,iCLs32,0)
      call WcrossW(zmstor,idSpfun('PGQ',1,1),iCLg21,iCLs32,1)
      call WtimesF(zmstor,tbet0,iCLs21,iCLs32,-1)
      call ScaleWt(zmstor,0.5D0,iCLs32)
C--   C_{L,g}^(3,2)
      iCLg32 = ig0+236
      call WcrossW(zmstor,idSpfun('PQG',1,1),iCLs21,iCLg32,0)
      call WcrossW(zmstor,idSpfun('PGG',1,1),iCLg21,iCLg32,1)
      call WtimesF(zmstor,tbet0,iCLg21,iCLg32,-1)
      call ScaleWt(zmstor,0.5D0,iCLg32)
C--   C_{L,+}^(3,2)
      iCLp32 = ig0+237
      call WcrossW(zmstor,idSpfun('PQQ',1,1),iCLp21,iCLp32,0)
      call WtimesF(zmstor,tbet0,iCLp21,iCLp32,-1)
      call ScaleWt(zmstor,0.5D0,iCLp32)
C--   C_{P,-}^(3,2)
      iCLm32 = ig0+238
      call WcrossW(zmstor,idSpfun('PQQ',1,1),iCLm21,iCLm32,0)
      call WtimesF(zmstor,tbet0,iCLm21,iCLm32,-1)
      call ScaleWt(zmstor,0.5D0,iCLm32)

      write(lun,'('' ZMFILLW: calculations completed'')')       

C--   Clear subroutine name
      call clrUmsg     

      return
      end
      
C     ======================================      
      double precision function beta0(iq,nf)
C     ======================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)
      
      idum  = iq !avoid compiler warning
      beta0 = 11.D0/2.D0 - nf/3.D0
      
      return
      end 
      
C     ======================================      
      double precision function tbet0(iq,nf)
C     ======================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)
      
      idum  = iq !avoid compiler warning
      tbet0 = 2.D0 * beta0(iq,nf)
      
      return
      end                   
      
C     ======================================      
      double precision function beta1(iq,nf)
C     ======================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)
      
      idum  = iq !avoid compiler warning
      beta1 = 51.D0/2.D0 - 19.D0*nf/6.D0
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void zmdumpw(int lun, string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fzmdumpwcpp FC_FUNC(zmdumpwcpp,ZMDUMPWCPP)
CXXHFW    void fzmdumpwcpp(int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zmdumpw(int lun, string fname)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fzmdumpwcpp(&lun,cfname,&ls);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ==================================
      subroutine zmdumpwCPP(lun,file,ls)
C     ==================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      character*(100) file

      if(ls.gt.100) stop 'zmdumpwCPP: file name size > 100 characters'

      call zmdumpw(lun,file(1:ls))

      return
      end


C     ============================
      subroutine zmdumpw(lun,file)
C     ============================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'

      character*(*) file

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('ZMDUMPW')
C--   Write (chvers is in zmstf.inc)
      call DumpTab(zmstor,1,lun,file,chvers)

C--   Clear subroutine name
      call clrUmsg     

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void zmreadw(int lun, string fname, int &nused, int &ierr);
C-----------------------------------------------------------------------
CXXHFW  #define fzmreadwcpp FC_FUNC(zmreadwcpp,ZMREADWCPP)
CXXHFW    void fzmreadwcpp(int*, char*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zmreadw(int lun, string fname, int &nused, int &ierr)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fzmreadwcpp(&lun,cfname,&ls,&nused,&ierr);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================================
      subroutine zmreadwCPP(lun,file,ls,nwords,ierr)
C     ==============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      character*(100) file

      if(ls.gt.100) stop 'zmreadwCPP: file name size > 100 characters'

      call zmreadw(lun,file(1:ls),nwords,ierr)

      return
      end

C     ========================================
      subroutine zmreadw(lun,file,nwords,ierr)
C     ========================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'

      character*(*) file

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('ZMREADW')
C--   Setup the table identifiers
      call zmwtids
C--   Read (chvers is in zmstf.inc)
      call ReadTab(zmstor,nzmstor,lun,file,chvers,1,jset,nwords,ierr)
      if(ierr.ne.0) return
C--   Store in common block zmstore.inc
      nzused = nwords      
      
C--   Set ZMSTF initialized
      if(izini.ne.12345) then
        izini  = 12345
        izpdf  = 1
        ascale = 1.D0
        bscale = 0.D0
        jscale = 0
      endif

C--   Store jset
      jsetw = jset

C--   Clear subroutine name
      call clrUmsg     

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void zmwfile(string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fzmwfilecpp FC_FUNC(zmwfilecpp,ZMWFILECPP)
CXXHFW    void fzmwfilecpp(char*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void zmwfile(string fname)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fzmwfilecpp(cfname,&ls);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================
      subroutine zmwfileCPP(file,ls)
C     ==============================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      character*(100) file

      if(ls.gt.100) stop 'zmwfileCPP: file name size > 100 characters'

      call zmwfile(file(1:ls))

      return
      end

C     ========================
      subroutine zmwfile(file)
C     ========================

C--   Maintain a weight file on disk

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'zmstf.inc'
      include 'zmstore.inc'
      include 'zmscale.inc'

      character*(*) file

C--   Set subroutine name (must be cleared before exit)
      call setUmsg('ZMWFILE')

C--   Get logical unit number
      lun = nxtlun(0)
C--   Read file
      call zmreadw(lun,file,nwords,ierr)
C--   Oeps! eror
      if(ierr.ne.0) then
        call zmfillw(nwords)
        call zmdumpw(lun,file)
      endif

C--   Clear subroutine name
      call clrUmsg

      return
      end

C     ==================
      subroutine zmwtids
C     ==================

C--   Setup the weight table identifiers
C--   This routine is called by zmfillw and zmreadw

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)

      include 'zmstf.inc'
      include 'zmwidee.inc'

      idwtLO = 101   !LO     delta(1-x)
      idC2G1 = 201   !NLO    F2G
      idC2Q1 = 102   !NLO    F2Q
      idCLG1 = 202   !NLO    FLG
      idCLQ1 = 103   !NLO    FLQ
      idC3Q1 = 104   !NLO    xF3
      idC2P2 = 203   !NNLO   F2ns+
      idC2M2 = 204   !NNLO   F2ns-
      idC2S2 = 205   !NNLO   F2si 
      idC2G2 = 206   !NNLO   F2gl
      idCLP2 = 207   !NNLO   FLns+
      idCLM2 = 208   !NNLO   FLns- 
      idCLS2 = 209   !NNLO   FLsi 
      idCLG2 = 210   !NNLO   FLgl
      idC3P2 = 211   !NNLO   F3ns+
      idC3M2 = 212   !NNLO   F3ns-
      idCLG3 = 213   !NNLO   FL'gl
      idCLS3 = 214   !NNLO   FL'si
      idCLN3 = 215   !NNLO   FL'ns

      return
      end
