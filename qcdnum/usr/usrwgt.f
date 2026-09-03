
C--   This is the file usrwgt.f containing the qcdnum weight routines

C--   subroutine fillwt(iselect,idmin,idmax,nwlast)
C--   subroutine fillwc(usub,idmin,idmax,nwlast)
C--   subroutine dmpwgtCPP(jset,lun,file,ls)
C--   subroutine dmpwgt(jset,lun,file)
C--   subroutine readwtCPP(lun,file,ls,idmin,idmax,nwlast,ierr)
C--   subroutine readwt(lun,file,ityp,idum,nwlast,ierr)
C--   subroutine wtfileCPP(itype,file,ls)
C--   subroutine wtfile(ityp,file)
C--   subroutine nwused(nwtot,nwuse,nwtab)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM weight routines from usrwgt.f                    */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM weight routines from usrwgt.f                      */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM weight routines from usrwgt.f                      */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C==   ===============================================================
C==   Weight calculation, dump and read =============================
C==   ===============================================================

C-----------------------------------------------------------------------
CXXHDR    void fillwt(int itype, int &idmin, int &idmax, int &nwds);
C-----------------------------------------------------------------------
CXXHFW  #define ffillwt FC_FUNC(fillwt,FILLWT)
CXXHFW    void ffillwt(int*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void fillwt(int itype, int &idmin, int &idmax, int &nwds)
CXXWRP  {
CXXWRP    ffillwt(&itype,&idmin,&idmax,&nwds);
CXXWRP  }
C-----------------------------------------------------------------------
            
C     =============================================
      subroutine fillwt(iselect,idmin,idmax,nwlast)
C     =============================================

C--   Fill standard weight tables
C--
C--   iselect (in)     1 = Unpolarised 
C--                    2 = Polarised
C--                    3 = Fragmentation functions (timelike)
C--                 else = Unpolarised
C--   idmin   (out)        first pdf identifier (always 0 = gluon)
C--   idmax   (out)        last  pdf identifier
C--   nwlast  (out)        last word used in the store < 0 no space

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)
      
      external sqcFilWU, sqcFilWP, sqcFilWF

      character*80 subnam
      data subnam /'FILLWT ( ITYPE, IDMIN, IDMAX, NWDS )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Yes/no re-initialize, thats the question
      if(.not.Lwtini7) call sqcIniWt

      if(iselect.eq.2) then
      
C--     Calculate polarised weight tables
        write(lunerr1,'(/
     +        '' FILLWT: start polarised weight calculations'')')
        write(lunerr1,'( '' Subgrids'',I5, 
     +        '' Subgrid points'',100I5)') nyg2,(nyy2(i),i=1,nyg2)
        call sqcFilWt(sqcFilWP,lunwgt1,2,nwlast,ierr)

        jselect = 2
        
      elseif(iselect.eq.3) then
      
C--     Calculate timelike weight tables (fragmentation functions)
        write(lunerr1,'(/
     +        '' FILLWT: start fragmentation weight calculations'')')
        write(lunerr1,'( '' Subgrids'',I5, 
     +        '' Subgrid points'',100I5)') nyg2,(nyy2(i),i=1,nyg2)
        call sqcFilWt(sqcFilWF,lunwgt1,3,nwlast,ierr)
        jselect = 3
      
      else
            
C--     Calculate unpolarised weight tables (default)
        write(lunerr1,'(/
     +        '' FILLWT: start unpolarised weight calculations'')')
        write(lunerr1,'( '' Subgrids'',I5, 
     +        '' Subgrid points'',100I5)') nyg2,(nyy2(i),i=1,nyg2)
        call sqcFilWt(sqcFilWU,lunwgt1,1,nwlast,ierr)
        jselect = 1
      
      endif
      
C--   Weight tables already existed so inform that nothing has been done
      if(ierr.eq.-1) then
        write(lunerr1,'( 
     +           '' Tables already exist --> nothing done'')')
      endif       
        
C--   Make sure that there is enough space to hold first word of next set
      nwfirst = abs(nwlast)+1

C--   Not enough space or iset count exceeded
      if(ierr.eq.-2 .or. ierr.eq.-3) call sqcMemMsg(subnam,nwfirst,ierr)

      write(lunerr1,'('' FILLWT: weight calculations completed''/)')
      
C--   Idmin and Idmax
      idmin  = 0
      idmax  = mnf0

C--   Update status bits
      call sqcSetflg(iset,idel,jselect)
      
      return
      end

C     ==========================================
      subroutine fillwc(usub,idmin,idmax,nwlast)
C     ==========================================

C--   Fill custom weight tables
C--
C--   usub    (in)         user defined subroutine     
C--   idmin   (out)        first pdf identifier (always 0 = gluon)
C--   idmax   (out)        last  pdf identifier
C--   nwlast  (out)        last word used in the store < 0 no space

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)
      
      character*60 emsg
      character*10 etxt
      
      external usub

      character*80 subnam
      data subnam /'FILLWC ( MYSUB, IDMIN, IDMAX, NWDS )'/

*mb   Fillwc is disabled for now
      call sqcErrMsg(subnam,
     +   'FILLWC disabled: use subroutine EVDGLAP for custom evolution')

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Yes/no re-initialize, thats the question
      if(.not.Lwtini7) call sqcIniWt

C--   Calculate custom weight tables
      write(lunerr1,'(/
     +      '' FILLWC: start custom weight calculations'')')
      write(lunerr1,'( '' Subgrids'',I5, 
     +      '' Subgrid points'',100I5)') nyg2,(nyy2(i),i=1,nyg2)
      call sqcFilWt(usub,lunwgt1,4,nwlast,ierr)
      
C--   Custom weight tables already exist so thats an error
      if(ierr.eq.-1) then
        call sqcErrMsg(subnam,'Custom tables already exist')
      endif
C--   Mxord not in range 1-3
      if(ierr.eq.-2) then
        call sqcErrMsg(subnam,'Maxord not in range [1-3]')
      endif             
        
C--   Make sure that there is enough space to hold first word of next set
      nwfirst = abs(nwlast)+1
      if(nwfirst.gt.nwf0) then
        call smb_itoch(nwfirst,etxt,ltxt)
        write(emsg,'(''Need at least '',A,
     &    '' words --> increase NWF0 '',
     &    ''in qcdnum.inc'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif
      write(lunerr1,'('' FILLWC: weight calculations completed''/)')
      
C--   Idmin and Idmax
      idmin  = 0
      idmax  = mnf0

C--   Update status bits
      call sqcSetflg(iset,idel,4)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void dmpwgt(int itype, int lun, string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fdmpwgtcpp FC_FUNC(dmpwgtcpp,DMPWGTCPP)
CXXHFW    void fdmpwgtcpp(int*, int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void dmpwgt(int itype, int lun, string fname)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fdmpwgtcpp(&itype,&lun,cfname,&ls);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ======================================
      subroutine dmpwgtCPP(jset,lun,file,ls)
C     ======================================

      implicit double precision (a-h,o-z)

      character*(100) file

      if(ls.gt.100) stop 'dmpwgtCPP: input file name > 100 characters'

      call dmpwgt(jset,lun,file(1:ls))

      return
      end

C     ================================
      subroutine dmpwgt(ityp,lun,file)
C     ================================

C--   Dump weights to a disk file (unformatted write)
C--
C--   ityp   (in)   1=unpol, 2=pol, 3=timelike, 4=custom

      implicit double precision (a-h,o-z)

      character*(*) file

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)
      
      character*5 etxt
      character*11 txt(4)
C--               12345678901
      data txt / 'unpolarised',
     +           'polarised  ',
     +           'time-like  ',
     +           'custom     '  /

      character*80 subnam
      data subnam /'DMPWGT ( ITYPE, LUN, FILE )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check logical unit number
      if(lun.le.0 .or. lun.eq.6) then
        call smb_itoch(lun,etxt,ltxt)
        call sqcErrMsg(subnam,
     +    'Invalid logical unit number lun = '//etxt(1:ltxt))      
      endif

C--   Check ityp
      call sqcIlele(subnam,'ISET',1,ityp,4,' ')

      leng = imb_lenoc(txt(ityp))
      write(lunerr1 ,'(/'' DMPWGT: dump '',A,'' weight tables'')')
     +  txt(ityp)(1:leng)

C--   Now dumpit 
      open(unit=lun,file=file,form='unformatted',
     +     status='unknown',err=500)
      call sqcDumpWt(lun,ityp,' ',ierr)
      close(lun)
      if(ierr.eq.1) goto 501
      if(ierr.eq.2) goto 502

      write(lunerr1 ,'(''         weights written to '',A/)') file
     
      return

  500 continue
C--   Open error
      call sqcErrMsg(subnam,'Cannot open output weight file')
      return
  501 continue
C--   Write error
      call sqcErrMsg(subnam,'Write error on output weight file')
      return
  502 continue
C--   No weight tables available
      leng = imb_lenoc(txt(ityp))
      call sqcErrMsg(subnam,
     +     'No '//txt(ityp)(1:leng)//' weight tables available')
      return      
  
      end

C-----------------------------------------------------------------------
CXXHDR    void readwt(int lun, string fname, int &idmin, int &idmax, int &nwds, int &ierr);
C-----------------------------------------------------------------------
CXXHFW  #define freadwtcpp FC_FUNC(readwtcpp,READWTCPP)
CXXHFW    void freadwtcpp(int*, char*, int*, int*, int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void readwt(int lun, string fname, int &idmin, int &idmax, int &nwds, int &ierr)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    freadwtcpp(&lun,cfname,&ls,&idmin,&idmax,&nwds,&ierr);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     =========================================================
      subroutine readwtCPP(lun,file,ls,idmin,idmax,nwlast,ierr)
C     =========================================================

      implicit double precision (a-h,o-z)

      character*(100) file

      if(ls.gt.100) stop 'readwtCPP: input file name > 100 characters'

      call readwt(lun,file(1:ls),idmin,idmax,nwlast,ierr)

      return
      end

C     ===================================================
      subroutine readwt(lun,file,idmin,idmax,nwlast,ierr)
C     ===================================================

C--   Read weights from a disk file (unformatted read)
C--
C--   lun     (in)   nput logical unit number
C--   file    (in)   input file name
C--   idmin   (out)  index of first pdf table (always 0 = gluon)
C--   idmax   (out)  index of last  pdf table
C--   nwlast  (out)  last word occupied in store < 0 not enough space 
C--   ierr    (out)  0 = all OK
C--                  1 = read error
C--                  2 = problem with QCDNUM version
C--                  3 = key mismatch
C--                  4 = x-mu2 grid not the same

      implicit double precision (a-h,o-z)

      character*(*) file

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)
      
      character*13 txt(4)
C--               1234567890123      
      data txt / 'unpolarised  ',
     +           'polarised    ',
     +           'fragmentation',
     +           'custom       '  /
      
      character*80 subnam
      data subnam /'READWT ( LUN, FNAM, ITYP, IDUM, NWDS, IERR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Yes/no re-initialize, thats the question
      if(.not.Lwtini7) call sqcIniWt
      
C--   Read weight tables
      write(lunerr1 ,'(/'' READWT: open file '',A)')
     +      file
      open(unit=lun,file=file,form='unformatted',status='old',
     +     err=500)
      call sqcReadWt(lun,' ',nwlast,iread,ierr)
      close(lun)

C--   First word of next set (next call to fillwt or readwt)
      nwnext = abs(nwlast)+1

C--   Not enough space
      if(nwnext.gt.nwf0) then
        call sqcMemMsg(subnam,nwnext,-2)
C--   Iset count exceeded
      elseif(ierr.eq.6) then
        call sqcMemMsg(subnam,nwnext,-3)
      endif

C--   Error reading table set
      if(ierr.ne.0) return

C--   Idmin and Idmax
      idmin  = 0
      idmax  = mnf0

C--   Update status bits
      if(iread.gt.0) then
        call sqcSetflg(iset,idel,iread)
        leng = imb_lenoc(txt(iread))
        write(lunerr1,'(''         read '', A,
     +     '' weight tables'')') txt(iread)(1:leng)
      elseif(iread.lt.0) then
          leng = imb_lenoc(txt(-iread))
          write(lunerr1,'(9X,A, '' tables already exist'', 
     +     '' --> nothing done'')') txt(-iread)(1:leng)
      else
        stop 'READWT : unknown weight type read in ---> STOP'
      endif
      write(lunerr1,'(/)')
      
      return

  500 continue
C--   Read error
      ierr = 1
      return

      end

C-----------------------------------------------------------------------
CXXHDR    void wtfile(int itype, string fname);
C-----------------------------------------------------------------------
CXXHFW  #define fwtfilecpp FC_FUNC(wtfilecpp,WTFILECPP)
CXXHFW    void fwtfilecpp(int*, char*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void wtfile(int itype, string fname)
CXXWRP  {
CXXWRP    int ls = fname.size();
CXXWRP    char *cfname = new char[ls+1];
CXXWRP    strcpy(cfname,fname.c_str());
CXXWRP    fwtfilecpp(&itype,cfname,&ls);
CXXWRP    delete[] cfname;
CXXWRP  }
C-----------------------------------------------------------------------

C     ===================================
      subroutine wtfileCPP(itype,file,ls)
C     ===================================

      implicit double precision (a-h,o-z)

      character*(100) file

      if(ls.gt.100) stop 'wtfileCPP: input file name > 100 characters'

      call wtfile(itype,file(1:ls))

      return
      end

C     =============================
      subroutine wtfile(itype,file)
C     =============================

C--   Maintain a weight file on disk

      implicit double precision (a-h,o-z)

      character*(*) file
      character*80 subnam
      data subnam /'WTFILE ( ITYPE, FILENAME )'/

      character*13 txt(4)
C--               1234567890123      
      data txt / 'unpolarised  ',
     +           'polarised    ',
     +           'fragmentation',
     +           'custom       '  /

C--   Check itype
      call sqcIlele(subnam,'ISET',1,itype,4,' ')

C--   Get logical unit number
      lun = nxtlun(0)

C--   Check itype is the same as on the file
      open(unit=lun,file=file,form='unformatted',status='old',err=500)
      call sqcReadWt(lun,'GIVETYPE',nwlast,jtype,ierr)
      close(unit=lun)
C--   Next call to readwt will take care of the errors
      if(ierr.ne.0) goto 500

      if(jtype.ne.itype) then
        leng = imb_lenoc(txt(itype))
        call sqcErrMsg(subnam,
     +  'File '//file//' does not contain '//txt(itype)(1:leng)//
     +  ' weight tables')
      endif

 500  continue
      call setUmsg('WTFILE')
      call readwt(lun,file,idmin,idmax,nw,ierr)
      if(ierr.ne.0) then
        call fillwt(itype,idmin,idmax,nw)
        call dmpwgt(itype,lun,file)
      endif
      call clrUmsg

      close(unit=lun)

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void nwused(int &nwtot, int &nwuse, int &ndummy);
C-----------------------------------------------------------------------
CXXHFW  #define fnwused FC_FUNC(nwused,NWUSED)
CXXHFW    void fnwused(int*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void nwused(int &nwtot, int &nwuse, int &ndummy)
CXXWRP  {
CXXWRP    fnwused(&nwtot,&nwuse,&ndummy);
CXXWRP  }
C-----------------------------------------------------------------------

C     ====================================      
      subroutine nwused(nwtot,nwuse,nwtab)      
C     ====================================

C--   Get store parameters

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      character*80 subnam
      data subnam /'NWUSED ( NWTOT, NWUSE, NWTAB )'/

      logical lqcIdExists

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

      nwtot  = nwf0
      nwuse  = iqcGetNumberOfWords(stor7,1)
      do jset = 1,mst0
        id = 1000*jset+501
        if(lqcIdExists(stor7,id)) nwtab = iqcGetTabLeng(stor7,id,2)
      enddo
      
      return
      end
      
