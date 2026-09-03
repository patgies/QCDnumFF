
C--   file usrstore.f containing local workspace managment routines

C--   subroutine MakeTab(w,nw,jtypes,new,jset,nwords)
C--   subroutine SetParW(w,jset,par,n)
C--   subroutine GetParW(w,jset,par,n)
C--   integer function iupdate(id,iset)
C--   subroutine dumptabCPP(w,jset,lun,fnam,lf,fkey,lk)
C--   subroutine DumpTab(w,jset,lun,file,key)
C--   subroutine readtabCPP(w,nw,lun,fnam,lf,fkey,lk,new,jset,nwd,ierr)
C--   subroutine ReadTab(w,nw,lun,file,key,new,jset,nwords,ierr)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  QCDNUM workspace routines from usrstore.f               */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  QCDNUM workspace routines from usrstore.f                 */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  QCDNUM workspace routines from usrstore.f                 */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    void maketab(double *w, int nw, int *jtypes, int npar,
CXXHDR                 int newt, int &jset, int &nwords);
C-----------------------------------------------------------------------
CXXHFW  #define fmaketab FC_FUNC(maketab,MAKETAB)
CXXHFW    void fmaketab(double*,int*,int*,int*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void maketab(double *w, int nw, int *jtypes, int npar,
CXXWRP               int newt, int &jset, int &nwords)
CXXWRP  {
CXXWRP    fmaketab(w, &nw, jtypes, &npar, &newt, &jset, &nwords);
CXXWRP  }
C-----------------------------------------------------------------------

C     ====================================================
      subroutine MakeTab(w,nw,jtypes,npar,new,jset,nwords)
C     ====================================================

C--   Partition the store into tables
C--
C--   w       (in) array dimensioned nw in the calling routine
C--   nw      (in) number of words in w
C--   jtypes  (in) jtypes(i) = requested # of type-i tables  i = 1,...,6
C--   npar    (in) space reserved for user parameters
C--   new     (in) 1=overwrite, otherwise add new set
C--   jset   (out) table set identifier, assigned by QCDNUM
C--                -1 attempt to book empty set
C--                -2 not enough space
C--                -3 iset count exceeded
C--   nwords (out) total number of words needed; < 0 not enough space

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension w(*), jtypes(*), itypes(mtyp0)

      character*80 subnam
      data subnam /'MAKETAB ( W, NW, ITYPES, NPAR, NEW, ISET, NWDS )'/

      character*80 emsg
      character*10 etxt

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Copy input array to avoid changing input
      do i = 1,6
        itypes(i) = jtypes(i)
      enddo
C--   Add subgrid pointer tables if there are pdf tables
      if(jtypes(5).ne.0) then
        itypes(7) = nsubt0
      else
        itypes(7) = 0
      endif

C--   Check user input
      call sqcIlele(subnam,'ITYPES(1)',0,abs(itypes(1)),99,' ')
      call sqcIlele(subnam,'ITYPES(2)',0,abs(itypes(2)),99,' ')
      call sqcIlele(subnam,'ITYPES(3)',0,abs(itypes(3)),99,' ')
      call sqcIlele(subnam,'ITYPES(4)',0,abs(itypes(4)),99,' ')
      call sqcIlele(subnam,'ITYPES(5)',0,abs(itypes(5)),99,' ')
      call sqcIlele(subnam,'ITYPES(6)',0,abs(itypes(6)),99,' ')
C--   Check new
      call sqcIlele(subnam,'NEW',0,new,1,' ')
C--   Do the work (mpar0 is the max number of evolution parameters)
      call sqcMakeTab(w,nw,itypes,mpar0,npar,new,jset,nwords)
C--   Check for errors
      if(jset.eq.-1) then
        call sqcErrMsg(subnam,'No tables to book')
      elseif(jset.eq.-2) then
        call smb_itoch(abs(nwords)+1,etxt,ltxt)
        write(emsg,'(''Increase NW to at least '',A,
     +               '' words'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      elseif(jset.eq.-3) then
        call smb_itoch(mst0,etxt,ltxt)
        write(emsg,'(''Setcount '',A,'' exceeded --> increase MST0 '',
     +               ''in qcdnum.inc and recompile'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void setparw(double *w, int jset, double *par, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fsetparw FC_FUNC(setparw,SETPARW)
CXXHFW    void fsetparw(double*,int*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void setparw(double *w, int jset, double *par, int n)
CXXWRP  {
CXXWRP    fsetparw(w, &jset, par, &n);
CXXWRP  }
C-----------------------------------------------------------------------

C     ================================
      subroutine SetParW(w,jset,par,n)
C     ================================

C--   Store parameters in table set jset
C--
C--   w        (in) : store dimensioned in the calling routine
C--   jset     (in) : input table set
C--   par(n)   (in) : input set of parameters
C--   n        (in) : number of parameters to be stored

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical lqcIsetExists

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension w(*),par(*)

      character*80 subnam
      data subnam /'SETPARW ( W, ISET, PAR, N )'/

      character*80 emsg
      character*10 etxt

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check store partitioned and jset exists
      if(.not.lqcIsetExists(w,jset)) then
        call smb_itoch(jset,etxt,ltxt)
        write(emsg,'(''W not partitioned or ISET = '',A,
     +               '' does not exist'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif

C--   Check enough words
      niw = iqcGetNumberOfUparam(w,jset)
      call sqcIlele(subnam,'N',1,n,niw,' ')

C--   Do the work
      ia0 = iqcFirstWordOfUparam(w,jset)-1
      do i = 1,n
        w(ia0+i) = par(i)
      enddo  

C--   Update status bits
      call sqcSetflg(iset,idel,0)     

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void getparw(double *w, int jset, double *par, int n);
C-----------------------------------------------------------------------
CXXHFW  #define fgetparw FC_FUNC(getparw,GETPARW)
CXXHFW    void fgetparw(double*,int*,double*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void getparw(double *w, int jset, double *par, int n)
CXXWRP  {
CXXWRP    fgetparw(w, &jset, par, &n);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ================================
      subroutine GetParW(w,jset,par,n)
C     ================================

C--   Read parameters from table set jset
C--
C--   w        (in) : store dimensioned in the calling routine
C--   jset     (in) : input table set
C--   par(n)  (out) : output set of parameters
C--   n        (in) : number of parameters to be retreived

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical lqcIsetExists

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension w(*),par(*)

      character*80 subnam
      data subnam /'GETPARW ( W, ISET, PAR, N )'/

      character*80 emsg
      character*10 etxt

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check store partitioned and jset exists
      if(.not.lqcIsetExists(w,jset)) then
        call smb_itoch(jset,etxt,ltxt)
        write(emsg,'(''W not partitioned or ISET = '',A,
     +               '' does not exist'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif

C--   Check enough words
      niw = iqcGetNumberOfUParam(w,jset)
      call sqcIlele(subnam,'N',1,n,niw,' ')

C--   Do the work
      ia0 = iqcFirstWordOfUparam(w,jset)-1
      do i = 1,n
        par(i) = w(ia0+i)
      enddo  

C--   Update status bits
      call sqcSetflg(iset,idel,0)     

      return
      end

C     =================================
      integer function iupdate(id,iset)
C     =================================

C--   Replace the iset field of a local workspace identifier

      implicit double precision (a-h,o-z)

      iupdate = id

C--   This is certainly not a local workspace identifier
      if(id.lt.1000) return

      idlocal = iqcGetLocalId(id)

C--   Check if id is a pdf in internal memory
      if(idlocal.ge.0 .and. idlocal.le.12) return

C--   Now replace the iset field
      iupdate = idlocal + 1000*iset

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void dumptab(double *w, int jset, int lun,
CXXHDR                 string fnam, string fkey);
C-----------------------------------------------------------------------
CXXHFW  #define fdumptabcpp FC_FUNC(dumptabcpp,DUMPTABCPP)
CXXHFW    void fdumptabcpp(double*,int*,int*,char*,int*,char*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void dumptab(double *w, int jset, int lun,
CXXWRP               string fnam, string fkey)
CXXWRP  {
CXXWRP    int lf = fnam.size();
CXXWRP    int lk = fkey.size();
CXXWRP    char *cfnam = new char[lf+1];
CXXWRP    char *cfkey = new char[lk+1];
CXXWRP    strcpy(cfnam,fnam.c_str());
CXXWRP    strcpy(cfkey,fkey.c_str());
CXXWRP    fdumptabcpp(w, &jset, &lun, cfnam, &lf, cfkey, &lk);
CXXWRP    delete[] cfnam;
CXXWRP    delete[] cfkey;
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================================
      subroutine dumptabCPP(w,jset,lun,fnam,lf,fkey,lk)
C     =================================================

      implicit double precision (a-h,o-z)

      dimension w(*)
      character*(100) fnam, fkey

      if(lf.gt.100) stop 'dumptabCPP: input file name > 100 characters'
      if(lk.gt.100) stop 'dumptabCPP: input key name > 100 characters'

      call DumpTab(w,jset,lun,fnam(1:lf),fkey(1:lk))

      return
      end
      
C     =======================================
      subroutine DumpTab(w,jset,lun,file,key)
C     =======================================

C--   Dump store to disk
C--
C--   w       (in)   store
C--   jset    (in)   table set
C--   lun     (in)   logical unit number
C--   file    (in)   output file name
C--   key     (in)   character string to be written on the file

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical lqcIsetExists

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension w(*)
      character*(*) file, key 

      character*80 subnam
      data subnam /'DUMPTAB ( W, ISET, LUN, FILE, KEY )'/

      character*80 emsg
      character*10 etxt

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check store partitioned and jset exists
      if(.not.lqcIsetExists(w,jset)) then
        call smb_itoch(jset,etxt,ltxt)
        write(emsg,'(''W not partitioned or ISET = '',A,
     +               '' does not exist'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif

C--   Open output file
      open(unit=lun,file=file,form='unformatted',status='unknown',
     +     err=500)
C--   Do the work
      call sqcDumpTab(w,jset,lun,key,ierr)
      close(lun)
      if(ierr.ne.0) goto 501

      write(lunerr1 ,'(/'' DUMPTAB: tables written to '',A/)')
     +      file 

C--   Update status bits
      call sqcSetflg(iset,idel,0)     

      return

  500 continue
C--   Open error
      call sqcErrMsg(subnam,'Cannot open output file')
      return
  501 continue
C--   Write error
      call sqcErrMsg(subnam,'Write error on output file')
      return

      end

C-----------------------------------------------------------------------
CXXHDR    void readtab(double *w, int nw, int lun, string fnam,
CXXHDR                 string fkey, int newt, int &jset,
CXXHDR                 int &nwords, int &ierr);
C-----------------------------------------------------------------------
CXXHFW  #define freadtabcpp FC_FUNC(readtabcpp,READTABCPP)
CXXHFW    void freadtabcpp(double*,int*,int*,char*,int*,char*,int*,
CXXHFW                     int*,int*,int*,int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void readtab(double *w, int nw, int lun, string fnam,
CXXWRP               string fkey, int newt, int &jset,
CXXWRP               int &nwords, int &ierr)
CXXWRP  {
CXXWRP    int lf = fnam.size();
CXXWRP    int lk = fkey.size();
CXXWRP    char *cfnam = new char[lf+1];
CXXWRP    char *cfkey = new char[lk+1];
CXXWRP    strcpy(cfnam,fnam.c_str());
CXXWRP    strcpy(cfkey,fkey.c_str());
CXXWRP    freadtabcpp(w, &nw, &lun, cfnam, &lf, cfkey, &lk,
CXXWRP                &newt, &jset, &nwords, &ierr);
CXXWRP    delete[] cfnam;
CXXWRP    delete[] cfkey;
CXXWRP  }
C-----------------------------------------------------------------------

C     =================================================================
      subroutine readtabCPP(w,nw,lun,fnam,lf,fkey,lk,new,jset,nwd,ierr)
C     =================================================================

      implicit double precision (a-h,o-z)

      dimension w(*)
      character*(100) fnam, fkey

      if(lf.gt.100) stop 'readtabCPP: input file name > 100 characters'
      if(lk.gt.100) stop 'readtabCPP: input key name > 100 characters'

      call ReadTab(w,nw,lun,fnam(1:lf),fkey(1:lk),new,jset,nwd,ierr)

      return
      end

C     ==========================================================
      subroutine ReadTab(w,nw,lun,file,key,new,jset,nwords,ierr)
C     ==========================================================

C--   Read store from disk
C--
C--   w       (in)  : store dimensioned in the calling routine
C--   nw      (in)  : dimension of w
C--   lun     (in)  : input logical unit number
C--   file    (in)  : input file name
C--   key     (in)  : key to match that on the file
C--   new     (in)  : 1=overwrite, otherwise add new set
C--   jset    (out) : table set number assigned by ReadTab
C--   nwords  (out) : number of words read in
C--   ierr    (out) :   0  : all OK
C--                     1  : read error
C--                     2  : problem with QCDNUM version
C--                     3  : problem with file stamp
C--                     4  : x-mu2 grid not the same 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,ichk,iset,idel)

      dimension w(*)
      character*(*) file, key

      character*80 subnam
      data subnam
     +       /'READTAB ( W, NW, LUN, FNAM, KEY, NEW, ISET, NWU, IERR )'/

      character*80 emsg
      character*10 etxt

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check new
      call sqcIlele(subnam,'NEW',0,new,1,' ')

C--   Open input file
      open(unit=lun,file=file,form='unformatted',status='old',
     +     err=500)

C--   Do the work
      call sqcReadTab(w,nw,lun,key,new,jset,nwords,ierr)
      close(lun)

C--   Check for errors
      if(ierr.eq.0) then
        write(lunerr1 ,'(/'' TABREAD: tables read in from '',A/)')
     +        file
      elseif(ierr.eq.5) then
        call smb_itoch(abs(nwords)+1,etxt,ltxt)
        write(emsg,'(''Increase NW to at least '',A,
     +               '' words'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      elseif(ierr.eq.6) then
        call smb_itoch(mst0,etxt,ltxt)
        write(emsg,'(''Setcount '',A,'' exceeded --> increase MST0 '',
     +               ''in qcdnum.inc and recompile'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)      

      return

  500 continue
C--   Read error
      ierr = 1
      return

      end

