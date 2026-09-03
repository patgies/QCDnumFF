
C--   This is the file hqstfs.f containing heavy quark stf routines
C--
C--   integer function ihqvers()
C--   subroutine hqwords(ntot,nused)  
C--   subroutine hqparms(qmas,aq,bq)
C--   double precision function hqqfrmu(qmu2)
C--   double precision function hqmufrq(q2)
C--   logical function LhqRvar(epsi)
C--   logical function LhqQvar(epsi)
C--
C--   subroutine hswitch(jset)
C--   
C--   subroutine hqStFun(istf,kcbt,def,x,q,f,n,ichk)
C--   double precision function as1fun(ix,iq,nf,ithrs)
C--   double precision function as2fun(ix,iq,nf,ithrs)
C--      
C--   subroutine hqSlowF(istf,icbt,def,x,q,f,n,ichk)
C--   double precision function dhqFij(ix,iq)
C--   double precision function dhqGetF(w,id0,icbt,def,ix,iq)

C-----------------------------------------------------------------------
CXXHDR
CXXHDR    /************************************************************/
CXXHDR    /*  HQSTF routines                                          */
CXXHDR    /************************************************************/
CXXHDR
C-----------------------------------------------------------------------
CXXHFW
CXXHFW  /**************************************************************/
CXXHFW  /*  HQSTF routines                                            */
CXXHFW  /**************************************************************/
CXXHFW
C-----------------------------------------------------------------------
CXXWRP
CXXWRP  /**************************************************************/
CXXWRP  /*  HQSTF routines                                            */
CXXWRP  /**************************************************************/
CXXWRP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
CXXHDR    int ihqvers();
C-----------------------------------------------------------------------
CXXHFW  #define fihqvers FC_FUNC(ihqvers,IHQVERS)
CXXHFW    int fihqvers();
C-----------------------------------------------------------------------
CXXWRP  //--------------------------------------------------------------
CXXWRP    int ihqvers()
CXXWRP    {
CXXWRP      return fihqvers();
CXXWRP    }
C-----------------------------------------------------------------------

C     ==========================
      integer function ihqvers()
C     ==========================

      implicit double precision(a-h,o-z)

      include 'hqstf.inc'

      ihqvers = iivers

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void hqwords(int &ntotal, int &nused);
C-----------------------------------------------------------------------
CXXHFW  #define fhqwords FC_FUNC(hqwords,HQWORDS)
CXXHFW    void fhqwords(int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void hqwords(int &ntotal, int &nused)
CXXWRP  {
CXXWRP    fhqwords(&ntotal,&nused);
CXXWRP  }
C-----------------------------------------------------------------------

C     ============================== 
      subroutine hqwords(ntot,nused)
C     ==============================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'hqstf.inc'
      include 'hqflags.inc'
      include 'hqstore.inc'
      
      ntot = nhqstor
      
      if(ihqini.eq.12345) then
        nused = nhused
      else
        nused = 0
      endif
      
      return        
      end

C-----------------------------------------------------------------------
CXXHDR    void hqparms(double *qmass, double &a, double &b);
C-----------------------------------------------------------------------
CXXHFW  #define fhqparms FC_FUNC(hqparms,HQPARMS)
CXXHFW    void fhqparms(double*, double*, double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void hqparms(double *qmass, double &a, double &b)
CXXWRP  {
CXXWRP    fhqparms(qmass,&a,&b);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================      
      subroutine hqparms(qmas,aq,bq)
C     ==============================

C--   Get mass and scale parameters

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'hqstf.inc'
      include 'hqflags.inc'
      include 'hqparms.inc'
      
      dimension qmas(*)
      
      if(ihqini.ne.12345) stop 
     +'HQPARMS: please first call HQFILLW or HQREADW'
      
      qmas(1) = hqmas(1)
      qmas(2) = hqmas(2)
      qmas(3) = hqmas(3)
      aq      = hqaaaa
      bq      = hqbbbb
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    double hqqfrmu(double qmu2);
C-----------------------------------------------------------------------
CXXHFW  #define fhqqfrmu FC_FUNC(hqqfrmu,HQQFRMU)
CXXHFW    double fhqqfrmu(double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double hqqfrmu(double qmu2)
CXXWRP  {
CXXWRP    return fhqqfrmu(&qmu2);
CXXWRP  }
C-----------------------------------------------------------------------

C     =======================================      
      double precision function hqqfrmu(qmu2)
C     =======================================

C--   Convert mu2 to Q2 = a mu2 + b

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'hqstf.inc'
      include 'hqflags.inc'
      include 'hqparms.inc'
      
      if(ihqini.ne.12345) stop 
     +'HQQFRMU: please first call HQFILLW or HQREADW'
      
      hqqfrmu = hqaaaa*qmu2 + hqbbbb
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    double hqmufrq(double Q2);
C-----------------------------------------------------------------------
CXXHFW  #define fhqmufrq FC_FUNC(hqmufrq,HQMUFRQ)
CXXHFW    double fhqmufrq(double*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  double hqmufrq(double Q2)
CXXWRP  {
CXXWRP    return fhqmufrq(&Q2);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     =====================================      
      double precision function hqmufrq(q2)
C     =====================================

C--   Convert Q2 to mu2 = (Q2 - b) / a

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'hqstf.inc'
      include 'hqflags.inc'
      include 'hqparms.inc'
      
      if(ihqini.ne.12345) stop
     +'HQMUFRQ: please first call HQFILLW or HQREADW'      
      
      hqmufrq = (q2-hqbbbb)/hqaaaa
      
      return
      end

C     ==============================      
      logical function LhqRvar(epsi)
C     ==============================

C--   True if renormalization scale is varied (within tolerance epsi)

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)

      include 'hqstf.inc'
      include 'hqflags.inc'
      include 'hqparms.inc'
      
      if(ihqini.ne.12345) stop 'LHQRVAR: HQSTF not initialized'
      
C--   Renormalization scale
      call GetABR(ar,br)
C--   Renor scale and fact scale equal or not thats the question
      if(abs(ar-1.D0).gt.epsi .or.
     +   abs(br     ).gt.epsi      ) then
        LhqRvar = .true.
      else
        LhqRvar = .false.
      endif      
      
      return
      end
      
C     ==============================      
      logical function LhqQvar(epsi)
C     ==============================

C--   True if Q2 scale is varied (within tolerance epsi)

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision(a-h,o-z)

      include 'hqstf.inc'
      include 'hqflags.inc'
      include 'hqparms.inc'
      
      if(ihqini.ne.12345) stop 'LHQQVAR: HQSTF not initialized'
      
      if(abs(hqaaaa-1.D0).gt.epsi .or.
     +   abs(hqbbbb     ).gt.epsi      ) then
        LhqQvar = .true.
      else
        LhqQvar = .false.
      endif      
      
      return
      end

C-----------------------------------------------------------------------
CXXHDR    void hswitch(int iset);
C-----------------------------------------------------------------------
CXXHFW  #define fhswitch FC_FUNC(hswitch,HSWITCH)
CXXHFW    void fhswitch(int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void hswitch(int iset)
CXXWRP  {
CXXWRP    fhswitch(&iset);
CXXWRP  }
C-----------------------------------------------------------------------
      
C     ========================      
      subroutine hswitch(jset)
C     ========================

C--   Pdf set idenitifier

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'hqstf.inc'
      include 'hqflags.inc'
      
      if(ihqini.ne.12345) stop
     +'HSWITCH: please first call HQFILLW or HQREADW'

      call getInt('mset',mset)
      
      if(jset.lt.1 .or. jset.gt.mset) stop
     + 'HSWITCH: iset not in range [1,mset]'
      if(jset.eq.2) stop
     + 'HSWITCH: cannot handle iset = 2 (polarised pdfs)'
      if(jset.eq.3) stop
     + 'HSWITCH: cannot handle iset = 3 (fragmentation functions)'
      if(jset.eq.4) stop
     + 'HSWITCH: cannot handle iset = 4 (custom/disabled) --> STOP'
       
      ihqpdf = jset

      return
      end

C-----------------------------------------------------------------------
CXXHDR    void hqstfun(int istf, int icbt, double *def, double *x, double *Q2, double *f, int n, int ichk);
C-----------------------------------------------------------------------
CXXHFW  #define fhqstfun FC_FUNC(hqstfun,HQSTFUN)
CXXHFW    void fhqstfun(int*, int*, double*, double*, double*, double*, int*, int*);
C-----------------------------------------------------------------------
CXXWRP//----------------------------------------------------------------
CXXWRP  void hqstfun(int istf, int icbt, double *def, double *x, double *Q2, double *f, int n, int ichk)
CXXWRP  {
CXXWRP   fhqstfun(&istf,&icbt,def,x,Q2,f,&n,&ichk);
CXXWRP  }
C-----------------------------------------------------------------------

C     ==============================================
      subroutine hqStFun(jstf,kcbt,def,x,q,f,n,ichk)
C     ==============================================

C--   Interpolate structure function
C--
C--   jstf       (in)  1=Fl 2=F2
C--   kcbt       (in)  1=c  2=b  3=t   <0 do not check nf
C--   def(-6:6)  (in)  quark linear combination
C--   x(n)       (in)  list of x values
C--   q(n)       (in)  list of Q2 (not mu2) values
C--   f(n)       (out) list of interpolated structure functions
C--   n          (in)  number of items in x, q, f 
C--   ichk       (in)  1 = check if x,Q2 is within grid boundaries

C--   Michiel Botje   h24@nikhef.nl 

      implicit double precision (a-h,o-z)
      
      include 'hqstf.inc'
      include 'hqparms.inc'
      include 'hqflags.inc'
      
      logical LhqRvar, LhqQvar
            
      common /passit/ dpass(-6:6),id0,jcbt

      dimension x(*), q(*), f(*)
      dimension def(-6:6)
      
C--   Set this the same as mpt0 in qcdnum.inc (not critical)      
      parameter ( nhqpnts = 5000 )
      dimension xx(nhqpnts),qmu(nhqpnts)
      
      dimension ioff(2)
      data ioff /0,3/
      dimension charge(3)
      save charge
      logical first
      save first
      data first /.true./

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,charge)
      
C--   Set subroutine name (must be cleared before exit)
      call setUmsg('HQSTFUN')
      if(first) then
        charge(1) = 4.D0/9.D0
        charge(2) = 1.D0/9.D0
        charge(3) = 4.D0/9.D0
        first     = .false.
      endif

C--   Get QCDNUM logical unit number, for error messages      
      call getInt('lunq',lun)

C--   Check if initialized
      if(ihqini.ne.12345) stop 
     +   'HQSTFUN: please first call HQFILLW or HQREADW --> STOP'

C--   Force evolution parameters of set ihqpdf
      call IdScope(0.D0,ihqpdf)

C--   Max number of interpolations
      call getint('mpt0',nmax)
      nmax = min(nmax,nhqpnts)

C--   Decode jstf
      iset = jstf/10
      if(iset.ne.0) call HSWITCH(iset)
      istf = jstf-10*iset

C--   Check input
      if(istf.lt.1 .or. istf.gt.2) stop
     +          'HQSTFUN: input ISTF not in range [1-2] --> STOP'
      icbt = abs(kcbt)
      if(icbt.lt.1 .or. icbt.gt.3) stop
     +          'HQSTFUN: input ICBT not in range [1-3] --> STOP'
C--   QCD order
      call getord(iord)
      if(iord.eq.3) stop
     +          'HQSTFUN: cannot handle NNLO --> STOP'
C--   Check flavour number scheme
      call GetCbt(nfix,q2c,q2b,q2t)
      if(nfix.eq.0) stop  !VFNS not allowed
     +          'HQSTFUN: pdfs evolved in the VFNS --> STOP'
      if(nfix.lt.0 .and. kcbt.gt.0) stop  !MFNS only when kcbt < 0
     +          'HQSTFUN: pdfs evolved in the MFNS --> STOP'
C--   Check number of flavors
      if(kcbt.eq.1 .and. nfix.ne.3) then
        write(lun,'(
     +  '' HQSTFUN: icbt = 1, but pdfs are not evolved with nf = 3''/
     +  ''          Set icbt to -1 to disable check on nf'')')
        stop
      elseif(kcbt.eq.2 .and. nfix.ne.4) then
        write(lun,'(
     +  '' HQSTFUN: icbt = 2, but pdfs are not evolved with nf = 4''/
     +  ''          Set icbt to -2 to disable check on nf'')')
        stop
      elseif(kcbt.eq.3 .and. nfix.ne.5) then
        write(lun,'(
     +  '' HQSTFUN: icbt = 3, but pdfs are not evolved with nf = 5''/
     +  ''          Set icbt to -3 to disable check on nf'')')
        stop
      elseif(nfix.eq.6) then
        write(lun,'('' HQSTFUN: pdfs evolved with nf = 6'')')
        stop  
      endif         
C--   Check Rscale and Qscale not both varied
      call GetVal('epsi',epsi)
      if(LhqRvar(epsi) .and. LhqQvar(epsi)) stop
     +          'HQSTFUN: cannot vary both Rscale and Q2 scale'      
C--   idtab0 is an array in hqparms.inc 
      id0  = idtab0(icbt+ioff(istf))
C--   No weights for this stf
      if(id0.eq.0) stop 
     +          'HQSTFUN: no weights available for this stf'
      
C--   Calculate Q2 limits (force Q2 > 0.5 GeV2)
      call grpars(nx,xmin,xmax,nq,qmin,qmax,iosp)
      q2min = max(hqaaaa*qmin+hqbbbb,0.5D0)
      q2max = hqaaaa*qmax+hqbbbb
C--   Check that Q2 is in range
      do i = 1,n     
        if((q(i).lt.q2min .or. q(i).gt.q2max) .and. ichk.ne.0) then
          write(lun,'('' HQSTFUN: Q2 = '',G13.5, 
     +    '' not in range [ '',G13.5,'','',G13.5,'' ]'')') 
     +    q(i),q2min,q2max
          stop
        endif
      enddo       
      
C--   Fill output array f in batches of nmax words        
      ipt = 0
      jj  = 0      
      do i = 1,n
        ipt     = ipt+1
        xx(ipt) = x(i)
C--     Convert Q2 to mu2
        qmu(ipt) = hqmufrq(q(i))
        if(ipt.eq.nmax) then
          call GetHqStf
     +         (icbt,charge,def,iord,id0,xx,qmu,f(jj*nmax+1),ipt,ichk)
          ipt = 0
          jj  = jj+1
        endif
      enddo
C--   Flush remaining ipt points
      if(ipt.ne.0) then
        call GetHqStf
     +         (icbt,charge,def,iord,id0,xx,qmu,f(jj*nmax+1),ipt,ichk)
      endif

C--   Clear subroutine name
      call clrUmsg
      
      return
      end           
      
C     ============================================================      
      subroutine GetHqStf(icbt,charge,def,iord,id0,x,qmu,f,n,ichk)
C     ============================================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'hqstf.inc'
      include 'hqparms.inc'
      include 'hqflags.inc'
      include 'hqstore.inc'
      
      external as1fun, as2fun
      
      dimension def(-6:6)
      
      dimension idwt(4)
      data idwt /0,0,0,0/
      
      dimension x(*), qmu(*), f(*)
      
      dimension charge(3),coef(0:12,3:6)

C--   Global identifier
      idg = 1000*jsetw+id0
C--   Pass interpolation points
      call fastIni(x,qmu,n,ichk)
      call fastClr(0)
C--   LO
C--   Select gluon
      do nf = 3,6
        coef(0,nf) = charge(icbt)
        do i = 1,12
          coef(i,nf) = 0.D0
        enddo
      enddo
      call fastSum(ihqpdf,coef,1)
      idwt(1) = idg
      call fastFxK(hqstor,idwt,1,2)
      call FastKin(2,as1fun)
      call fastCpy(2,3,1)
C--   NLO
      if(iord.eq.2) then
C--     Gluon
        idwt(1) = idg+1
        call fastFxK(hqstor,idwt,1,2)
        call FastKin(2,as2fun)
        call fastCpy(2,3,1)
        idwt(1) = idg+2
        call fastFxK(hqstor,idwt,1,2)
        call FastKin(2,as2fun)
        call fastCpy(2,3,1)
C--     Singlet
        do nf = 3,6
          coef(0,nf) = 0.D0
          coef(1,nf) = charge(icbt)
          do i = 2,12
            coef(i,nf) = 0.D0
          enddo
        enddo
        call fastSum(ihqpdf,coef,1)
        idwt(1) = idg+3
        call fastFxK(hqstor,idwt,1,2)
        call fastKin(2,as2fun)
        call fastCpy(2,3,1)
        idwt(1) = idg+4
        call fastFxK(hqstor,idwt,1,2)
        call fastKin(2,as2fun)
        call fastCpy(2,3,1)

C--     Proton
        do nf = 3,6
          coef(0,nf) = 0
          call efromqq(def,coef(1,nf),nf)
        enddo
        call fastSum(ihqpdf,coef,1)
        idwt(1) = idg+5
        call fastFxK(hqstor,idwt,1,2)
        call fastKin(2,as2fun)
        call fastCpy(2,3,1)
        idwt(1) = idg+6
        call fastFxK(hqstor,idwt,1,2)
        call fastKin(2,as2fun)
        call fastCpy(2,3,1)
      endif
C--   Interpolate
      call fastFxq(3,f,n)
      
      return
      end                  
                                         
C     ================================================ 
      double precision function as1fun(ix,iq,nf,ithrs)
C     ================================================

C--   Returns as/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'hqflags.inc'
      
      jx    = ix      !Avoid compiler warning
      jf    = nf      !Avoid compiler warning
      jthrs = ithrs   !Avoid compiler warning 
      
      as1fun = altabn(ihqpdf,iq,1,ierr)
      
      return
      end
      
C     ================================================ 
      double precision function as2fun(ix,iq,nf,ithrs)
C     ================================================

C--   Returns (as/2pi)^2

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'hqflags.inc'
      
      jx    = ix      !Avoid compiler warning
      jf    = nf      !Avoid compiler warning
      jthrs = ithrs   !Avoid compiler warning 
      
      as2fun = altabn(ihqpdf,iq,2,ierr)
      
      return
      end      

C--   ----------------------------------------------------- 
C--   The following routines are slow calculations used for
C--   prototyping; they may also be easier to understand
C--   -----------------------------------------------------
      
C     ==============================================
      subroutine hqSlowF(istf,icbt,def,x,q,f,n,ichk)
C     ==============================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      include 'hqparms.inc'
      include 'hqflags.inc'
      
      external dhqFij
      
      common /passit/ dpass(-6:6),id0,jcbt
      
      dimension x(*), q(*), f(*)
      dimension def(-6:6)
      dimension ioff(2)
      data ioff /0,3/
      
C--   Set subroutine name (must be cleared before exit)
      call setUmsg('HQSLOWF')
C--   Check if initialized
      if(ihqini.ne.12345) stop
     +  'HQSLOWF: please first call HQFILLW or HQREADW'
C--   Check input
      if(istf.lt.1 .or. istf.gt.3) stop
     +          'HQSLOWF: input ISTF not in range [1-3]'
      if(icbt.lt.1 .or. icbt.gt.3) stop
     +          'HQSLOWF: input ICBT not in range [1-3]'
C--   Force evolution parameters of set ihqpdf
      call IdScope(0.D0,ihqpdf)
C--   Copy to common block; idtab0 is an array hqparms.inc 
      id0  = idtab0(icbt+ioff(istf))
      jcbt = icbt
      do i = -6,6
        dpass(i) = def(i)
      enddo
C--   No weights for this stf
      if(id0.eq.0) stop 
     +          'HQSLOWF: no weights available for this stf'          
C--   Get structure function            
      call stfunxq(dhqFij,x,q,f,n,ichk)
C--   Clear subroutine name
      call clrUmsg
      
      return
      end                                        

C     =======================================
      double precision function dhqFij(ix,iq)
C     =======================================

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z) 

      include 'hqstf.inc'
      include 'hqstore.inc'
            
      common /passit/ dpass(-6:6),id0,jcbt

      dhqFij = 0.D0

      if(ix.eq.-1) then
        continue
      else
        dhqFij = dhqGetF(hqstor,id0,jcbt,dpass,ix,iq)
      endif
      
      return
      end      

C     =======================================================
      double precision function dhqGetF(w,id0,icbt,def,ix,iq)
C     =======================================================

C--   Input w          store filled with tables
C--         id0        index of first table C^0_{k,G}
C--         icbt       1=c, 2=b, 3=t
C--         def(-6:6)  quark flavor decomposition
C--         ix,iq      grid point 

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      include 'hqstf.inc'
      include 'hqparms.inc'
      include 'hqflags.inc'
      include 'hqstore.inc'
      
      logical first
      save first
      data first /.true./
      dimension charge(3),def(-6:6)
      save charge

      dimension w(*), evec(12)

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(first,charge)

      bla = def(1) !avoid compiler warning

C--   Global identifier
      idw = 1000*jsetw+id0

C--   Should be 1=unpolarised, 4=custom or 5-9=external
      jset = ihqpdf
      
C--   Initialize      
      if(first) then
        charge(1) = 4.D0/9.D0
        charge(2) = 1.D0/9.D0
        charge(3) = 4.D0/9.D0
        first     = .false.
      endif
      
C--   QCDNUM
      call getval('null',qnull)
      call getval('epsi',epsi)
C--   QCD order
      call getord(iord)

C--   No NNLO, thank you
      if(iord.eq.3) then
        dhqGetF = qnull
        return
      endif          
      
C--   Initialize
      F1  = 0.D0    !LO contribution
      F2  = 0.D0    !NLO gluon
      F3  = 0.D0    !NLO singlet
      F4  = 0.D0    !NLO proton
      as1 = 0.D0    !alfas
      as2 = 0.D0    !alfas2
      idg = ipdftab(jset,0)
      ids = ipdftab(jset,1)

C--   LO
      if(iord.le.2) then
        F1      = FcrossK(w,idw,jset,idg,ix,iq)
        as1     = altabn(ihqpdf,iq,1,ierr)
      endif          
      
C--   NLO
      if(iord.eq.2) then
        qmu2 = qfrmiq(iq)
        as2  = altabn(ihqpdf,iq,2,ierr)
        nfl  = nfrmiq(0,iq,ithrs)
        call efromqq(def,evec,nfl)
        F2   =      FcrossK(w,idw+1,jset,idg,ix,iq)
        F2   = F2 + FcrossK(w,idw+2,jset,idg,ix,iq)
        F3   =      FcrossK(w,idw+3,jset,ids,ix,iq)
        F3   = F3 + FcrossK(w,idw+4,jset,ids,ix,iq)
        do ipdf = 1,nfl
          jpdf0 = ipdftab(jset,ipdf)
          jpdf6 = ipdftab(jset,ipdf+6)
          F4 = F4 + evec(ipdf  )*FcrossK(w,idw+5,jset,jpdf0,ix,iq)
          F4 = F4 + evec(ipdf  )*FcrossK(w,idw+6,jset,jpdf0,ix,iq)
          F4 = F4 + evec(ipdf+6)*FcrossK(w,idw+5,jset,jpdf6,ix,iq)
          F4 = F4 + evec(ipdf+6)*FcrossK(w,idw+6,jset,jpdf6,ix,iq)
        enddo          
      endif
      
C--   Assemble structure function      
      dhqGetF = charge(icbt)*(as1*F1 + as2*(F2+F3)) + as2*F4      
      
      return
      end                                      
