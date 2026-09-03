
C--   This is the file srcConvolution.f containing convolution routines
C--
C--   double precision function dqcFcrossC(w,id,idi,idt,jy,it)
C--   double precision function dqcFcrossK(w,idw,wf,idf,iy,it)
C--   double precision function dqcFcrossF(ww,idw,wa,ida,wb,idb,iy,it)

C==   ==============================================================
C==   Convolution ==================================================
C==   ==============================================================

C     ========================================================
      double precision function dqcFcrossC(w,id,idi,idt,jy,it)
C     ========================================================

C--   Dummy routine

      implicit double precision (a-h,o-z)

      dimension w(*)

      dum        = w(1)
      idum       = id
      idum       = idi
      idum       = idt
      idum       = jy
      idum       = it
      dqcFcrossC = 0.D0

      return
      end
      
C     =========================================================
      double precision function dqcFcrossK(ww,idw,wf,idf,iy,it)
C     =========================================================

C--   Calculate convolution F cross K

C--   ww    (in)  :  store with weight tables
C--   idw   (in)  :  identifier of K in global format
C--   wf    (in)  :  store with pdf tables
C--   idf   (in)  :  identifier of F in global format
C--   iy,it (in)  :  grid point           
      
      implicit double precision(a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'

      dimension ww(*),wf(*),coef(mxx0)

      iz = izfit5( it)
      nf = itfiz5(-iz)

C--   Get spline coefficients
      call sqcGetSplA(wf,idf,iy,iz,ig,iyg,coef)
      
C--   Weight table base address
      iwk = iqcGaddr(ww,1,abs(it),nf,ig,idw)-1

C--   Convolute               
      fxk = 0.D0
      do j = 1,iyg
        fxk = fxk + coef(j)*ww(iwk+iyg-j+1)
      enddo
                      
      dqcFcrossK = fxk
      
      return
      end

C     ================================================================
      double precision function dqcFcrossF(ww,idw,wa,ida,wb,idb,iy,it)
C     ================================================================

C--   Calculate convolution f_a cross f_b

C--   ww    (in)  :  store with weight tables
C--   idw   (in)  :  identifier table with fxf weights in global format
C--   wa    (in)  :  store with pdfs
C--   ida   (in)  :  identifier of pdf f_a in global format
C--   wb    (in)  :  store with pdfs
C--   idb   (in)  :  identifier of pdf f_b in global format
C--   iy,it (in)  :  grid point           
      
      implicit double precision(a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'point5.inc'
      include 'qstor7.inc'
      
      dimension ww(*),wa(*),wb(*),coefa(mxx0),coefb(mxx0)

      iz = izfit5( it)
      nf = itfiz5(-iz)

C--   Get spline coefficients
      call sqcGetSplA(wa,ida,iy,iz,ig,iyg,coefa)
      call sqcGetSplA(wb,idb,iy,iz,ig,iyg,coefb)
      
C--   Weight table base address
      iwx = iqcGaddr(ww,1,abs(it),nf,ig,idw)-1

C--   Convolute               
      fxf = 0.D0
      do j = 1,iyg
        Aj = coefa(j)
        do k = 1,iyg-j+1
          Bk  = coefb(k) 
          fxf = fxf + Aj*Bk*ww(iwx+iyg-j-k+2)
        enddo
      enddo
                      
      dqcFcrossF = fxf
      
      return
      end

