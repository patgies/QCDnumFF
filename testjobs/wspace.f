
C=======================================================================
      program wspace
C=======================================================================
C--
C--   Examples described in the WSTORE write-up
C--
C=======================================================================

      implicit double precision (a-h,o-z)
      integer wtable

      parameter( nw = 20000 )
      parameter( nt = 10    )
      dimension w(nw)

      dimension k3(0:3)

C--   Inline pointer function
      IAijk(i,j,k) = k3(0)+k3(1)*i+k3(2)*j+k3(3)*k+ixq

      write(6,'(/A,I8,A)') ' You are using WSTORE version ',
     &                                        iws_Version(),' ---------'

      iaS = iws_WsInit ( w, nw, nt, 'Please increase parameter NW' )

      write(6,'(/A)') ' Example of a table-set ------------------------'
      iaS = mytab(w)
      call tabinfo(w,iaS,iax,iaq,ixq,k3)
      write(6,'(A,2F8.0)') ' Limits of x-bin 10 :', w(iax+10),w(iax+11)
      write(6,'(A,2F8.0)') ' Limits of q-bin  5 :', w(iaq+5),w(iaq+6)
      write(6,'(A, F8.0)') ' Value of T(10,5,4) :', w(IAijk(10,5,4))

      write(6,'(/A)') ' Use home-made pointer function iP3 ------------'
      write(6,'(A, F8.0)') ' Value of T(10,5,4) :', w(iP3(w,ixq,10,5,4))

      write(6,'(/A)') ' Fast 3-dim loop construct ---------------------'
      i1 = 11
      i2 = 13
      j1 = 15
      j2 = 17
      k1 = 3
      k2 = 4
      call fastloop(w,ixq,i1,i2,j1,j2,k1,k2)

      write(6,'(/A)') ' Convolution integral --------------------------'
      iaW = wtable(w)                        !delta function
      ix  = 10
      iq  = 20
      nf  =  4
      jaF = iP3(w,ixq,1,iq,nf) !Address of F(1,iq,nf)
      jaW = iP3(w,iaW,1,ix,nf) !Address of W(1,ix,nf)
      FxW = dmb_VdotV(w(jaF),w(jaW),ix)      !Convolution at (ix,iq,nf)
      write(6,'(A, F8.0)') ' T(10,20,4) x W     :', FxW

      write(6,'(/A)') ' Initialise table ------------------------------'
      call IniTab(w,ixq,-1.D0)
      call fastloop(w,ixq,i1,i2,j1,j2,k1,k2)
    
      end

C     =========================
      integer function mytab(w)
C     =========================

      implicit double precision (a-h,o-z)
      dimension w(*), imi(3), ima(3)
      data imi/1,1,3/, ima/50,25,6/

C--   Create table-set
      nt       = iws_Ntags(w)
      if(nt.lt.3) stop 'MYTAB: increase parameter NT to at least 3'
      ias      = iws_NewSet(w)
      ima(1)   = ima(1)+1
      ima(2)   = ima(2)+1
      iax      = iws_WTable(w,imi(1),ima(1),1)
      iaq      = iws_WTable(w,imi(2),ima(2),1)
      ima(1)   = ima(1)-1
      ima(2)   = ima(2)-1
      ixq      = iws_WTable(w,imi,ima,3)
      iat      = iws_IaFirstTag(w,ias)
      w(iat)   = dble(iax-ias)
      w(iat+1) = dble(iaq-ias)
      w(iat+2) = dble(ixq-ias)
      mytab    = ias

C--   Fill x-table
      ia = iws_BeginTbody(w,iax)-1
      do i = imi(1),ima(1)+1
        ia    = ia+1
        w(ia) = dble(i)
      enddo
      ja = iws_EndTbody(w,iax)
      if(ia.ne.ja) stop 'MYTAB: something wrong with x-indexing'
C--   Fill q-table
      ia = iws_BeginTbody(w,iaq)-1
      do i = imi(2),ima(2)+1
        ia    = ia+1
        w(ia) = dble(i)
      enddo
      ja = iws_EndTbody(w,iaq)
      if(ia.ne.ja) stop 'MYTAB: something wrong with q-indexing'
C--   Fill xq-table
      ia = iws_BeginTbody(w,ixq)-1
      do k = imi(3),ima(3)
        do j = imi(2),ima(2)
          do i = imi(1),ima(1)
            ia    = ia+1
            w(ia) = dble(10000*i+100*j+k)
          enddo
        enddo
      enddo
      ja = iws_EndTbody(w,ixq)
      if(ia.ne.ja) stop 'MYTAB: something wrong with xq-indexing'

      return
      end

C     ========================================
      subroutine tabinfo(w,ias,iax,iaq,ixq,k3)
C     ========================================

      implicit double precision (a-h,o-z)
      dimension w(*), k3(0:3)

      iat   = iws_IaFirstTag(w,ias)
      itx   = int(w(iat))+ias
      itq   = int(w(iat+1))+ias
      ixq   = int(w(iat+2))+ias
      iax   = iws_BeginTbody(w,itx)-1
      iaq   = iws_BeginTbody(w,itq)-1
      ikk   = iws_IaKARRAY(w,ixq)
      k3(0) = int(w(ikk))
      k3(1) = int(w(ikk+1))
      k3(2) = int(w(ikk+2))
      k3(3) = int(w(ikk+3))

      return
      end

C     ========================
      subroutine K3(w, ia, kk)
C     ========================
      dimension kk(0:4)
      double precision w(*)
      iak   = iws_IaKARRAY(w,ia)
      kk(0) = int(w(iak))
      kk(1) = int(w(iak+1))
      kk(2) = int(w(iak+2))
      kk(3) = int(w(iak+3))
      kk(4) = iws_FingerPrint(w,ia)
      return
      end

C     ====================================
      integer function iP3(w, ia, i, j, k)
C     ====================================
      double precision w(*)
      dimension kk(0:4)
      save kk
      ifp = iws_FingerPrint(w,ia)
      if(kk(4).ne.ifp) call K3(w, ia, kk)
      ip = kk(0)+kk(1)*i+kk(2)*j+kk(3)*k
      iP3 = ip + ia
      return
      end


C     ==============================================
      subroutine fastloop(w,iatab,i1,i2,j1,j2,k1,k2)
C     ==============================================

      implicit double precision (a-h,o-z)

      dimension w(*), kk(0:4)

      call K3(w,iatab,kk)
      inci = kk(1)
      incj = kk(2)
      inck = kk(3)
      ia   = iP3(w,iatab,i1,j1,k1)
      do i = i1,i2
        ja = ia
        do j = j1,j2
          ka = ja
          do k = k1,k2
            write(6,'(A,3I3,A,F8.0)') ' T(',i,j,k,') = ',w(ka)
            ka = ka + inck
          enddo
          ja = ja + incj
        enddo
        ia = ia + inci
      enddo

      return
      end

C     ==========================
      integer function wtable(w)
C     ==========================

      implicit double precision (a-h,o-z)
      dimension w(*), imi(3), ima(3)
      data imi/1,1,3/, ima/50,50,6/
      dimension k3(0:3)

C--   Inline pointer function
      iW3(i,j,k) = k3(0) + k3(1)*i + k3(2)*j + k3(3)*k + iaw

C--   Create table-set
      nt       = iws_Ntags(w)
      if(nt.lt.1) stop 'WTABLE: increase parameter NT to at least 1'
      ias      = iws_NewSet(w)
      iaw      = iws_WTable(w,imi,ima,3)
      iat      = iws_IaFirstTag(w,ias)
      w(iat)   = dble(iaw-ias)
C--   Pointer coefficients
      iak      = iws_IaKARRAY(w,iaw)
      k3(0)    = int(w(iak))
      k3(1)    = int(w(iak+1))
      k3(2)    = int(w(iak+2))
      k3(3)    = int(w(iak+3))
C--   Return weight table address
      wtable   = iaw

C--   Fill weight table with delta-function
      do k = imi(3),ima(3)
        do j = imi(2),ima(2)
          w(iW3(j,j,k)) = 1.D0
        enddo
      enddo

      return
      end

C     =============================
      subroutine IniTab(w, ia, val)
C     =============================

      double precision w(*), val

      i1 = iws_BeginTbody(w,ia)
      i2 = iws_EndTbody(w,ia)
      do i = i1,i2
        w(i) = val
      enddo

      return
      end

