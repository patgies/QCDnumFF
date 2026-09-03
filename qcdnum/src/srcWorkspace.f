
C--  This is the file srcWorkSpace.f with memory management routines
C--
C--   integer function iqcGaddr(ww,i,j,k,l,m)
C--   integer function iqcGCadr(ww,i,j,k,l,m)
C--   integer function iqcG1ijk(ww,i,j,k)
C--   integer function iqcG2ijkl(ww,i,j,k,l)
C--   integer function iqcG3ijkl(ww,i,j,k,l)
C--   integer function iqcG4ijklm(ww,i,j,k,l,m)
C--   integer function iqcG5ijk(ww,i,j,k)
C--   integer function iqcG6ij(ww,i,j)
C--   integer function iqcG7ij(ww,i,j)
C--   integer function iqcGSij(ww,i,j)
C--   integer function iqcGSi(ww,i)
C--
C--   integer function iqcWaddr(wa,i,j,k,l,m)
C--   integer function iqcWCadr(wa,i,j,k,l,m)
C--   integer function iqcW1ijk(wa,i,j,k)
C--   integer function iqcW2ijkl(wa,i,j,k,l)
C--   integer function iqcW3ijkl(wa,i,j,k,l)
C--   integer function iqcW4ijklm(wa,i,j,k,l,m)
C--   integer function iqcW5ijk(wa,i,j,k)
C--   integer function iqcW6ij(wa,i,j)
C--   integer function iqcW7ij(wa,i,j)
C--   integer function iqcWSij(wa,i,j)
C--   integer function iqcWSi(wa,i)
C--
C--   integer function iqcFirstWordOfSet(ww,kset)
C--   integer function iqcGetNumberOfWords(ww)
C--   integer function iqcGetNumberOfSets(ww)
C--   integer function iqcGetNumberOfParams(ww,kset)
C--   integer function iqcFirstWordOfParams(ww,kset)
C--   integer function iqcGetNumberOfUparam(ww,kset)
C--   integer function iqcFirstWordOfUparam(ww,kset)
C--   integer function iqcGetNumberOfTables(ww,kset,itype)
C--   integer function iqcGetNumberOfTabsWa(wa,itype)
C--   integer function iqcSgnNumberOfTables(ww,kset,itype)
C--   integer function iqcSgnNumberOfTabsWa(wa,itype)
C--   subroutine sqcGetLimits(ww,id,imin,imax,jmax)
C--   subroutine sqcGetLimsWa(wa,id,imin,imax,jmax)
C--   subroutine sqcGetLimSpl(ww,id,iosp1,iosp2)
C--   subroutine sqcGetLimSpa(wa,id,iosp1,iosp2)
C--   logical function lqcIsDouble(ww,id)
C--   integer function iqcFirstWordOfTable(ww,id)
C--   integer function iqcGetTabLeng(ww,id,ndim)
C--   integer function iqcGetTbLenWa(wa,id,ndim)
C--   integer function iqcGetSetNumber(id)
C--   integer function iqcGetLocalId(id)
C--   logical function lqcWpartitioned(ww)
C--   logical function lqcIsetExists(ww,kset)
C--   logical function lqcIdExists(ww,id)
C--
C--   subroutine sqcValidate(ww,id)
C--   subroutine sqcInvalidate(ww,id)
C--   logical function lqcIsFilled(ww,id)
C--   subroutine sqcSetMin6(ww,id,itmin)
C--   subroutine sqcGetMin6(ww,id,itmin)
C--
C--   subroutine sqcTcopyType5(ww,id,it1,it2)
C--   subroutine sqcCopyType5(w1,id1,w2,id2)
C--   subroutine sqcCopyType6(w1,id1,w2,id2)
C--   subroutine sqcCopyType7(w1,id1,w2,id2)
C--
C--   subroutine sqcMakeTab(ww,nw,itypes,npar,nusr,new,kset,nwords)
C--   subroutine sqcBookSet(wa,nw,itypes,mpar,musr,nwords,ierr)

C--   ==================================================================
C--   Layout of QCDNUM dynamic memory
C--   ==================================================================
C--
C--   A dynamic memory is a large linear array, declared in the calling
C--   routine as w(nw).  This dynamic memory is structured as follows:
C--
C--   |header|set1|set2|.....|setn|
C--
C--   The max number of sets in memory is given by mst0 in qcdnum.inc
C--
C--   Layout of the header:
C--
C--   1            : init code 654321
C--   2            : total number of words used
C--   3            : number of checksums (mchk0)
C--   next mchk0   : checksums
C--   4+mchk0      : current number of sets
C--   5+mchk0      : max number of sets (mst0)
C--   next mst0    : w(5+mchk0+kset) = start address of kset = 1,...,mst0
C--   6+mchk0+mst0 : start address of set mst0+1 (dummy address)
C--   7+mchk0+mst0 : first word of first set
C--
C--   Pointer information in each set is relative to first word of set
C--
C--   -->       address_in_w = address_in_kset + w(5+mchk0+kset)-1
C--
C--   Layout of a set:
C--
C--   word  1          : partition code 123456
C--         2          : number of words used for this set
C--         3          : number of parameter info words (npar)
C--         4          : number of user info words (nusr)
C--         5          : number of table types (mtyp0)
C--         next npar  : parameter info words
C--         next nusr  : user info words
C--         next mtyp0 : pointer to first word of each type (0 = absent)
C--         next 20    : imin,imax,karr,ntab of first type 6-dim tables
C--         next  7    : jmin,jmax,ksat of first 2-dim satellite table
C--         next words : tables of the first type + satellite table
C--         next 20    : imin,imax,karr,ntab of second type 6-dim tables
C--         next  7    : jmin,jmax,ksat of second 2-dim satellite table
C--         next words : tables of the second type + satellite table
C--         etc.
C--
C--   Each table has an associated satellite table where information
C--   is stored for each id, e.g. validation flag, cuts, etc.
C--
C--   Indexing of tables: always 6 indices with one or more dummies
C--                   1        2        3       4          5      6
C--   type 1 : w( iy[1,ny] it[1,1]  nf[3,3] ig[1,ng] id[101,199] iosp )
C--   type 2 : w( iy[1,ny] it[1,1]  nf[3,6] ig[1,ng] id[201,299] iosp )
C--   type 3 : w( iy[1,ny] it[1,nt] nf[3,3] ig[1,ng] id[301,399] iosp )
C--   type 4 : w( iy[1,ny] it[1,nt] nf[3,6] ig[1,ng] id[401,499] iosp )
C--   type 5 : w( iy[1,ny] it[1,nt] nf[3,3] ig[1,1]  id[501,599] iosp )
C--   type 6 : w( iy[1,1]  it[1,nt] nf[3,3] ig[1,1]  id[601,699] iosp )
C--   type 7 : w( iy[1,1]  it[1,nt] nf[3,3] ig[1,1]  id[701,799] iosp )
C--
C--   Satellite tables: always 2 indices with perhaps dummy first index
C--                   1         2
C--   type 1 : w( ii[1,1] id[101,199] ) 1=validation
C--   type 2 : w( ii[1,1] id[201,299] ) 1=validation
C--   type 3 : w( ii[1,1] id[301,399] ) 1=validation
C--   type 4 : w( ii[1,1] id[401,499] ) 1=validation
C--   type 5 : w( ii[1,4] id[501,599] ) 1=validation 2=xmi 3=qmi 4=qma
C--   type 6 : w( ii[1,2] id[601,699] ) 1=validation 2=itmin
C--   type 7 : w( ii[1,1] id[701,799] ) 1=validation
C--
C--   The range of the iosp index depends on the current interpolation
C--   order ioy2 in the qgrid2 common block, and on if the table is a
C--   splitting function or a coefficient function
C--
C--   Current spline order ioy2        |  2  |  3  |
C--                                    +-----+-----+
C--   iosp range splitting function    | 2-2 | 2-3 |
C--   iosp range coefficient function  | 2-2 | 3-3 |

C==   ==================================================================
C==   Global address functions  ========================================
C==   ==================================================================

C--   A global address function operates on the store and returns the
C--   address with respect to the first word of the store. Identifiers
C--   must be entered in GLOBAL format: id_global = 1000*kset + id_local
C--   
C--   It follows that w should NOT be aligned in a subroutine call:
C--
C--           call SUB(ww,...)
C--
C--           subroutine SUB(ww,...)
C--           dimension ww(*)
C--           iword = iqcG5ijk(ww,iy,it,id)         !absolute address
C--           value = ww(iword)                     !pdftable(iy,it,id)

C     =======================================
      integer function iqcGaddr(ww,i,j,k,l,m)
C     =======================================

C--                        i  j  k  l  m
C--   iwabs = iqcGaddr(ww,iy,it,nf,ig,id)
C--
C--   On input, the store ww must not be aligned
C--   Works on tables of any type
C--   Returns 0 if ww not initialised or if itype = m/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(m.lt.0) stop 'iqcGadr m < 0'
      kset        =     m/1000
      id          =     m-1000*kset
      iwset       =     iqcFirstWordOfSet(ww,kset)
      iqcGaddr    =     iqcWaddr(ww(iwset),i,j,k,l,id)
      if(iqcGaddr.eq.0) return
      iqcGaddr    =     iqcGaddr+iwset-1

      return
      end 

C     =======================================
      integer function iqcGCadr(ww,i,j,k,l,m)
C     =======================================

C--   As iqcGaddr but with array boundary check 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(m.lt.0) stop 'iqcGCadr m < 0'
      kset        =     m/1000
      id          =     m-1000*kset
      iwset       =     iqcFirstWordOfSet(ww,kset)
      iqcGCadr    =     iqcWCadr(ww(iwset),i,j,k,l,id)
      if(iqcGCadr.eq.0) return
      iqcGCadr    =     iqcGCadr+iwset-1

      return
      end 

C     ===================================
      integer function iqcG1ijk(ww,i,j,k)
C     ===================================

C--                                     i  j  k
C--   Indices:  iaddress = iqcG1ijk(ww,iy,ig,id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(k.lt.0) stop 'iqcG1ijk k < 0'
      kset        =     k/1000
      id          =     k-1000*kset
      iwset       =     iqcFirstWordOfSet(ww,kset)
      iqcG1ijk    =     iqcW1ijk(ww(iwset),i,j,id)
      if(iqcG1ijk.eq.0) return
      iqcG1ijk    =     iqcG1ijk+iwset-1

      return
      end

C     ======================================
      integer function iqcG2ijkl(ww,i,j,k,l)
C     ======================================

C--                                      i  j  k  l
C--   Indices:  iaddress = iqcG2ijkl(ww,iy,nf,ig,id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(l.lt.0) stop 'iqcG2ijkl l < 0'
      kset         =     l/1000
      id           =     l-1000*kset
      iwset        =     iqcFirstWordOfSet(ww,kset)
      iqcG2ijkl    =     iqcW2ijkl(ww(iwset),i,j,k,id)
      if(iqcG2ijkl.eq.0) return
      iqcG2ijkl    =     iqcG2ijkl+iwset-1

      return
      end 
      
C     ======================================
      integer function iqcG3ijkl(ww,i,j,k,l)
C     ======================================

C--                                      i  j  k  l
C--   Indices:  iaddress = iqcG3ijkl(ww,iy,it,ig,id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(l.lt.0) stop 'iqcG3ijkl l < 0'
      kset         =     l/1000
      id           =     l-1000*kset
      iwset        =     iqcFirstWordOfSet(ww,kset)
      iqcG3ijkl    =     iqcW3ijkl(ww(iwset),i,j,k,id)
      if(iqcG3ijkl.eq.0) return
      iqcG3ijkl    =     iqcG3ijkl+iwset-1

      return
      end       

C     =========================================
      integer function iqcG4ijklm(ww,i,j,k,l,m)
C     =========================================

C--                                       i  j  k  l  m
C--   Indices:  iaddress = iqcG4ijklm(ww,iy,it,nf,ig,id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(m.lt.0) stop 'iqcG4ijklm m < 0'
      kset          =     m/1000
      id            =     m-1000*kset
      iwset         =     iqcFirstWordOfSet(ww,kset)
      iqcG4ijklm    =     iqcW4ijklm(ww(iwset),i,j,k,l,id)
      if(iqcG4ijklm.eq.0) return
      iqcG4ijklm    =     iqcG4ijklm+iwset-1

      return
      end 

C     ===================================
      integer function iqcG5ijk(ww,i,j,k)
C     ===================================

C--                                     i  j  k
C--   Indices:  iaddress = iqcG5ijk(ww,iy,it,id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(k.lt.0) stop 'iqcG5ijk k < 0'
      kset        =     k/1000
      id          =     k-1000*kset
      iwset       =     iqcFirstWordOfSet(ww,kset)
      iqcG5ijk    =     iqcW5ijk(ww(iwset),i,j,id)
      if(iqcG5ijk.eq.0) return
      iqcG5ijk    =     iqcG5ijk+iwset-1

*mbC--   Call address function with array boundary check
*mb      itest = iqcGCadr(ww,i,j,3,1,k)
*mbC--   Check that the call to iqcGCadr is OK
*mb      if(itest.ne.iqcG5ijk) stop 'Gadr .ne. G5ijk'

      return
      end     

C     ================================
      integer function iqcG6ij(ww,i,j)
C     ================================

C--                                    i  j
C--   Indices:  iaddress = iqcG6ij(ww,it,id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(j.lt.0) stop 'iqcG6ij j < 0'
      kset       =     j/1000
      id         =     j-1000*kset
      iwset      =     iqcFirstWordOfSet(ww,kset)
      iqcG6ij    =     iqcW6ij(ww(iwset),i,id)
      if(iqcG6ij.eq.0) return
      iqcG6ij    =     iqcG6ij+iwset-1

      return
      end

C     ================================
      integer function iqcG7ij(ww,i,j)
C     ================================

C--                                    i  j
C--   Indices:  iaddress = iqcG7ij(ww,it,id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(j.lt.0) stop 'iqcG7ij j < 0'
      kset       =     j/1000
      id         =     j-1000*kset
      iwset      =     iqcFirstWordOfSet(ww,kset)
      iqcG7ij    =     iqcW7ij(ww(iwset),i,id)
      if(iqcG7ij.eq.0) return
      iqcG7ij    =     iqcG7ij+iwset-1

      return
      end

C     ================================
      integer function iqcGSij(ww,i,j)
C     ================================

C--   Satellite table                 i  j
C--   Indices:  iaddress = iqcGSij(ww,i,id)

      implicit double precision (a-h,o-z)

      dimension ww(*)

      if(j.lt.0) stop 'iqcGSij j < 0'
      kset       =     j/1000
      id         =     j-1000*kset
      iwset      =     iqcFirstWordOfSet(ww,kset)
      iqcGSij    =     iqcWSij(ww(iwset),i,id)
      if(iqcGSij.eq.0) return
      iqcGSij    =     iqcGSij+iwset-1

      return
      end

C     =============================
      integer function iqcGSi(ww,i)
C     =============================

C--   Satellite table                 i
C--   Indices:  iaddress = iqcGSi(ww,id)

      implicit double precision (a-h,o-z)

      dimension ww(*)

      if(i.lt.0) stop 'iqcGSi i < 0'
      kset       =     i/1000
      id         =     i-1000*kset
      iwset      =     iqcFirstWordOfSet(ww,kset)
      iqcGSi     =     iqcWSi(ww(iwset),id)
      if(iqcGSi.eq.0)  return
      iqcGSi     =     iqcGSi+iwset-1

      return
      end

C==   ==================================================================
C==   Local address functions  =========================================
C==   ==================================================================

C--   A local address function operates on a table-set and returns the
C--   address with respect to the first word of this set. Thus it is
C--   necessary to somehow align the store onto the first word of a set.
C--
C--   The table identifier must be given in 3-digit local format
C--   
C--   Here is an example of alignment INSIDE a subroutine:
C--
C--           call SUB1(ww,kset,...)
C--
C--           subroutine SUB1(ww,kset,...)
C--           dimension ww(*)
C--           iwset = iqcFirstWordOfSet(ww,kset)    !first word of kset
C--           iwrel = iqcW5ijk(ww(iwset),iy,it,id)  !relative address
C--           iwabs = iwrel+iwset-1                 !absolute address
C--           value = ww(iwabs)                     !pdftable(iy,it,id)
C--
C--   Here is an example of alignment OUTSIDE a subroutine:
C--
C--           iwset = iqcFirstWordOfSet(ww,kset)    !first word of kset
C--           call SUB2(ww(iwset),...)
C--
C--           subroutine SUB2(wa,...)
C--           dimension wa(*)
C--           iword = iqcW5ijk(wa,iy,it,id)         !local address
C--           value = wa(iword)                     !pdftable(iy,it,id)
C--
C--   The function iqcFirstWordOfSet returns the address of kset when w
C--   is not aligned, and returns 1 when w is aligned (on any set):
C--
C--           ifirst   = iqcFirstWordOfSet(ww,kset)   !address of kset
C--           jfirst   = iqcFirstWordOfSet(ww(ifirst),kset) !always 1
C--
C--   Note from this that SUB1 is more robust since it does not matter
C--   if ww is aligned or not; BOTH calls below will work
C--
C--           iwset = iqcFirstWordOfSet(ww,kset)
C--           call SUB1(ww,kset,...)               !do not align w
C--           call SUB1(ww(iwset),kset,...)        !align w
C--
C--   SUB2 is faster but less robust since ww MUST be aligned on entry
C--
C--   Words relative to iw (see routines below)
C--
C--   iw + 0  = min iy
C--      + 1  = max iy
C--      + 2  = min it
C--      + 3  = max it
C--      + 4  = min nf
C--      + 5  = max nf
C--      + 6  = min ig
C--      + 7  = max ig
C--      + 8  = min id
C--      + 9  = max id
C--      +10  = min io
C--      +11  = max io
C--      +12  = karr(0)
C--      +13  = karr(1)   (*i = iy)
C--      +14  = karr(2)   (*j = it)
C--      +15  = karr(3)   (*k = nf)
C--      +16  = karr(4)   (*l = ig)
C--      +17  = karr(5)   (*m = id)
C--      +18  = karr(6)   (*n = io)
C--      +19  = signed number of tables of itype
C--      +20  = min index 1 of satellite table
C--      +21  = max index 1 of satellite table
C--      +22  = min id      of satellite table
C--      +23  = max id      of satellite table
C--      +24  = ksat(0)
C--      +25  = ksat(1)   (*i = index 1 satellite table)
C--      +26  = ksat(2)   (*j = id      satellite table)

C     =======================================
      integer function iqcWaddr(wa,i,j,k,l,m)
C     =======================================

C--                        i  j  k  l  m
C--   iwrel = iqcWaddr(wa,iy,it,nf,ig,id)
C--                       13 14 15 16 17
C--
C--   On input, the store wa must be aligned on a table set  
C--   Works on tables of any type
C--   Returns 0 if wa not initialised or if itype = m/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*)

      iqcWaddr = 0
      if(int(wa(1)).ne.123456)          return
      ityp = m/100
      if(ityp.le.0 .or. ityp.gt.mtyp0)  return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                       return

      iqcWaddr = int(wa(iw+12))   +
     +           int(wa(iw+13))*i + int(wa(iw+14))*j +
     +           int(wa(iw+15))*k + int(wa(iw+16))*l +
     +           int(wa(iw+17))*m + int(wa(iw+18))*ioy2
     
      return
      end 

C     =======================================
      integer function iqcWCadr(wa,i,j,k,l,m)
C     =======================================

C--   As iqcWaddr but with array boundary check 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*)

      iqcWCadr = 0
      if(int(wa(1)).ne.123456)                             stop
     +  'iqcWCadr: store not partitioned'
      ityp = m/100
      if(ityp.le.0 .or. ityp.gt.mtyp0)                     stop
     +  'iqcWCadr: impossible table type'
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                                          stop
     +  'iqcWCadr: table type not in store'
      if(i.lt.int(wa(iw   )) .or. i.gt.int(wa(iw+1 ))) stop
     +  'iqcWCadr: index 1 (i) out of range'
      if(j.lt.int(wa(iw+2 )) .or. j.gt.int(wa(iw+3 ))) stop
     +  'iqcWCadr: index 2 (j) out of range'
      if(k.lt.int(wa(iw+4 )) .or. k.gt.int(wa(iw+5 ))) stop
     +  'iqcWCadr: index 3 (k) out of range'
      if(l.lt.int(wa(iw+6 )) .or. l.gt.int(wa(iw+7 ))) stop
     +  'iqcWCadr: index 4 (l) out of range'
      if(m.lt.int(wa(iw+8 )) .or. m.gt.int(wa(iw+9 ))) stop
     +  'iqcWCadr: index 5 (m) out of range'
      iqcWCadr = int(wa(iw+12))   +
     +           int(wa(iw+13))*i + int(wa(iw+14))*j +
     +           int(wa(iw+15))*k + int(wa(iw+16))*l +
     +           int(wa(iw+17))*m + int(wa(iw+18))*ioy2
     
      return
      end 

C     ===================================
      integer function iqcW1ijk(wa,i,j,k)
C     ===================================

C--                        i  j  k
C--   iwrel = iqcW1ijk(wa,iy,ig,id)
C--                       13 16 17
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on tables of type-1
C--   Returns 0 if wa not initialised or if itype = k/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*)

      iqcW1ijk = 0
      if(int(wa(1)).ne.123456)          return
      ityp = k/100
      if(ityp.le.0 .or. ityp.gt.mtyp0)  return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                       return

      iqcW1ijk = int(wa(iw+12))   +
     +           int(wa(iw+13))*i +
     +           int(wa(iw+16))*j +
     +           int(wa(iw+17))*k +
     +           int(wa(iw+18))*ioy2
     
      return
      end

C     ======================================
      integer function iqcW2ijkl(wa,i,j,k,l)
C     ======================================

C--                         i  j  k  l
C--   iwrel = iqcW2ijkl(wa,iy,nf,ig,id)
C--                        13 15 16 17
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on tables of type-2
C--   Returns 0 if wa not initialised or if itype = l/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*)

      iqcW2ijkl = 0
      if(int(wa(1)).ne.123456)           return
      ityp = l/100
      if(ityp.le.0 .or. ityp.gt.mtyp0)   return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                        return

      iqcW2ijkl = int(wa(iw+12))   +
     +            int(wa(iw+13))*i +
     +            int(wa(iw+15))*j +
     +            int(wa(iw+16))*k +
     +            int(wa(iw+17))*l +
     +            int(wa(iw+18))*ioy2
     
      return
      end 
      
C     ======================================
      integer function iqcW3ijkl(wa,i,j,k,l)
C     ======================================

C--                         i  j  k  l
C--   iwrel = iqcW3ijkl(wa,iy,it,ig,id)
C--                        13 14 16 17
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on tables of type-3
C--   Returns 0 if wa not initialised or if itype = l/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*)

      iqcW3ijkl = 0
      if(int(wa(1)).ne.123456)           return
      ityp = l/100
      if(ityp.le.0 .or. ityp.gt.mtyp0)   return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                        return

      iqcW3ijkl = int(wa(iw+12))   +
     +            int(wa(iw+13))*i +
     +            int(wa(iw+14))*j +
     +            int(wa(iw+16))*k +
     +            int(wa(iw+17))*l +
     +            int(wa(iw+18))*ioy2
     
      return
      end       

C     =========================================
      integer function iqcW4ijklm(wa,i,j,k,l,m)
C     =========================================

C--                          i  j  k  l  m
C--   iwrel = iqcW4ijklm(wa,iy,it,nf,ig,id)
C--                         13 14 15 16 17
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on tables of type-4
C--   Returns 0 if wa not initialised or if itype = m/100 does not exist

      implicit double precision (a-h,o-z)

      dimension wa(*)

      iqcW4ijklm = iqcWaddr(wa, i, j, k, l, m)
     
      return
      end 

C     ===================================
      integer function iqcW5ijk(wa,i,j,k)
C     ===================================

C--                        i  j  k
C--   iwrel = iqcW5ijk(wa,iy,it,id)
C--                       13 14 17
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on tables of type-5
C--   Returns 0 if wa not initialised or if itype = k/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*)


      iqcW5ijk  = 0
      if(int(wa(1)).ne.123456)           return
      ityp = k/100
C--   Version without boundary check
*      if(ityp.le.0 .or. ityp.gt.mtyp0)   return
*      npar = int(wa(3))
*      nusr = int(wa(4))
*      iw   = int(wa(ityp+5+npar+nusr))
*      if(iw.eq.0)                        return
C--   Version with boundary check
      if(ityp.ne.5)                                    stop
     +  'iqcW5ijk: not table type 5'
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                                      stop
     +  'iqcW5ijk: table type 5 not in store'
      if(i.lt.int(wa(iw   )) .or. i.gt.int(wa(iw+1 ))) stop
     +  'iqcW5ijk: index 1 (iy) out of range'
      if(j.lt.int(wa(iw+2 )) .or. j.gt.int(wa(iw+3 ))) stop
     +  'iqcW5ijk: index 2 (it) out of range'
      if(k.lt.int(wa(iw+8 )) .or. k.gt.int(wa(iw+9 ))) stop
     +  'iqcW5ijk: index 3 (id) out of range'
C--   End of version with boundary check

      iqcW5ijk  = int(wa(iw+12))   +
     +            int(wa(iw+13))*i +
     +            int(wa(iw+14))*j +
     +            int(wa(iw+17))*k

      return
      end     

C     ================================
      integer function iqcW6ij(wa,i,j)
C     ================================

C--                       i  j
C--   iwrel = iqcW6ij(wa,it,id)
C--                      14 17
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on tables of type-6
C--   Returns 0 if wa not initialised or if itype = j/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*)
      
      iqcW6ij   = 0
      if(int(wa(1)).ne.123456)           return
      ityp = j/100
C--   Version without boundary check
*      if(ityp.le.0 .or. ityp.gt.mtyp0)   return
*      npar = int(wa(3))
*      nusr = int(wa(4))
*      iw   = int(wa(ityp+5+npar+nusr))
*      if(iw.eq.0)                        return
C--   Version with boundary check
      if(ityp.ne.6)                 stop
     +  'iqcW6ij: not table type 6'
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                                      stop
     +  'iqcW6ij: table type 6 not in store'
      if(i.lt.int(wa(iw+2 )) .or. i.gt.int(wa(iw+3 ))) stop
     +  'iqcW6ij: index 1 (it) out of range'
      if(j.lt.int(wa(iw+8 )) .or. j.gt.int(wa(iw+9 ))) stop
     +  'iqcW6ij: index 2 (id) out of range'
C--   End of version with boundary check

      iqcW6ij   = int(wa(iw+12))   +
     +            int(wa(iw+14))*i +
     +            int(wa(iw+17))*j

      return
      end

C     ================================
      integer function iqcW7ij(wa,i,j)
C     ================================

C--                       i  j
C--   iwrel = iqcW7ij(wa,it,id)
C--                      14 17
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on tables of type-7
C--   Returns 0 if wa not initialised or if itype = j/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*)

      iqcW7ij   = 0
      if(int(wa(1)).ne.123456)           return
      ityp = j/100
C--   Version without boundary check
*      if(ityp.le.0 .or. ityp.gt.mtyp0)   return
*      npar = int(wa(3))
*      nusr = int(wa(4))
*      iw   = int(wa(ityp+5+npar+nusr))
*      if(iw.eq.0)                        return
C--   Version with boundary check
      if(ityp.ne.7)                 stop
     +  'iqcW7ij: not table type 7'
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                                      stop
     +  'iqcW7ij: table type 7 not in store'
      if(i.lt.int(wa(iw+2 )) .or. i.gt.int(wa(iw+3 ))) stop
     +  'iqcW7ij: index 1 (it) out of range'
      if(j.lt.int(wa(iw+8 )) .or. j.gt.int(wa(iw+9 ))) stop
     +  'iqcW7ij: index 2 (id) out of range'
C--   End of version with boundary check

      iqcW7ij   = int(wa(iw+12))   +
     +            int(wa(iw+14))*i +
     +            int(wa(iw+17))*j

      return
      end

C     ================================
      integer function iqcWSij(wa,i,j)
C     ================================

C--                       i  j
C--   iwrel = iqcWSij(wa, i,id)
C--                      25 26
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on satellite tables
C--   Returns 0 if wa not initialised or if itype = j/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension wa(*)
      
      iqcWSij   = 0
      if(int(wa(1)).ne.123456)           return
      ityp = j/100
C--   Version without boundary check
*      if(ityp.le.0 .or. ityp.gt.mtyp0)   return
*      npar = int(wa(3))
*      nusr = int(wa(4))
*      iw   = int(wa(ityp+5+npar+nusr))
*      if(iw.eq.0)                        return
C--   Version with boundary check
      if(ityp.le.0 .or. ityp.gt.mtyp0)                 stop
     +  'iqcWSij: wrong table type '
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                                      stop
     +  'iqcWSij: satellite table not in store'
      if(i.lt.int(wa(iw+20)) .or. i.gt.int(wa(iw+21))) stop
     +  'iqcWSij: index 1 (i) out of range'
      if(j.lt.int(wa(iw+22)) .or. j.gt.int(wa(iw+23))) stop
     +  'iqcWSij: index 2 (id) out of range'
C--   End of version with boundary check

      iqcWSij   = int(wa(iw+24))   +
     +            int(wa(iw+25))*i +
     +            int(wa(iw+26))*j

      return
      end

C     =============================
      integer function iqcWSi(wa,i)
C     =============================

C--                      i
C--   iwrel = iqcWSi(wa,id)
C--                     26
C--
C--   On input, the store wa must be aligned on a table set
C--   Works on satellite tables
C--   Returns 0 if wa not initialised or if itype = j/100 does not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension wa(*)
      
      iqcWSi    = 0
      if(int(wa(1)).ne.123456)           return
      ityp = i/100
C--   Version without boundary check
*      if(ityp.le.0 .or. ityp.gt.mtyp0)   return
*      npar = int(wa(3))
*      nusr = int(wa(4))
*      iw   = int(wa(ityp+5+npar+nusr))
*      if(iw.eq.0)                        return
C--   Version with boundary check
      if(ityp.le.0 .or. ityp.gt.mtyp0)                 stop
     +  'iqcWSi: wrong table type'
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                                      stop
     +  'iqcWSi: satellite table not in store'
      if(i.lt.int(wa(iw+22)) .or. i.gt.int(wa(iw+23))) stop
     +  'iqcWSij: index 1 (id) out of range'
C--   End of version with boundary check

      iqcWSi   =  int(wa(iw+24))   +
     +            int(wa(iw+26))*i

      return
      end

C==   ==================================================================
C==   General information ==============================================
C==   ==================================================================

C     ===========================================
      integer function iqcFirstWordOfSet(ww,kset)
C     ===========================================

C--   Return position of first word of table-set kset
C--   If ww is aligned, then return 1
C--
C--   NB: it is not checked that kset is in range

      implicit double precision (a-h,o-z)

      dimension ww(*)

      if(iqcGetNumberOfSets(ww).eq.0)
     + stop 'iqcFirstWordOfSet no table sets in ww'
      if(kset.lt.1 .or. kset.gt.iqcGetNumberOfSets(ww))
     + stop 'iqcFirstWordOfSet wrong kset'

      if(int(ww(1)).eq.654321) then
        mchk              = int(ww(3))
        iqcFirstWordOfSet = int(ww(5+mchk+kset))
      else
        iqcFirstWordOfSet = 1
      endif

      return
      end

C     ========================================
      integer function iqcGetNumberOfWords(ww)
C     ========================================

C--   Get number of words used
C--
C--   If ww not aligned then return the total number of words used
C--   If w is aligned onto kset then return the number of words in kset
C--   If w is not partitioned then the function returns 0
C--
C--       nwtot = iqcGetNumberOfWords(ww)  !total number of words
C--       iwset = iqcFirstWordOfSet(ww,kset)
C--       nwset = iqcGetNumberOfWords(ww(iwset)) !# words in kset

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(int(ww(1)).eq.654321 .or. int(ww(1)).eq.123456) then
        iqcGetNumberOfWords = int(ww(2))
      else
        iqcGetNumberOfWords = 0
      endif

      return
      end

C     =======================================
      integer function iqcGetNumberOfSets(ww)
C     =======================================

C--   Get number of sets in a store
C--
C--   If ww not aligned then return the number of sets
C--   If w is aligned onto kset then return 1
C--   If w is not partitioned then the function returns 0
C--
C--       nsets = iqcGetNumberOfSets(ww)        !total number of sets
C--       iwset = iqcFirstWordOfSet(ww,kset)
C--       nsets = iqcGetNumberOfSets(ww(iwset)) !always 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(int(ww(1)).eq.654321) then
        mchk               = int(ww(3))
        iqcGetNumberOfSets = int(ww(4+mchk))
      elseif(int(ww(1)).eq.123456) then
        iqcGetNumberOfSets = 1
      else
        iqcGetNumberOfSets = 0
      endif

      return
      end

C     ==============================================
      integer function iqcGetNumberOfParams(ww,kset)
C     ==============================================

C--   Return number of evolution parameter words in kset

      implicit double precision (a-h,o-z)

      dimension ww(*)

      if(iqcGetNumberOfSets(ww).eq.0)
     + stop 'iqcGetNumberOfParams no table sets in ww'
      if(kset.lt.1 .or. kset.gt.iqcGetNumberOfSets(ww))
     + stop 'iqcGetNumberOfParams wrong kset'

      ifirst               = iqcFirstWordOfSet(ww,kset)
      iqcGetNumberOfParams = int(ww(ifirst+2))

      return
      end

C     ==============================================
      integer function iqcFirstWordOfParams(ww,kset)
C     ==============================================

C--   Return address of first evolution parameter in kset

      implicit double precision (a-h,o-z)

      dimension ww(*)

      if(iqcGetNumberOfSets(ww).eq.0)
     + stop 'iqcFirstWordOfParams no table sets in ww'
      if(kset.lt.1 .or. kset.gt.iqcGetNumberOfSets(ww))
     + stop 'iqcFirstWordOfParams wrong kset'

      ifirst               = iqcFirstWordOfSet(ww,kset)
      iqcFirstWordOfParams = ifirst+5

      return
      end

C     ==============================================
      integer function iqcGetNumberOfUparam(ww,kset)
C     ==============================================

C--   Return number of user parameter words in kset

      implicit double precision (a-h,o-z)

      dimension ww(*)

      if(iqcGetNumberOfSets(ww).eq.0)
     + stop 'iqcGetNumberOfUparam no table sets in ww'
      if(kset.lt.1 .or. kset.gt.iqcGetNumberOfSets(ww))
     + stop 'iqcGetNumberOfUparam wrong kset'

      ifirst               = iqcFirstWordOfSet(ww,kset)
      iqcGetNumberOfUparam = int(ww(ifirst+3))

      return
      end

C     ==============================================
      integer function iqcFirstWordOfUparam(ww,kset)
C     ==============================================

C--   Return address of first user parameter in kset

      implicit double precision (a-h,o-z)

      dimension ww(*)

      if(iqcGetNumberOfSets(ww).eq.0)
     + stop 'iqcFirstWordOfUParam no table sets in ww'
      if(kset.lt.1 .or. kset.gt.iqcGetNumberOfSets(ww))
     + stop 'iqcFirstWordOfUparam wrong kset'

      ifirst               = iqcFirstWordOfSet(ww,kset)
      npar                 = int(ww(ifirst+2))
      iqcFirstWordOfUparam = ifirst+5+npar

      return
      end


C     ====================================================
      integer function iqcGetNumberOfTables(ww,kset,itype)
C     ====================================================

C--   Get number of tables of type itype in store ww
C--   Undefined if ww not partitioned or if kset,itype does not exist
C--
C--   ww    (in)   store (not aligned)
C--   kset  (in)   table set identifier
C--   itype (in)   table type

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension ww(*)

      if(iqcGetNumberOfSets(ww).eq.0)
     + stop 'iqcGetNumberOfTables no table sets in ww'
      if(kset.lt.1 .or. kset.gt.iqcGetNumberOfSets(ww))
     + stop 'iqcGetNumberOfTables wrong kset'
      if(itype.lt.1 .or. itype.gt.mtyp0)
     + stop 'iqcGetNumberOfTables wrong itype'

      ifirst = iqcFirstWordOfSet(ww,kset)
      iqcGetNumberOfTables = iqcGetNumberOfTabsWa(ww(ifirst),itype)

      return
      end

C     ===============================================
      integer function iqcGetNumberOfTabsWa(wa,itype)
C     ===============================================

C--   Get number of tables of type id/100 in aligned store wa
C--   Returns zero if wa not partitioned or if itype does not exist
C--
C--   wa       (in)   store aligned on a table set
C--   itype    (in)   table type

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension wa(*)

      iqcGetNumberOfTabsWa = 0
      if(int(wa(1)).ne.123456)            return
      if(itype.le.0 .or. itype.gt.mtyp0)  return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(itype+5+npar+nusr))
      if(iw.eq.0)                         return
      idmin = int(wa(iw+8))
      idmax = int(wa(iw+9))
      iqcGetNumberOfTabsWa = idmax-idmin+1

      return
      end

C     ====================================================
      integer function iqcSgnNumberOfTables(ww,kset,itype)
C     ====================================================

C--   Get signed number of tables of type itype in store ww
C--   Undefined if ww not partitioned or if kset,itype does not exist
C--
C--   ww    (in)   store (not aligned)
C--   kset  (in)   table set identifier
C--   itype (in)   table type

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension ww(*)

      if(iqcGetNumberOfSets(ww).eq.0)
     + stop 'iqcSgnNumberOfTables no table sets in ww'
      if(kset.lt.1 .or. kset.gt.iqcGetNumberOfSets(ww))
     + stop 'iqcSgnNumberOfTables wrong kset'
      if(itype.lt.1 .or. itype.gt.mtyp0)
     + stop 'iqcSgnNumberOfTables wrong itype'

      ifirst = iqcFirstWordOfSet(ww,kset)
      iqcSgnNumberOfTables = iqcSgnNumberOfTabsWa(ww(ifirst),itype)

      return
      end

C     ===============================================
      integer function iqcSgnNumberOfTabsWa(wa,itype)
C     ===============================================

C--   Get signed number of tables of type id/100 in aligned store wa
C--   Returns zero if wa not partitioned or if itype does not exist
C--
C--   wa       (in)   store aligned on a table set
C--   itype    (in)   table type

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension wa(*)

      iqcSgnNumberOfTabsWa = 0
      if(int(wa(1)).ne.123456)            return
      if(itype.le.0 .or. itype.gt.mtyp0)  return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(itype+5+npar+nusr))
      if(iw.eq.0)                         return
      iqcSgnNumberOfTabsWa = int(wa(iw+19))
     
      return
      end 

C     =============================================
      subroutine sqcGetLimits(ww,id,imin,imax,jmax)
C     =============================================

C--   Get index limits of table id in store ww
C--   Returns zero imin,imax if ww is not partitioned or id does not exist
C--
C--   ww       (in)   store (not aligned)
C--   id       (in)   identifier in global format can be < 0
C--   imin(6)  (out)  array with lower limits
C--   imax(6)  (out)  array with upper limits
C--   jmax     (out)  upper limit of first index satellite table(j,id)
C--
C--   Indexing in imin, imax :   1  2  3  4  5  6
C--                             iy it nf ig id iosp

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*),imin(*),imax(*)

      if(id.ge.0) then
        kset   =  id/1000
        jd     =  id - 1000*kset
        ifirst =  iqcFirstWordOfSet(ww,kset)
        call sqcGetLimsWa(ww(ifirst),jd,imin,imax,jmax)
      else
        kset   = -id/1000
        jd     = -id - 1000*kset
        ifirst =  iqcFirstWordOfSet(stor7,kset)
        call sqcGetLimsWa(stor7(ifirst),jd,imin,imax,jmax)
      endif

      return
      end

C     =============================================
      subroutine sqcGetLimsWa(wa,id,imin,imax,jmax)
C     =============================================

C--   Get index limits of table id in aligned store wa
C--   Returns zero imin,imax if wa is not partitioned or id does not exist
C--
C--   wa       (in)   store aligned on a table set
C--   id       (in)   local identifier
C--   imin(6)  (out)  array with lower limits
C--   imax(6)  (out)  array with upper limits
C--   jmax     (out)  upper limit of first index of satellite table(j,id)
C--
C--   Indexing in imin, imax :   1  2  3  4  5  6
C--                             iy it nf ig id iosp

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension wa(*),imin(*),imax(*)

      do i = 1,6
        imin(i) = 0
        imax(i) = 0
      enddo
      if(int(wa(1)).ne.123456)          return
      ityp = id/100
      if(ityp.le.0 .or. ityp.gt.mtyp0)  return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                       return

      imin(1) = int(wa(iw   ))          !iy_min
      imax(1) = int(wa(iw+1 ))          !iy_max
      imin(2) = int(wa(iw+2 ))          !it_min
      imax(2) = int(wa(iw+3 ))          !it_max
      imin(3) = int(wa(iw+4 ))          !nf_min
      imax(3) = int(wa(iw+5 ))          !nf_max
      imin(4) = int(wa(iw+6 ))          !ig_min
      imax(4) = int(wa(iw+7 ))          !ig_max
      imin(5) = int(wa(iw+8 ))          !id_min
      imax(5) = int(wa(iw+9 ))          !id_max
      imin(6) = int(wa(iw+10))          !iosp_min
      imax(6) = int(wa(iw+11))          !iosp_max
      jmax    = int(wa(iw+21))          !max of first index satellite
     
      return
      end

C     ==========================================
      subroutine sqcGetLimSpl(ww,id,iosp1,iosp2)
C     ==========================================

C--   Get the (hidden) spline order limits of id in store ww
C--   Returns zero if ww is not partitioned or id does not exist
C--
C--   ww       (in)   store (not aligned)
C--   id       (in)   identifier in global format can be < 0
C--   iosp1    (out)  lower iosp limit
C--   iosp2    (out)  upper iosp limit
C--
C--   Indexing in imin, imax :   1  2  3  4  5  6
C--                             iy it nf ig id iosp

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(id.ge.0) then
        kset   =  id/1000
        jd     =  id - 1000*kset
        ifirst =  iqcFirstWordOfSet(ww,kset)
        call sqcGetLimSpa(ww(ifirst),jd,iosp1,iosp2)
      else
        kset   = -id/1000
        jd     = -id - 1000*kset
        ifirst =  iqcFirstWordOfSet(stor7,kset)
        call sqcGetLimSpa(stor7(ifirst),jd,iosp1,iosp2)
      endif

      return
      end

C     ==========================================
      subroutine sqcGetLimSpa(wa,id,iosp1,iosp2)
C     ==========================================

C--   Get the (hidden) spline order limits of id in aligned store wa
C--   Returns zero if wa is not partitioned or id does not exist
C--
C--   wa       (in)   store aligned on a table set
C--   id       (in)   local identifier
C--   iosp1    (out)  lower iosp limit
C--   iosp2    (out)  upper iosp limit
C--
C--   Indexing in imin, imax :   1  2  3  4  5  6
C--                             iy it nf ig id iosp

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension wa(*)

      iosp1 = 0
      iosp2 = 0

      if(int(wa(1)).ne.123456)          return
      ityp = id/100
      if(ityp.le.0 .or. ityp.gt.mtyp0)  return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                       return

      iosp1 = int(wa(iw+10))          !iosp_min
      iosp2 = int(wa(iw+11))          !iosp_max
     
      return
      end

C     ===================================
      logical function lqcIsDouble(ww,id)
C     ===================================

C--   True if id is a duplicate table

      implicit double precision (a-h,o-z)

      dimension ww(*)

      call sqcGetLimSpl(ww,id,iosp1,iosp2)

      if(iosp1.eq.iosp2) then
        lqcIsDouble = .false.
      else
        lqcIsDouble = .true.
      endif

      return
      end

C     ===========================================
      integer function iqcFirstWordOfTable(ww,id)
C     ===========================================

C--   Return position of first word of table
C--   Undefined if ww is not partitioned or if id does not exist
C--
C--   ww       (in)   store (not aligned)
C--   id       (in)   table id in global format can be < 0

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*),imin(6),imax(6)

      if(id.ge.0) then
        call sqcGetLimits(ww,id,imin,imax,jmax)
        iqcFirstWordOfTable =
     +    iqcGaddr(ww,imin(1),imin(2),imin(3),imin(4),id)
      else
        call sqcGetLimits(stor7,-id,imin,imax,jmax)
        iqcFirstWordOfTable =
     +    iqcGaddr(stor7,imin(1),imin(2),imin(3),imin(4),-id)
      endif

      return
      end

C     ==========================================
      integer function iqcGetTabLeng(ww,id,ndim)
C     ==========================================

C--   Get length of table id in words
C--   Returns zero if ww is not partitioned or if id does not exist
C--
C--   ww       (in)   store (not aligned)
C--   id       (in)   identifier in global format can be < 0
C--   ndim     (in)   number of dimensions iy,it,nf,ig,id,iosp
C--                                         1  2  3  4  5  6

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(ndim.lt.1 .or. ndim.gt.6)
     + stop 'iqcGetTabLeng wrong ndim'

      if(id.ge.0) then
        kset          =  id/1000
        jd            =  id - 1000*kset
        ifirst        =  iqcFirstWordOfSet(ww,kset)
        iqcGetTabLeng =  iqcGetTbLenWa(ww(ifirst),jd,ndim)
      else
        kset          = -id/1000
        jd            = -id - 1000*kset
        ifirst        =  iqcFirstWordOfSet(stor7,kset)
        iqcGetTabLeng =  iqcGetTbLenWa(stor7(ifirst),jd,ndim)
      endif
     
      return
      end 

C     ==========================================
      integer function iqcGetTbLenWa(wa,id,ndim)
C     ==========================================

C--   Get length of table id in words
C--   Returns zero if wa is not partitioned or if id does not exist
C--
C--   wa       (in)   store aligned on a table set
C--   id       (in)   identifier
C--   ndim     (in)   number of dimensions iy,it,nf,ig,id,iosp
C--                                         1  2  3  4  5  6

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension wa(*)

      iqcGetTbLenWa = 0
      if(int(wa(1)).ne.123456)          return
      ityp = id/100
      if(ityp.le.0 .or. ityp.gt.mtyp0)  return
      npar = int(wa(3))
      nusr = int(wa(4))
      iw   = int(wa(ityp+5+npar+nusr))
      if(iw.eq.0)                       return
      idim = max(ndim,1)
      idim = min(ndim,6)
      jdim = 2*idim-1
      leng = 1
      do i = 1,jdim,2
        imin = int(wa(iw+i-1))
        imax = int(wa(iw+i  ))
        leng = leng*(imax-imin+1)
      enddo
      iqcGetTbLenWa = leng
     
      return
      end

C     ====================================
      integer function iqcGetSetNumber(id)
C     ====================================

C--   Get table set numder from global id

      implicit double precision (a-h,o-z)

      iqcGetSetNumber = abs(id)/1000

      return
      end

C     ==================================
      integer function iqcGetLocalId(id)
C     ==================================

C--   Get local id from global id

      implicit double precision (a-h,o-z)

      kset          = abs(id)/1000
      iqcGetLocalId = abs(id) - 1000*kset

      return
      end


C     ====================================
      logical function lqcWPartitioned(ww)
C     ====================================

C--   True if ww partitioned (works also when ww is aligned)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      dimension ww(*)

      if(int(ww(1)).eq.654321 .or. int(ww(1)).eq.123456) then
        lqcwpartitioned = .true.
      else
        lqcwpartitioned = .false.
      endif

      return
      end

C     =======================================
      logical function lqcIsetExists(ww,kset)
C     =======================================

C--   True if kset exists in ww (not aligned)

      implicit double precision (a-h,o-z)

      logical lqcWPartitioned

      dimension ww(*)

      if(.not.lqcWPartitioned(ww)) then
        lqcIsetExists = .false.
      elseif(kset.lt.1.or.kset.gt.iqcGetNumberOfSets(ww)) then
        lqcIsetExists = .false.
      else
        lqcIsetExists = .true.
      endif

      return
      end

C     ==============================================
      logical function lqcItypeExists(ww,kset,itype)
C     ==============================================

C--   True if itype exists in kset of ww (not aligned)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      logical lqcIsetExists

      dimension ww(*)

      if(itype.lt.1 .or. itype.gt.mtyp0)
     + stop 'lqcItypeExists wrong itype'

      if(.not.lqcIsetExists(ww,kset)) then
        lqcItypeExists = .false.
      elseif(iqcGetNumberOfTables(ww,kset,itype).eq.0) then
        lqcItypeExists = .false.
      else
        lqcItypeExists = .true.
      endif

      return
      end

C     ===================================
      logical function lqcIdExists(ww,id)
C     ===================================

C--   True if id (in global format) exists in ww (not aligned)
C--   id can be < 0 to address stor7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      logical lqcIsetExists

      dimension ww(*)

      if(id.ge.0) then
        kset  =  id/1000
        if(.not.lqcIsetExists(ww,kset)) then
          lqcIdExists = .false.
          return
        endif
        idloc   =  id - 1000*kset
        itype   =  idloc/100
        ntable  =  idloc - 100*itype
        ntables =  iqcGetNumberOfTables(ww,kset,itype)
        if(ntable.lt.1 .or. ntable.gt.ntables) then
          lqcIdExists = .false.
        else
          lqcIdExists = .true.
        endif
      else
        kset  = -id/1000
        if(.not.lqcIsetExists(stor7,kset)) then
          lqcIdExists = .false.
          return
        endif
        idloc   = -id - 1000*kset
        itype   =  idloc/100
        ntable  =  idloc - 100*itype
        ntables =  iqcGetNumberOfTables(stor7,kset,itype)
        if(ntable.lt.1 .or. ntable.gt.ntables) then
          lqcIdExists = .false.
        else
          lqcIdExists = .true.
        endif
      endif

      return
      end

C==   ==================================================================
C==   Operations on satellite tables ===================================
C==   ==================================================================

C     =============================
      subroutine sqcValidate(ww,id)
C     =============================

C--   Validate (flag as non-empty) table id (global format) in store ww
C--   Acts as a do-nothing if ww is not partitioned or id does not exist
C--   id can be < 0 to address stor7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      logical lqcIdExists

      dimension ww(*)

      if(id.ge.0) then
        if(.not.lqcIdExists(ww,id))         return
        ia = iqcGSij(ww,1,id)
        ww(ia) = dble(ipver6)
      else
        if(.not.lqcIdExists(stor7,-id))     return
        ia = iqcGSij(stor7,1,-id)
        stor7(ia) = dble(ipver6)
      endif

      return
      end

C     ===============================
      subroutine sqcInvalidate(ww,id)
C     ===============================

C--   Invalidate (flag as empty) table id (global format) in store ww
C--   Acts as a do-nothing if ww is not partitioned or id does not exist
C--   id can be < 0 to address stor7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      logical lqcIdExists

      dimension ww(*)

      if(id.ge.0) then
        if(.not.lqcIdExists(ww,id))         return
        ia = iqcGSij(ww,1,id)
        ww(ia) = 0.D0
      else
        if(.not.lqcIdExists(stor7,-id))     return
        ia = iqcGSij(stor7,1,-id)
        stor7(ia) = 0.D0
      endif

      return
      end

C     ===================================
      logical function lqcIsFilled(ww,id)
C     ===================================

C--   True if id exists and is not invalidated (i.e. id is non-empty)
C--   id in global format can be < 0 to address stor7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      logical lqcIdExists

      dimension ww(*)

      if(id.ge.0) then
        if(lqcIdExists(ww,id)) then
          ia = iqcGSij(ww,1,id)
          if(int(ww(ia)) .eq. 0) then
            lqcIsFilled = .false.
          else
            lqcIsFilled = .true.
          endif
        else
          lqcIsFilled = .false.
        endif
      else
        if(lqcIdExists(stor7,-id)) then
          ia = iqcGSij(stor7,1,-id)
          if(int(stor7(ia)) .eq. 0) then
            lqcIsFilled = .false.
          else
            lqcIsFilled = .true.
          endif
        else
          lqcIsFilled = .false.
        endif
      endif

      return
      end

C     ==================================
      subroutine sqcSetMin6(ww,id,itmin)
C     ==================================

C--   Set itmin cut on type-6 table (id in global format) in store ww
C--   Acts as a do-nothing if ww is not partitioned or id does not exist
C--   or id is not type-6
C--   id can be < 0 to address stor7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      logical lqcIdExists

      dimension ww(*)

      if(int(iqcGetLocalId(id)/100).ne.6)   return

      if(id.ge.0) then
        if(.not.lqcIdExists(ww,id))         return
        ia = iqcGSij(ww,1,id)
        ww(ia+1) = dble(itmin)
      else
        if(.not.lqcIdExists(stor7,-id))     return
        ia = iqcGSij(stor7,1,-id)
        stor7(ia+1) = itmin
      endif

      return
      end

C     ==================================
      subroutine sqcGetMin6(ww,id,itmin)
C     ==================================

C--   Get itmin cut on type-6 table (id in global format) in store ww
C--   Acts as a do-nothing if ww is not partitioned or id does not exist
C--   or id is not type-6
C--   id can be < 0 to address stor7

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'

      logical lqcIdExists

      dimension ww(*)

      itmin = 0

      if(int(iqcGetLocalId(id)/100).ne.6)   return

      if(id.ge.0) then
        if(.not.lqcIdExists(ww,id))         return
        ia = iqcGSij(ww,1,id)
        itmin = int(ww(ia+1))
      else
        if(.not.lqcIdExists(stor7,-id))     return
        ia = iqcGSij(stor7,1,-id)
        itmin = int(stor7(ia+1))
      endif

      return
      end

C==   ==================================================================
C==   Operations on tables =============================================
C==   ==================================================================

C     =======================================
      subroutine sqcTcopyType5(ww,id,it1,it2)
C     =======================================

C--   Copy it1 to it2 in a type5 (pdf) table
C--
C--   ww      (in) store
C--   id      (in) global identifier of type5 table
C--   it1,it2 (in) copy bin it1 to it2

      implicit double precision (a-h,o-z)

      dimension ww(*), mi(6), ma(6)

      logical lqcIdExists

      if(.not.lqcIdExists(ww,id))
     +    stop 'sqcTcopyType5: id does not exist'
      call sqcGetLimits(ww,id,mi,ma,jmax)
      ia1 = iqcG5ijk(ww,mi(1),it1,id)-1
      ia2 = iqcG5ijk(ww,mi(1),it2,id)-1
      do iy = mi(1),ma(1)
        ia1     = ia1+1
        ia2     = ia2+1
        ww(ia2) = ww(ia1)
      enddo

      return
      end

C     ======================================
      subroutine sqcCopyType5(w1,id1,w2,id2)
C     ======================================

C--   Copy type5 table (+ satellite) id1 in w1 to id2 in w2
C--
C--   w1   (in) : workspace 1
C--   id1  (in) : type5 table identifier in w1 (global format)
C--   w2   (in) : workspace 2
C--   id2  (in) : type5 table identifier in w2 (global format)

      implicit double precision (a-h,o-z)

      dimension w1(*), w2(*), mi(6), ma(6)

      logical lqcIdExists

      if(.not.lqcIdExists(w1,id1)) stop 'sqcCopyType5: nonexisting id1'
      if(.not.lqcIdExists(w2,id2)) stop 'sqcCopyType5: nonexisting id2'
      call sqcGetLimits(w1,id1,mi,ma,jmax)
C--   Copy pdf table
      ia1 = iqcG5ijk(w1,mi(1),mi(2),id1)
      ia2 = iqcG5ijk(w1,ma(1),ma(2),id1)
      ibb = iqcG5ijk(w2,mi(1),mi(2),id2)-1
      do ia = ia1,ia2
        ibb     = ibb+1
        w2(ibb) = w1(ia)
      enddo
C--   Copy satellite table
      ia1 = iqcGSij(w1,1,id1)-1
      ia2 = iqcGSij(w2,1,id2)-1
C--   Copy jmax words
      do i = 1,jmax
        ia1     = ia1+1
        ia2     = ia2+1
        w2(ia2) = w1(ia1)
      enddo

      return
      end

C     ======================================
      subroutine sqcCopyType6(w1,id1,w2,id2)
C     ======================================

C--   Copy type6 table (+satellite) id1 in w1 to id2 in w2
C--
C--   w1   (in) : workspace 1
C--   id1  (in) : type6 table identifier in w1 (global format)
C--   w2   (in) : workspace 2
C--   id2  (in) : type6 table identifier in w2 (global format)

      implicit double precision (a-h,o-z)

      dimension w1(*), w2(*), mi(6), ma(6)

      logical lqcIdExists

      if(.not.lqcIdExists(w1,id1)) stop 'sqcCopyType6: nonexisting id1'
      if(.not.lqcIdExists(w2,id2)) stop 'sqcCopyType6: nonexisting id2'
      call sqcGetLimits(w1,id1,mi,ma,jmax)
C--   Copy alfas table
      ia1 = iqcG6ij(w1,mi(2),id1)
      ia2 = iqcG6ij(w1,ma(2),id1)
      ibb = iqcG6ij(w2,mi(2),id2)-1
      do ia = ia1,ia2
        ibb     = ibb+1
        w2(ibb) = w1(ia)
      enddo
C--   Copy satellite table
      ia1 = iqcGSij(w1,1,id1)-1
      ia2 = iqcGSij(w2,1,id2)-1
C--   Copy jmax words
      do i = 1,jmax
        ia1     = ia1+1
        ia2     = ia2+1
        w2(ia2) = w1(ia1)
      enddo

      return
      end

C     ======================================
      subroutine sqcCopyType7(w1,id1,w2,id2)
C     ======================================

C--   Copy type7 table (+satellite) id1 in w1 to id2 in w2
C--
C--   w1   (in) : workspace 1
C--   id1  (in) : type7 table identifier in w1 (global format)
C--   w2   (in) : workspace 2
C--   id2  (in) : type7 table identifier in w2 (global format)

      implicit double precision (a-h,o-z)

      dimension w1(*), w2(*), mi(6), ma(6)

      logical lqcIdExists

      if(.not.lqcIdExists(w1,id1)) stop 'sqcCopyType7: nonexisting id1'
      if(.not.lqcIdExists(w2,id2)) stop 'sqcCopyType7: nonexisting id2'
      call sqcGetLimits(w1,id1,mi,ma,jmax)
C--   Copy pointer table
      ia1 = iqcG7ij(w1,mi(2),id1)
      ia2 = iqcG7ij(w1,ma(2),id1)
      ibb = iqcG7ij(w2,mi(2),id2)-1
      do ia = ia1,ia2
        ibb     = ibb+1
        w2(ibb) = w1(ia)
      enddo
C--   Copy satellite table
      ia1 = iqcGSi(w1,id1)-1
      ia2 = iqcGSi(w2,id2)-1
C--   Copy jmax words
      do i = 1,jmax
        ia1     = ia1+1
        ia2     = ia2+1
        w2(ia2) = w1(ia1)
      enddo

      return
      end

C==   ==================================================================
C==   Book tables ======================================================
C==   ==================================================================

C     =============================================================
      subroutine sqcMakeTab(ww,nw,itypes,npar,nusr,new,kset,nwords)
C     =============================================================

C--   Add a set of tables to the store ww
C--
C--   ww        (in)  array dimensioned nw in the calling routine
C--   nw        (in)  dimension of ww
C--   itypes(7) (in)  number of tables of type (i) >0 splitting function
C--                                                =0 no tables type (i)
C--                                                <0 coeff function
C--   npar      (in)  size of evolution parameter space
C--   nusr      (in)  size of user parameter space
C--   new       (in)  = 1 overwrite existing sets, otherwise add new set
C--   kset      (out) set identifier <= mst0
C--                    -1 attempt to create set with 0 tables
C--                    -2 not enough space
C--                    -3 iset limit mst0 exceeded
C--   nwords    (out) total number of words used < 0 if not enough space

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'                 !needed for spline order ioy2

      save iwversion              !iwversion stamps different workspaces
      data iwversion/0/

!Below is an opemMP directive to make the routine threadsafe for xFitter
!$OMP THREADPRIVATE(iwversion)

      dimension ww(*), itypes(mtyp0)

      if(int(ww(1)).ne.654321 .or. new .eq.1) then
C--     Initialize
        do i = 1,nw
          ww(i) = 0.D0
        enddo
        kset      = 1
        ifirstw   = 7 + mchk0 + mst0
        iwversion = iwversion+1
      elseif(int(ww(4+mchk0)).lt.mst0) then
C--     iset count not yet exceeded
        kset    = int(ww(4+mchk0))+1
        ifirstw = int(ww(5+mchk0+kset))
      else
C--     iset count exceeded
        kset   = -3
        nwords =  0
        return
      endif

C--   Check enough space
      if(ifirstw.gt.nw) then
        nwords = -ifirstw+1
        kset   = -2
        return
      endif

C--   Go ...
      nwleft = nw-ifirstw+1
      call sqcBookSet(ww(ifirstw),nwleft,itypes,npar,nusr,nwu,ierr)

      if(ierr.eq.1) then
C--     nothing to book
        kset   = -1
        nwords = -ifirstw+1
        return
      elseif(ierr.eq.2) then
C--     not enough space
        kset   = -2
        nwords = -ifirstw-nwu+1
        return
      endif

C--   All OK now set all the header words
      nwords = ifirstw+nwu-1
      ww(1)  = dble(654321)
      ww(2)  = dble(nwords)
      ww(3)  = mchk0
      ww(4)            = dble(iwversion)  !version# in 1st checksum word
      ww(4+mchk0)      = dble(kset)
      ww(5+mchk0)      = dble(mst0)
      ww(5+mchk0+kset) = dble(ifirstw)
      ww(6+mchk0+kset) = dble(nwords+1)

C--   Flag all tables as unfilled
      do itype = 1,mtyp0
        ntab = iqcGetNumberOfTables(ww,kset,itype)
        do i = 1,ntab
          id = 1000*kset + 100*itype + i
          call sqcInvalidate(ww,id)
        enddo
      enddo

      return
      end

C     =========================================================
      subroutine sqcBookSet(wa,nw,itypes,mpar,musr,nwords,ierr)
C     =========================================================

C--   Create a set of tables with w(1) the first word of the set
C--
C--   wa        (in)  store, aligned on first free word
C--   nw        (in)  number of words available in wa
C--   itypes(6) (in)  number of tables of type (i) >0 splitting function
C--                                                =0 no tables type (i)
C--                                                <0 coeff function
C--   mpar      (in)  size of evolution parameter space
C--   musr      (in)  size of user parameter space
C--   nwords    (out) number of words used/required
C--   ierr      (out) 0 all OK
C--                   1 no valid itypes encountered (empty set)
C--                   2 not enough space               

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'    !needed for spline order ioy2

      dimension wa(*), itypes(mtyp0)
      dimension imin(6),imax(6),karr(0:6)
      dimension jmin(2),jmax(2),ksat(0:2)

C--   Initialize store
      npar = max(mpar,1)
      nusr = max(musr,1)
      do i = 1,nw
        wa(i) = 0.D0
      enddo
      do i = 1,6
        imin(i) = 0
        imax(i) = 0
        karr(i) = 0
      enddo
      karr(0) = 0
      do i = 1,2
        jmin(i) = 0
        jmax(i) = 0
        ksat(i) = 0
      enddo
      ksat(0) = 0
      ifirstw = 6 + npar + nusr + mtyp0
      nwheadr = 20  !6 imin + 6 imax + 7 karr + 1 ntables
      nwheadr = nwheadr + 7  !2 jmin + 2 jmax + 3 ksat for sattelites
      ierr    = 0
      jwrite  = 0
C--   Loop over table types
      do ityp = 1,mtyp0
        iwrite = 0
        if(ityp.eq.1 .and. 
     +    abs(itypes(1)).ge.1 .and. abs(itypes(1)).le.99) then
C--       Type-1 table (6-dim)
          imin(1) =  1
          imax(1) =  nyy2(0)
          imin(2) =  1
          imax(2) =  1
          imin(3) =  3
          imax(3) =  3
          imin(4) =  1
          imax(4) =  nyg2
          imin(5) =  101
          imax(5) =  100+abs(itypes(1))
          if(itypes(1).gt.0) then
C--         Splitting function table
            imin(6) = 2
            imax(6) = ioy2
          else
C--         Coefficient function table
            imin(6) = ioy2
            imax(6) = ioy2
          endif
          ntables =  itypes(1)
C--       Type-1 satellite table  (2-dim)
          jmin(1) =  1
          jmax(1) =  1
          jmin(2) =  imin(5)
          jmax(2) =  imax(5)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.2 .and.
     +    abs(itypes(2)).ge.1 .and. abs(itypes(2)).le.99) then
C--       Type-2 table (6-dim)
          imin(1) =  1
          imax(1) =  nyy2(0)
          imin(2) =  1
          imax(2) =  1
          imin(3) =  3
          imax(3) =  6
          imin(4) =  1
          imax(4) =  nyg2
          imin(5) =  201
          imax(5) =  200+abs(itypes(2))
          if(itypes(2).gt.0) then
C--         Splitting function table
            imin(6) = 2
            imax(6) = ioy2
          else
C--         Coefficient function table
            imin(6) = ioy2
            imax(6) = ioy2
          endif
          ntables =  itypes(2)
C--       Type-2 satellite table  (2-dim)
          jmin(1) =  1
          jmax(1) =  1
          jmin(2) =  imin(5)
          jmax(2) =  imax(5)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.3 .and.
     +    abs(itypes(3)).ge.1 .and. abs(itypes(3)).le.99) then
C--       Type-3 table (6-dim)
          imin(1) =  1
          imax(1) =  nyy2(0)
          imin(2) =  1
          imax(2) =  ntt2
          imin(3) =  3
          imax(3) =  3
          imin(4) =  1
          imax(4) =  nyg2
          imin(5) =  301
          imax(5) =  300+abs(itypes(3))
          if(itypes(3).gt.0) then
C--         Splitting function table
            imin(6) = 2
            imax(6) = ioy2
          else
C--         Coefficient function table
            imin(6) = ioy2
            imax(6) = ioy2
          endif
          ntables =  itypes(3)
C--       Type-3 satellite table (2-dim)
          jmin(1) =  1
          jmax(1) =  1
          jmin(2) =  imin(5)
          jmax(2) =  imax(5)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.4 .and.
     +    abs(itypes(4)).ge.1 .and. abs(itypes(4)).le.99) then
C--       Type-4 table (6-dim)
          imin(1) =  1
          imax(1) =  nyy2(0)
          imin(2) =  1
          imax(2) =  ntt2
          imin(3) =  3
          imax(3) =  6
          imin(4) =  1
          imax(4) =  nyg2
          imin(5) =  401
          imax(5) =  400+abs(itypes(4))
          if(itypes(4).gt.0) then
C--         Splitting function table
            imin(6) = 2
            imax(6) = ioy2
          else
C--         Coefficient function table
            imin(6) = ioy2
            imax(6) = ioy2
          endif
          ntables =  itypes(4)
C--       Type-4 satellite table (2-dim)
          jmin(1) =  1
          jmax(1) =  1
          jmin(2) =  imin(5)
          jmax(2) =  imax(5)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.5 .and.
     +    abs(itypes(5)).ge.1 .and. abs(itypes(5)).le.99) then
C--       Type-5 table (6-dim)
          imin(1) =  0
          imax(1) =  nstory2
          imin(2) = -mxg0-2            !store subgrid start/endvalues
          imax(2) =  ntt2+mnf0         !enough space to hold the z-grid
          imin(3) =  3
          imax(3) =  3
          imin(4) =  1
          imax(4) =  1
          imin(5) =  501
          imax(5) =  500+abs(itypes(5))
          imin(6) =  ioy2
          imax(6) =  ioy2
          ntables =  itypes(5)
C--       Type-5 satellite table (2-dim)
          jmin(1) =  1
          jmax(1) =  5
          jmin(2) =  imin(5)
          jmax(2) =  imax(5)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.6 .and.
     +    abs(itypes(6)).ge.1 .and. abs(itypes(6)).le.99) then
C--       Type-6 table (6-dim)
          imin(1) =  1
          imax(1) =  1
          imin(2) =  1
          imax(2) =  ntt2+mnf0
          imin(3) =  3
          imax(3) =  3
          imin(4) =  1
          imax(4) =  1
          imin(5) =  601
          imax(5) =  600+abs(itypes(6))
          imin(6) =  ioy2
          imax(6) =  ioy2
          ntables =  itypes(6)
C--       Type-6 satellite table (2-dim)
          jmin(1) =  1
          jmax(1) =  2
          jmin(2) =  imin(5)
          jmax(2) =  imax(5)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.7 .and.
     +    abs(itypes(7)).ge.1 .and. abs(itypes(7)).le.99) then
C--       Type-7 table (6-dim)
          imin(1) =  1
          imax(1) =  1
          imin(2) = -ntt2-mnf0         !enough space to hold the z-grid
          imax(2) =  ntt2+mnf0         !enough space to hold the z-grid
          imin(3) =  3
          imax(3) =  3
          imin(4) =  1
          imax(4) =  1
          imin(5) =  701
          imax(5) =  700+abs(itypes(7))
          imin(6) =  ioy2
          imax(6) =  ioy2
          ntables =  itypes(7)
C--       Type-7 satellite table (2-dim)
          jmin(1) =  1
          jmax(1) =  1
          jmin(2) =  imin(5)
          jmax(2) =  imax(5)
          iwrite  =  1
          jwrite  =  1
        endif
        if(iwrite.eq.1) then
C--       Table partition
          call smb_bkmat(imin,imax,karr,6,ifirstw+nwheadr,iend)
C--       Satellite table
          call smb_bkmat(jmin,jmax,ksat,2,iend+1,ilastw)
          if(ilastw.le.nw) then
C--         Store table partition definition (if enough space)
            iw = ifirstw-1
            do i = 1,6
              iw     = iw+1
              wa(iw) = dble(imin(i))
              iw     = iw+1
              wa(iw) = dble(imax(i))
            enddo
            do i = 0,6
              iw     = iw+1
              wa(iw) = dble(karr(i))
            enddo
            iw = iw+1
            wa(iw) = dble(ntables)
C--         Partition definition of satellite table
            do i = 1,2
              iw     = iw+1
              wa(iw) = dble(jmin(i))
              iw     = iw+1
              wa(iw) = dble(jmax(i))
            enddo
            do i = 0,2
              iw     = iw+1
              wa(iw) = dble(ksat(i))
            enddo
C--         Store base address
            wa(ityp+5+npar+nusr) = dble(ifirstw)
          endif
C--       Start of next set of tables
          ifirstw = ilastw+1
        endif
C--   End of loop over table types
      enddo

C--   No table type found
      if(jwrite.eq.0) then
         ierr   =  1
         nwords =  0
C--   Enough words? 
      elseif(ilastw.le.nw) then
        nwords = ilastw
        wa(1)  = dble(123456)
        wa(2)  = dble(nwords)
        wa(3)  = dble(npar)
        wa(4)  = dble(nusr)
        wa(5)  = dble(mtyp0)
      else
        nwords = ilastw
        ierr   = 2
      endif

      return
      end

