
C--   This you can set as needed ---------------------------------------

      parameter( mebuf0 = 50 )               !Size of message buffer

C--   Do not touch anything below --------------------------------------
C--
C--   Address notation:  IT                   Table coordinates
C--                      IS = IS0 + IT        Table-set coordinates
C--                      IW = IW0 + IS        Workspace coordinates
C--                      IA = IA0 + IW        Array indices (IA0 = root)

C--   General ----------------------------------------------------------
      parameter( iroot0       = 1          )
      parameter( iWSversion0  = 20211117   )
      parameter( iCWorkspace0 = 920210714  )
      parameter( iCWTableSet0 = 987654321  )
      parameter( iCWTable0    = 123456789  )

C--   Workspace --------------------------------------------------------
      parameter( mdim0        = 25         )
      parameter( nwHeader0    = 15         )
C--   All headers -- never redefine the numbers 0 -- 5
      parameter(ICword0 =  0 )       !control word
      parameter(IwAddr0 =  1 )       !IW of object (distance to root)
      parameter(NFTabl0 =  2 )       !distance to next table marker
      parameter(NBTabl0 =  3 )       !distance to previous table marker
      parameter(NFTset0 =  4 )       !distance to next tset marker
      parameter(NBTset0 =  5 )       !distance to previous tset marker
      parameter(IFprnt0 =  6 )       !fingerprint
      parameter(NObjec0 =  7 )       !number of objects
      parameter(IObjec0 =  8 )       !object number
      parameter(INwUse0 =  9 )       !Words used
C--   Workspace header
      parameter(IWlstS0 = 10 )       !IW last table-set
      parameter(IWlstT0 = 11 )       !IW last table
      parameter(INwMax0 = 12 )       !Wspace total size
      parameter(IDrain0 = 13 )       !Drain word
      parameter(IZnull0 = 14 )       !Null word
C--   Table-set header
      parameter(INhead0 = 10 )       !Header size
      parameter(INtags0 = 11 )       !Tag field size
      parameter(IHskip0 = 12 )       !Hskip
      parameter(ISlstT0 = 13 )       !IS last table
C--   Table header
      parameter(ItMeta0 = 10 )       !IT metadata 1st word
      parameter(ItImin0 = 11 )       !IT IMIN array
      parameter(ItImax0 = 12 )       !IT IMAX array
      parameter(ItBody0 = 13 )       !IT begin table body
      parameter(ItBend0 = 14 )       !IT end table body

C--   Istore -----------------------------------------------------------
      parameter( nwHeader1    = 12         )
C--   All headers -- never redefine the numbers 0 -- 3
      parameter(ICword1 =  0 )       !control word
      parameter(IwAddr1 =  1 )       !IW of object (distance to root)
      parameter(NFTabl1 =  2 )       !distance to next table marker
      parameter(NBTabl1 =  3 )       !distance to previous table marker
      parameter(IFprnt1 =  4 )       !fingerprint
      parameter(INwUse1 =  5 )       !Words used
C--   Istore header
      parameter(NObjec1 =  6 )       !number of arrays
      parameter(IWlstT1 =  7 )       !IW last array
      parameter(INwMax1 =  8 )       !Istore total size
      parameter(INhead1 =  9 )       !Header size
C--   Array header
      parameter(IObjec1 =  6 )       !Array number
      parameter(KnulTb1 =  7 )       !K0 pointer coefficient
      parameter(IminTb1 =  8 )       !Imin value
      parameter(ImaxTb1 =  9 )       !Imax value
      parameter(ItBody1 = 10 )       !IT begin array body
      parameter(ItBend1 = 11 )       !IT end array body

