      common /pstor8/ pars8(nwp0)
      common /pcntr8/ keycnt8(mpl0,mgrp0), keyadr8(mpl0)
      common /pbits8/ ipbits8

C--   Flagbits in ipbits8
      parameter( infbit8 = 1 )
      parameter( iasbit8 = 2 )
      parameter( izcbit8 = 3 )
      parameter( ipsbit8 = 4 )

C--   # parameters copied to pdf sets
      parameter( npcopy8 = 41 )

C--   Order
      parameter( idiord8    =  1 )
C--   Alfas
      parameter( idalfq08   =  2 )
      parameter( idq0alf8   =  3 )
C--   Thresholds
      parameter( idnfix8    =  4 )
      parameter( idqthrc8   =  5 )
      parameter( idqthrb8   =  6 )
      parameter( idqthrt8   =  7 )
      parameter( idtthrc8   =  8 )
      parameter( idtthrb8   =  9 )
      parameter( idtthrt8   = 10 )
      parameter( idrthrc8   = 11 )
      parameter( idrthrb8   = 12 )
      parameter( idrthrt8   = 13 )
C--   Scale
      parameter( idaar8     = 14 )
      parameter( idbbr8     = 15 )
C--   Cuts
      parameter( idiymac8   = 16 )
      parameter( iditmic8   = 17 )
      parameter( iditmac8   = 18 )
C--   Up to here  parameters are checked for equality
      parameter( idaslim8   = 19 )
      parameter( idievtyp8  = 20 )
      parameter( idnfheavy8 = 21 )
      parameter( idnfmin8   = 22 )
      parameter( idnfmax8   = 23 )
      parameter( iditmin8   = 24 )
      parameter( idizmic8   = 25 )
      parameter( idizmac8   = 26 )
      parameter( idnzz8     = 27 )
      parameter( idz138     = 28 )
      parameter( idz148     = 29 )
      parameter( idz158     = 30 )
      parameter( idz168     = 31 )
      parameter( idz238     = 32 )
      parameter( idz248     = 33 )
      parameter( idz258     = 34 )
      parameter( idz268     = 35 )
C--   Version numbers
      parameter( idiover8   = 36 )
      parameter( idiaver8   = 37 )
      parameter( iditver8   = 38 )
      parameter( idisver8   = 39 )
      parameter( idicver8   = 40 )
      parameter( idipver8   = 41 )
C--   Up to here parameters are copied to pdf sets

C--   Parameter group start index
      parameter (iparo8 =  1 )                                    !order
      parameter (ipara8 =  2 )                                     !alfa
      parameter (ipart8 =  4 )                               !thresholds
      parameter (ipars8 = 14 )                                    !scale
      parameter (iparc8 = 16 )                                     !cuts
      parameter (iparp8 =  1 )                                      !all

C--   Parameter group number of items
      parameter (nparo8 =  1 )                                    !order
      parameter (npara8 =  2 )                                     !alfa
      parameter (npart8 = 10 )                               !thresholds
      parameter (npars8 =  2 )                                    !scale
      parameter (nparc8 =  3 )                                     !cuts
      parameter (nparp8 = 18 )                                      !all


