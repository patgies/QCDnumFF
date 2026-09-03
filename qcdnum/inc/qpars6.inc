      logical Lscopechek6, Lqswrite6
      common /lpars6/ Lscopechek6, Lqswrite6
      common /qpard6/ qthrs6(4:6),tthrs6(4:6),rthrs6(4:6),
     &                q0alf6, alfq06, aar6, bbr6, aslim6
      common /qpari6/ nfmin6, nfmax6, itmin6, iord6, nfix6, ievtyp6
      common /qthrs6/ ntsubg6, nfsubg6(4),it1sub6(4),it2sub6(4),
     +                iz1sub6(4),iz2sub6(4),itlim6(10),izlim6(13),
     +                izminf6(3:6),izmanf6(3:6)
      common /qvers6/ ipver6, itver6
      common /qinfty/ qlimd6, qlimu6
      common /steer6/ niter6, itlmc6, idbug6
      common /epsval/ aepsi6, dflim6, qnull6, gepsi6
      common /evcomm/ nnopt6(0:mord0), inopt6
      common /scope6/ iscopeslot6, iscopeuse6, iscopekey6
      common /plist6/ parlist6(mpar0,0:mpl0),keycount6(mpl0)
      common /plifo6/ parbf6(mpar0,npbuf0)
      common /qstor6/ qstore6(mqs0)
      common /pbits6/ ipbits6

C--   Flagbits in ipbits6
      parameter( infbit6 = 1 )
      parameter( iasbit6 = 2 )
      parameter( izcbit6 = 3 )
      parameter( ipsbit6 = 4 )

C--   # parameters used for checking if thresholds are equal
      parameter( npcthr6 = 10)
C--   # parameters used for checking if parameters are equal
      parameter( npcomp6 = 18)
C--   # parameters copied to pdf sets
      parameter( npcopy6 = 35 )

      parameter( idnfix6   =  1 )
      parameter( idqthrc6  =  2 )
      parameter( idqthrb6  =  3 )
      parameter( idqthrt6  =  4 )
      parameter( idtthrc6  =  5 )
      parameter( idtthrb6  =  6 )
      parameter( idtthrt6  =  7 )
      parameter( idrthrc6  =  8 )
      parameter( idrthrb6  =  9 )
      parameter( idrthrt6  = 10 )
C--   Up to here thresholds are checked for equality
      parameter( idiymac2  = 11 )
      parameter( iditmic2  = 12 )
      parameter( iditmac2  = 13 )
      parameter( idiord6   = 14 )
      parameter( idalfq06  = 15 )
      parameter( idq0alf6  = 16 )
      parameter( idaar6    = 17 )
      parameter( idbbr6    = 18 )
C--   Up to here  parameters are checked for equality
      parameter( idaslim6  = 19 )
      parameter( idipver6  = 20 )
      parameter( iditver6  = 21 )
      parameter( idievtyp6 = 22 )
      parameter( idnfmin6  = 23 )
      parameter( idnfmax6  = 24 )
      parameter( iditmin6  = 25 )
      parameter( idizmic2  = 26 )
      parameter( idizmac2  = 27 )
      parameter( idz136    = 28 )
      parameter( idz146    = 29 )
      parameter( idz156    = 30 )
      parameter( idz166    = 31 )
      parameter( idz236    = 32 )
      parameter( idz246    = 33 )
      parameter( idz256    = 34 )
      parameter( idz266    = 35 )
C--   Up to here parameters are copied to pdf sets

