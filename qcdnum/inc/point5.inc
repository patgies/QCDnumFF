
C--   Pointers in the current pdf set
      common /qpntr5/ ksetw5
C--   iz15(6) : 1=nfmin, 2=izmicut, 3-6=iz1(nf)
C--   iz25(6) : 1=nfmax, 2=izmacut, 3-6=iz2(nf)
      common /ilims5/ nzz5, iz15(6), iz25(6), iymac5, ixmic5,
     +                itmic5, itmac5, izmic5, izmac5, itmin5
      common /dlims5/ xmic5, ymac5, tmic5, tmac5, qmic5, qmac5
      common /qvers5/ ipver5 , itver5
      common /qdpar5/ aar5, bbr5
      common /qipar5/ iord5, nfix5
      common /qsubg5/ izfit5(-mqq0-7:mqq0+7), itfiz5(-mqq0-7:mqq0+7)
