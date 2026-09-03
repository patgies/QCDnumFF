
C--   Adjust these parameters to your need -----------------------------

      parameter ( nw0   = 900000 )      !internal workspace size
      parameter ( maxn0 =   1000 )      !max number of nodes
      parameter ( deps0 =  1.D-9 )      !comparison tolerance
      parameter ( nusr0 =    100 )      !in-spline user store size
      
C--   Leave these alone ------------------------------------------------

C--   Release date (SPLINT code version number)
      parameter ( IvDate0   = 20220308  )
      
C--   File key (SPLINT file format version number)
      parameter ( KeyVal0   = 20210919  )

C--   Spline marker
      parameter ( MarkSp0   = 185218521 )
      
C--   Workspace tag size
      parameter ( Ntags0    = 18 )   !tag field size
      
C--   Workspace tag field indices
      parameter ( Iversion0 = 0 )    !Splint version
      parameter ( IaUbegin0 = 1 )    !IA of user store
      parameter ( NwUstore0 = 2 )    !size of user store
      parameter ( IaSfirst0 = 3 )    !IA of first spline
      
C--   Spline-set tag field indices
      parameter ( ImarkS0  =  0 )    !spline marker
      parameter ( IrOnly0  =  1 )    !Readonly flag
      parameter ( IdimSp0  =  2 )    !spline dimension (+-1,2)
      parameter ( IrsCut0  =  3 )    !roots cut
      parameter ( NedegU0  =  4 )    !degree extrap u
      parameter ( NedegV0  =  5 )    !degree extrap v
      parameter ( IsUtab0  =  6 )    !IS of u-nodes
      parameter ( NwUtab0  =  7 )    !number of u-nodes
      parameter ( IsVtab0  =  8 )    !IS of v-nodes
      parameter ( NwVtab0  =  9 )    !number of v-nodes
      parameter ( Nodesa0  = 10 )    !number of active nodes
      parameter ( IsFtab0  = 11 )    !IS 1-d tab of F coefficients
      parameter ( IsBtab0  = 12 )    !IS 1-d tab of B coefficients
      parameter ( IsCtab0  = 13 )    !IS 1-d tab of C coefficients
      parameter ( IsDtab0  = 14 )    !IS 1-d tab of D coefficients
      parameter ( IsUser0  = 15 )    !IS first word of user table
      parameter ( IsFBCD0  = 16 )    !IS of 2-d coefficient table Fij
      parameter ( IsCCij0  = 17 )    !IS of 2-d coefficient copy Cij



                        
