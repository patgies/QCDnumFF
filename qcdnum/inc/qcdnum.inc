
C--   --------- you can adjust these to your needs ---------------------

      parameter ( mxg0  = 5       )  !max # x-subgrids
      parameter ( mxx0  = 320     )  !max # x-gridpoints
      parameter ( mqq0  = 170     )  !max # Q2-gridpoints
      parameter ( mst0  = 30      )  !max # sets in a workspace .ge. 10
      parameter ( mce0  = 20      )  !max # coupled evolutions
      parameter ( mbf0  = 10      )  !max # fast engine scratch buffers
      parameter ( mky0  = 50      )  !max # of keycards must be .ge. 13
      parameter ( mqs0  = 500     )  !tot # words in qstore
      parameter ( nwf0  = 2000000 )  !tot # of words in pdf store

C--   --------- adjust only when error message -------------------------

      parameter ( msl0  = 100     )  !max C++ input string length

C--   --------- do not touch the values below --------------------------

      parameter ( mord0  = 3      ) !max perturbative order
      parameter ( iom0   = 3      ) !max spline order in x or Q2
      parameter ( mnf0   = 12     ) !# q + qbar flavors
      parameter ( mbp0   = 3      ) !# of status words
      parameter ( nwrk0  = 5      ) !# scratch pdf tables (..,-3,-2,-1)
      parameter ( mpp0   = 7      ) !# Pij tables (QQ,QG,GQ,GG,PL,MI,VA)
      parameter ( maa0   = 3      ) !matching matrix dimension
      parameter ( mpdf0  = 90     ) !max # pdfs in a pdf set
      parameter ( mset0  = mst0-6 ) !max # pdf sets in memory
      parameter ( mstp0  = 4      ) !# Pij types unp/pol/frag/custom
      parameter ( mtyp0  = 7      ) !# of table types
      parameter ( mchk0  = 3      ) !# of checksums must be .ge. 3
      parameter ( iwg0   = 5      ) !store gridversion in ww(iwg0)
      parameter ( mpar0  = 45     ) !# of evolution parameters
      parameter ( mgrp0  = 6      ) !# of parameter groups
      parameter ( mpl0   = mst0   ) !# of parameter lists
      parameter ( mpt0   = 5000   ) !max # interpolations in one call
      parameter ( npbuf0 = 5      ) !depth of param LIFO buffer .le. 6
      parameter ( nsubt0 = 2      ) !# of t-subgrid pointer tables

C--   --------- size of parameter store (generous estimate) ----------

      parameter (
     +  nwp0 = (mqq0*(2*nsubt0+2*mord0+2)+mpar0+100)*mpl0+100
     +          )
