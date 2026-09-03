
C--   15-05-12   v01-01 - Allow MFNS when icbt < 0
C--   10-04-12   v01-00 - Set nhqstor = 6000000    in hqstf.inc
C--                     - Remove nhqpnts parameter in hqstf.inc
C--                     - Arbitray long list of interpolation points
C--                       allowed in HQSTFUN.
C--   19-07-11   v00-04 - Force Q2 > 0.25 in coefficient functions 
C--                     - Set nhqstor = 900000 in hqstf.inc 
C--   30-03-11   v00-03 - Print progress in weight calculation
C--   16-04-10   v00-02 - Introduce HSWITCH routine
C--   10-12-09   v00-01 - Introduce Q2 limits
C--   21-07-09   v00-00 - Initial version 

C--   HQSTF version stamp
      character*22 chvers
C--                         1234567890123456789012
      parameter ( chvers = 'HQSTF v01-01  15-05-12' )
      

