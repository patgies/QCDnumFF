
C--   file zmcij.f containing the interfaces to cij routines
C--   the cij coefficient functions are stored in ../cij

C--   double precision function dzmAchi(qmu2)
C--   double precision function dzmConst1(x,qmu2,nf)
C--   double precision function dzmC2G(x,qmu2,nf)
C--   double precision function dzmC2Q(x,qmu2,nf)
C--   double precision function dzmCLG(x,qmu2,nf)
C--   double precision function dzmCLQ(x,qmu2,nf)
C--   double precision function dzmD3Q(x,qmu2,nf)
C--   double precision function dzmC2NN2A(x,qmu2,nf)
C--   double precision function dzmC2NS2B(x,qmu2,nf)
C--   double precision function dzmC2NN2C(x,qmu2,nf)
C--   double precision function dzmC2S2A(x,qmu2,nf)
C--   double precision function dzmC2G2A(x,qmu2,nf)
C--   double precision function dzmC2G2C(x,qmu2,nf)
C--   double precision function dzmCLNN2A(x,qmu2,nf)
C--   double precision function dzmCLNN2C(x,qmu2,nf)
C--   double precision function dzmCLNC2A(x,qmu2,nf)
C--   double precision function dzmCLNC2C(x,qmu2,nf)
C--   double precision function dzmCLS2A(x,qmu2,nf)
C--   double precision function dzmCLG2A(x,qmu2,nf)
C--   double precision function dzmC3NP2A(x,qmu2,nf)
C--   double precision function dzmC3NS2B(x,qmu2,nf)
C--   double precision function dzmC3NP2C(x,qmu2,nf)
C--   double precision function dzmC3NM2A(x,qmu2,nf)
C--   double precision function dzmC3NM2C(x,qmu2,nf)

C--   ----------------------------------------------------------------
C--   Rescaling variable chi = x    Jawohl
C--   ----------------------------------------------------------------

C     =======================================
      double precision function dzmAchi(qmu2)
C     =======================================

C--   Returns coefficient 'a' of scaling function chi = a*x
C--   For zero-mass coefficient functions chi = x and a = 1

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy   = qmu2     !avoid compiler warning
      dzmAchi = 1.D0

      return
      end
      
C--   ----------------------------------------------------------------
C--   Coefficient of LO delta(1-x) terms
C--   ----------------------------------------------------------------      

C     ==============================================
      double precision function dzmConst1(x,qmu2,nf)
C     ==============================================

C--   Returns the value of 1.D0

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

C--   Avoid compiler warnings
      dummy     = x
      dummy     = qmu2
      idumm     = nf
      dzmConst1 = 1.D0

      return
      end
      
C--   ----------------------------------------------------------------
C--   NLO F2
C--   ----------------------------------------------------------------      

C     ===========================================
      double precision function dzmC2G(x,qmu2,nf)
C     ===========================================

C--   Returns the value of ../cij/C2G.f   -- C2G(x,nf) regular

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy  = qmu2        !avoid compiler warning
      dzmC2G = CEEJ2G(x,nf)

      return
      end

C     ===========================================
      double precision function dzmC2Q(x,qmu2,nf)
C     ===========================================

C--   Returns the value of ../cij/C2Q.f  -- C2Q(x)     singular

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy  = qmu2        !avoid compiler warning
      dzmC2Q = CEEJ2Q(x,nf)

      return
      end
      
C--   ----------------------------------------------------------------
C--   NLO FL  or  LO FL'
C--   ----------------------------------------------------------------      

C     ===========================================
      double precision function dzmCLG(x,qmu2,nf)
C     ===========================================

C--   Returns the value of ../cij/CLG.f  -- CLG(x,nf)  regular

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy  = qmu2        !avoid compiler warning      
      dzmCLG = CEEJLG(x,nf)

      return
      end

C     ===========================================
      double precision function dzmCLQ(x,qmu2,nf)
C     ===========================================

C--   Returns the value of ../cij/CLQ.f  -- CLQ(x)     regular

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy  = qmu2        !avoid compiler warning
      dzmCLQ = CEEJLQ(x,nf)

      return
      end
      
C--   ----------------------------------------------------------------
C--   NLO xF3
C--   ----------------------------------------------------------------      

C     ===========================================
      double precision function dzmD3Q(x,qmu2,nf)
C     ===========================================

C--   Returns the value of minus ../cij/D3Q.f  -- D3Q(x)  regular 
C--   Minus sign is necessary because D3Q is subtracted from C2Q

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy  =  qmu2        !avoid compiler warning
      dzmD3Q = -DEEJ3Q(x,nf)

      return
      end
      
C--   ----------------------------------------------------------------
C--   NNLO F2
C--   ----------------------------------------------------------------      

C     ==============================================
      double precision function dzmC2NN2A(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C2NN2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC2NN2A = C2NN2A(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmC2NS2B(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C2NS2B(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC2NS2B = C2NS2B(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmC2NN2C(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C2NN2C(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC2NN2C = C2NN2C(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmC2NC2A(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C2NC2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC2NC2A = C2NC2A(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmC2NC2C(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C2NC2C(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC2NC2C = C2NC2C(x,nf)/4.D0

      return
      end

C     =============================================
      double precision function dzmC2S2A(x,qmu2,nf)
C     =============================================

C--   Returns the value of ../cij/C2S2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy    = qmu2        !avoid compiler warning
      dzmC2S2A = C2S2A(x,nf)/4.D0

      return
      end

C     =============================================
      double precision function dzmC2G2A(x,qmu2,nf)
C     =============================================

C--   Returns the value of ../cij/C2G2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy    = qmu2        !avoid compiler warning
      dzmC2G2A = C2G2A(x,nf)/4.D0

      return
      end

C     =============================================
      double precision function dzmC2G2C(x,qmu2,nf)
C     =============================================

C--   Returns the value of ../cij/C2G2C(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy    = qmu2        !avoid compiler warning
      dzmC2G2C = C2G2C(x,nf)/4.D0

      return
      end
      
C--   ----------------------------------------------------------------
C--   NNLO FL  or  NLO FL'
C--   ----------------------------------------------------------------      

C     ==============================================
      double precision function dzmCLNN2A(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/CLNN2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmCLNN2A = CLNN2A(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmCLNN2C(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/CLNN2C(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      idumm     = nf
      dzmCLNN2C = CLNN2C(x)/4.D0

      return
      end

C     ==============================================
      double precision function dzmCLNC2A(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/CLNC2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmCLNC2A = CLNC2A(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmCLNC2C(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/CLNC2C(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      idumm     = nf
      dzmCLNC2C = CLNC2C(x)/4.D0

      return
      end

C     =============================================
      double precision function dzmCLS2A(x,qmu2,nf)
C     =============================================

C--   Returns the value of ../cij/CLS2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy    = qmu2        !avoid compiler warning
      dzmCLS2A = CLS2A(x,nf)/4.D0

      return
      end

C     =============================================
      double precision function dzmCLG2A(x,qmu2,nf)
C     =============================================

C--   Returns the value of ../cij/CLG2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy    = qmu2        !avoid compiler warning
      dzmCLG2A = CLG2A(x,nf)/4.D0

      return
      end
      
C--   ----------------------------------------------------------------
C--   NNLO xF3
C--   ----------------------------------------------------------------      

C     ==============================================
      double precision function dzmC3NP2A(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C3NP2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC3NP2A = C3NP2A(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmC3NS2B(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C3NS2B(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC3NS2B = C3NS2B(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmC3NP2C(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C3NP2C(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC3NP2C = C3NP2C(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmC3NM2A(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C3NM2A(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC3NM2A = C3NM2A(x,nf)/4.D0

      return
      end

C     ==============================================
      double precision function dzmC3NM2C(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/C3NM2C(x,nf)
C--   Divide by 4 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmC3NM2C = C3NM2C(x,nf)/4.D0

      return
      end
      
C--   ----------------------------------------------------------------
C--   NNLO FL'  3-loop terms
C--   ----------------------------------------------------------------      

C     =============================================
      double precision function dzmCLG3A(x,qmu2,nf)
C     =============================================

C--   Returns the value of ../cij/CLG3A(x,nf)
C--   Divide by 8 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy    = qmu2        !avoid compiler warning
      dzmCLG3A = CLG3A(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dzmCLS3A(x,qmu2,nf)
C     =============================================

C--   Returns the value of ../cij/CLS3A(x,nf)
C--   Divide by 8 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy    = qmu2        !avoid compiler warning
      dzmCLS3A = CLS3A(x,nf)/8.D0

      return
      end
      
C     ==============================================
      double precision function dzmCLNP3A(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/CLNP3A(x,nf)
C--   Divide by 8 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmCLNP3A = CLNP3A(x,nf)/8.D0

      return
      end
      
C     ==============================================
      double precision function dzmCLNP3C(x,qmu2,nf)
C     ==============================================

C--   Returns the value of ../cij/CLNP3C(x,nf)
C--   Divide by 8 because of alfa/4pi instead of alfa/2pi

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)

      dummy     = qmu2        !avoid compiler warning
      dzmCLNP3C = CLNP3C(x,nf)/8.D0

      return
      end
      
