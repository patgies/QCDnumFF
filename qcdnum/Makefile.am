noinst_LTLIBRARIES = libqcdnum.la

SRC =  src/srcAlfa.f src/srcConvolution.f src/srcEvolNN.f src/srcEvolve.f src/srcEvolveOld.f src/srcFastEngine.f src/srcFillPijAij.f src/srcGlobalid.f src/srcGrid.f src/srcInterpolation.f src/srcMatch.f src/srcParamStore.f src/srcPdfFill.f src/srcPdfInOut.f src/srcPdfNN.f src/srcPdfTables.f src/srcPijFunctions.f src/srcQcards.f src/srcQcdInit.f src/srcQcdSpline.f src/srcQcdUtil.f src/srcTboxWeights.f src/srcThresholds.f src/srcTransform.f src/srcWorkspace.f

PIJ =  pij/ome.f pij/pij_lo.f pij/pij_nlo.f pij/vogelsang.f pij/wgplg.f pij/xa2hgp.f pij/xpij2p.f pij/xpns2p.f

USR =  usr/cxxwrappers.f usr/obsolete.f usr/usrGlobalId.f usr/usrchecks.f usr/usrcvol.f usr/usrerr.f usr/usrevnn.f usr/usrevol.f usr/usrfast.f usr/usrgrd.f usr/usrini.f usr/usrparams.f usr/usrpdf.f usr/usrqcards.f usr/usrstore.f usr/usrwgt.f

EXTRA_DIST = inc

libqcdnum_la_SOURCES = $(SRC) $(PIJ) $(USR)
