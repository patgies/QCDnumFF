      logical Lfill7, Lwtini7
      common /lflag7/ Lfill7(-1:mset0), Lwtini7
      common /qstor7/ stor7(nwf0), bufy7(mxx0), buft7(mqq0+10)
      common /steer7/ isetf7(-1:mset0) ,itypf7(-1:mset0),
     +                ifrst7(-1:mset0) ,ilast7(-1:mset0),
     +                ikeyf7(-1:mset0)
      common /ilist7/ idPij7(mpp0,mord0,mstp0),
     +                mxord7(mstp0),isetp7(mstp0)
      common /wlist7/ idWijk7(2,2,mord0,4,mstp0), idEijk7(2,2,mord0),
     +                idAijk7(maa0,maa0,mord0+1,mstp0)
      common /flist7/ idPdfi7(2,12)
      common /incmt7/ inciy7, inciz7, incid7, incpt7
      common /pstor7/ pars7(nwp0)
