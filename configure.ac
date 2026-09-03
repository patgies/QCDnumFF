AC_PREREQ([2.69])
AC_INIT([QCDNUM],[18-00-00],[h24@nikhef.nl],[qcdnum])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_MAINTAINER_MODE([disable])

AC_PROG_CXX
AC_CONFIG_HEADERS([CxxHead/QCDNUM/FortranWrappers.h])

## Set Fortran compiler behaviour
if test "x$FCFLAGS" == "x"; then
  FCFLAGS="-O3  -Iinc"
#  FCFLAGS="-c -Wall -O2 -Iinc"
fi
# Try to respect users' Fortran compiler variables
if test "x$FC" == "x"; then
  if test "x$F77" == "x"; then
    FC="$GFORTRAN"
  else
    FC="$F77"
  fi
fi
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported:\
 use a different --prefix argument])
fi

## Hack to make libtool use the FC Fortran compiler via F77, since
## libtool < 2.x doesn't understand --tag=FC
F77=$FC
AM_FFLAGS=$AM_FCFLAGS
FFLAGS="$FFLAGS $FCFLAGS"
AC_SUBST(FCFLAGS)
AC_SUBST(FFLAGS)
AC_SUBST(F77)
AC_SUBST(AM_FFLAGS)

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

## Set final FCFLAGS, CXXFLAGS and CPPFLAGS
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(top_srcdir)/CxxHead -I\$(top_builddir)/CxxHead"
AM_CPPFLAGS=["$AM_CPPFLAGS -DDATA_PATH="$datadir" -DQCDNUM_VERSION="$PACKAGE_VERSION" "]
AC_SUBST(AM_CPPFLAGS)

#AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall "])
AC_SUBST(AM_CXXFLAGS)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([CxxHead/Makefile CxxHead/QCDNUM/Makefile])
AC_CONFIG_FILES([CxxWrap/Makefile])
AC_CONFIG_FILES([mbutil/Makefile])
AC_CONFIG_FILES([wstore/Makefile])
AC_CONFIG_FILES([splint/Makefile])
AC_CONFIG_FILES([qcdnum/Makefile])
AC_CONFIG_FILES([zmstf/Makefile])
AC_CONFIG_FILES([hqstf/Makefile])
AC_CONFIG_FILES([lib/Makefile])
AC_CONFIG_FILES([bin/Makefile bin/qcdnum-config])

AC_OUTPUT
