#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
m4_define(dragon_version_major, 3)
m4_define(dragon_version_minor, 1)
m4_define(dragon_version_micro, 0)

AC_PREREQ(2.59)
AC_INIT([DRAGON],[dragon_version_major().dragon_version_minor().dragon_version_micro()],[luca.maccione@lmu.de])
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE
LT_INIT
AC_SUBST([LIBTOOL_DEPS])
AC_CONFIG_MACRO_DIR([m4])
ACLOCAL_AMFLAGS='-I m4'

DRAGON_VERSION_MAJOR=dragon_version_major()
DRAGON_VERSION_MINOR=dragon_version_minor()
DRAGON_VERSION_MICRO=dragon_version_micro()
AC_SUBST(DRAGON_VERSION_MAJOR)
AC_SUBST(DRAGON_VERSION_MINOR)
AC_SUBST(DRAGON_VERSION_MICRO)

AC_COPYRIGHT([
Copyright 2010.
All Rights Reserved.
This file is part of the open source package DRAGON which is distributed
under the Gnu LGPL.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(dragon.cc)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX([icpc g++ clang++])
AC_PROG_CC([icc gcc clang])
AC_PROG_F77([ifort gfortran g77 g95])
AC_PROG_INSTALL

# Use libraries in tools
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)
AC_LTDL_DLLIB

AC_ARG_ENABLE([optimization], [AS_HELP_STRING([--disable-optimization],[disable optimization [default=enabled] ])],
[AC_MSG_RESULT(no)]
,
[AC_MSG_RESULT(yes)]
[CXXFLAGS="${CXXFLAGS} -O3"]
)

AC_ARG_ENABLE([openmp], [AS_HELP_STRING([--disable-openmp],[disable parallelization with OpenMP [default=enabled] ])],
[AC_MSG_RESULT(no)]
[CXXFLAGS="${CXXFLAGS} -U_OPENMP"]
,
[AC_MSG_RESULT(yes)]
[
if test "$CXX" == "icpc"
then
CXXFLAGS="${CXXFLAGS} -openmp -mkl"
LDFLAGS_INI="${LDFLAGS}"
fi

if test "$CXX" == "g++"
then
CXXFLAGS="${CXXFLAGS} -fopenmp"
LDFLAGS_INI="${LDFLAGS} -lgomp"
fi
]
)

# Add mixing f77 with C feature
AC_F77_LIBRARY_LDFLAGS

# Checks for external libraries.
AC_CHECK_HEADER( [string], [CPPFLAGS="$CPPFLAGS -DTIXML_USE_STL" LDFLAGS_INI="${LDFLAGS_INI} -DTIXML_USE_STL"] )
#LDFLAGS_INI="$LDFLAGS -DTIXML_USE_STL"

#CPPFLAGS="$CPPFLAGS -DTIXML_USE_STL"

# -- GSL Test --
AC_ARG_WITH(
    [gsl-path],
    AS_HELP_STRING(
        [--with-gsl-path=DIR],
        [gsl-config binary is in DIR]
    )
)

AC_PATH_PROG(
    [GSL_CONFIG],
    [gsl-config],
    [notfound],
    [$with_gsl_path:$PATH]
)	
if test "$GSL_CONFIG" == "notfound"
  then
    AC_MSG_ERROR([gsl-config not found. See ./configure --help to know how to set a proper gsl path])
  fi
CPPFLAGS="$CPPFLAGS `$GSL_CONFIG --cflags`"
LDFLAGS_GSL="$LDFLAGS_INI `$GSL_CONFIG --libs`"
#  fi
LDFLAGS="$LDFLAGS_GSL"

AC_CHECK_HEADERS(
    [gsl/gsl_vector.h gsl/gsl_linalg.h gsl/gsl_errno.h gsl/gsl_spline.h],
    [AC_MSG_RESULT([check for gsl include files... yes])],
    [AC_MSG_ERROR([Cannot find some needed GSL include files])]
)

#AC_LINK_IFELSE(
#[AC_LANG_SOURCE([[
##include "gsl/gsl_vector.h"
#
#int main () {
#gsl_vector* x = gsl_vector_calloc(2) ;
#gsl_vector_free(x);
#return 0;
#}
#]]
#)],
#[AC_MSG_RESULT([check for gsl link... yes])],
#[AC_MSG_ERROR([Cannot link with GSL lib. Check GSL installation.])]
#)
# -- ROOT Test --
AC_ARG_WITH(
[root],
AS_HELP_STRING(
[--with-root(=DIR)],
[using Root output option, assuming Root installed in DIR (default is NO)]
)
)

AM_CONDITIONAL(ROOTCOMPIL,test "x$with_root" != "x")

if test "x$with_root" != "x" 
  then
  # Using root
    AC_PATH_PROG(
        [ROOT_CONFIG],
        [root-config],
        [notfound],
        [$with_root/bin:$ROOTSYS/bin:$PATH]
    )
    if test "$ROOT_CONFIG" == "notfound"
      then
        AC_MSG_ERROR([root-config not found. See ./configure --help to know how to set a proper root path])
    else
        CPPFLAGS="$CPPFLAGS `$ROOT_CONFIG --cflags` -DHAVE_ROOT"
        LDFLAGS_ROOT=" `$ROOT_CONFIG --libs` -DHAVE_ROOT "
        AC_DEFINE([HAVE_ROOT],[1], [Define if ROOT is present])

    fi
fi

AC_ARG_WITH(
[numcpu],
AS_HELP_STRING(
[--with-numcpu],
[Setting max number of CPU to be used]
)
)
if test "x$with_numcpu" != "x"
then
   AC_DEFINE_UNQUOTED([OMP_NUM_THREADS],[$with_numcpu], [Define the maximum number of CPU to be used])
fi

# -- CFITSIO Test --
LDFLAGS="$LDFLAGS_INI"

AC_ARG_WITH(
    [cfitsio],
    AS_HELP_STRING(
        [--with-cfitsio=DIR],
        [cfitsio library files are in DIR/include and DIR/lib]
    )
)
if test "x$with_cfitsio" != "x" 
  then
    CPPFLAGS="$CPPFLAGS -I$with_cfitsio/include"
    LDFLAGS_CFITSIO="$LDFLAGS -L$with_cfitsio/lib -lcfitsio"
else
    AC_MSG_ERROR([CFITSIO path not provided. See ./configure --help to see how to set a path properly.])
fi

AC_CHECK_HEADER(
[fitsio.h],
[],
[AC_MSG_ERROR([Cannot find fitsio.h. See ./configure --help to know how to set a proper cfitsio include directory])]
)

# Test CFITSIO_VERSION >= 3.006
#AC_RUN_IFELSE(
#AC_LANG_SOURCE([[
##include <iostream>
##include "fitsio.h"
#int main() {
#int r;
#r = !( CFITSIO_VERSION >= 3. ) ;
#if (r) std::cerr << "cfitsio version = " << CFITSIO_VERSION << std::endl;
#return r ;
#}
#]]),
#[AC_MSG_RESULT([checking for cfistio version... yes])],
#[AC_MSG_ERROR([Bad cfitsio version. Need >= 3.])]
#)


LDFLAGS="$LDFLAGS_CFITSIO"
AC_CHECK_LIB(
[cfitsio],
[ffvers],
[],
[AC_MSG_ERROR(
[Cannot find cfitsio library. See ./configure --help to know how to set a proper path]
)]
)

LDFLAGS="$LDFLAGS_INI $LDFLAGS_GSL $LDFLAGS_CFITSIO $LDFLAGS_ROOT $FLIBS"

LDFLAGS="$LDFLAGS_INI"


AC_ARG_WITH(
[numcpu],
AS_HELP_STRING(
[--with-numcpu],
[Setting max number of CPU to be used]
)
)
if test "x$with_numcpu" != "x"
then
   AC_DEFINE_UNQUOTED([NUMTHREADS],[$with_numcpu], [Define the maximum number of CPU to be used])
fi



# -- DarkSUSY Test --
AC_ARG_WITH(
    [ds],
    AS_HELP_STRING(
        [--with-ds=DIR],
        [DarkSUSY files are in DIR/include and DIR/lib.]
    )
)
if test "x$with_ds" != "x" 
  then
    FFLAGS="$FFLAGS -I$with_ds/include -DHAVE_DS"
    CPPFLAGS="${CPPFLAGS} -I$with_ds/include -DHAVE_DS"
    LDFLAGS_DS="$LDFLAGS -L$with_ds/lib -ldarksusy -lgfortran -lFH -lHB -DHAVE_DS"
    AC_DEFINE([HAVE_DS],[1], [Define if DS is present])
    AC_CHECK_HEADER(
        [dsge.h],
        [],
        [AC_MSG_ERROR([Cannot find dsge.h. See ./configure --help to know how to set a proper DarkSusy include directory])]
    )
else
    AC_MSG_WARN([DarkSUSY path not provided. See ./configure --help to see how to set a path properly.])
fi


LDFLAGS="$LDFLAGS_DS"

LDFLAGS="$LDFLAGS_INI $LDFLAGS_GSL $LDFLAGS_CFITSIO $LDFLAGS_ROOT $LDFLAGS_DS $FLIBS"


# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([iostream fstream sstream string vector cmath])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
#	AC_FUNC_MALLOC
AC_FUNC_STAT
#AC_CHECK_FUNCS([floor pow sqrt strchr strtoul])

# Checks for programs to build docs
#AC_CHECK_PROG(LATEX, latex, latex)
#AC_CHECK_PROG(DVIPS, dvips, dvips)
#AC_CHECK_PROG(DVIPDF, dvipdf, dvipdf)
AC_CHECK_PROG(DOXYGEN, doxygen, doxygen)
AC_DEFINE([TIXML_USE_STL],[1], [Define if TinyXML can use stl containers])

#AM_CONDITIONAL(LATEXTEST,test x$LATEX = xlatex)
#AM_CONDITIONAL(DVIPSTEST,test x$DVIPS = xdvips)
#AM_CONDITIONAL(DVIPDFTEST,test x$DVIPDF = xdvipdf)
AM_CONDITIONAL(DOXYGENTEST,test x$DOXYGEN = xdoxygen)
# global CPP flags
AC_SUBST(CPPFLAGS)

# global CXX flags
AC_SUBST(CXXFLAGS)

# Make makefiles
#
AC_CONFIG_FILES([
Makefile
cparamlib/Makefile
])
AC_OUTPUT

#	AC_CONFIG_FILES([
#			Makefile
#			docs/Makefile
#			])
#	AC_OUTPUT

