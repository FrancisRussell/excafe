#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(excafe, 0.1, fpr02@doc.ic.ac.uk)
AC_CONFIG_SRCDIR([navier_stokes/solver.cpp])
AC_CONFIG_HEADER([include/excafe/config.h])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_LIBTOOL
AC_LANG_CPLUSPLUS
CXXFLAGS="-O2 -Wall -pedantic -Wno-long-long"

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# C++ compiler for Excafe.
AC_DEFINE_UNQUOTED(EXCAFE_CXX_COMPILER, "${CXX}", [C++ compiler that Excafe will use at run-time.])

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Check how to compile MPI
ACX_MPI([], AC_MSG_ERROR([Unable to find MPI.]))
CC="$MPICC"
CXX="$MPICXX"

# Check for OpenMP
#AX_OPENMP(CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS")

# Check for Boost
AX_BOOST_BASE([1.33])
if test "$want_boost" = no; then
  AC_MSG_ERROR([Boost is a required dependency and cannot be disabled.])
fi
CXXFLAGS="$CXXFLAGS $BOOST_CPPFLAGS"

AX_BOOST_UNIT_TEST_FRAMEWORK

AC_CHECK_LIB(triangle, triangulate, [], AC_MSG_ERROR([Unable to find libtriangle.]))

AX_WITH_LIBRARY(gmp, gmp.h, gmp)
CXXFLAGS="$CXXFLAGS $GMP_CFLAGS"

AX_WITH_LIBRARY(petsc, petsc.h, petsc)
CXXFLAGS="$CXXFLAGS $PETSC_CFLAGS"

AX_WITH_LIBRARY(papi, papi.h, papi, optional)
CXXFLAGS="$CXXFLAGS $PAPI_CFLAGS"
AM_CONDITIONAL([HAVE_PAPI], [test $HAVE_PAPI -ne 0])

PKG_CHECK_MODULES([CLN], [cln])
CXXFLAGS="$CXXFLAGS $CLN_CFLAGS"
AC_SUBST(CLN_LIBS)

PKG_CHECK_MODULES([GINAC], [ginac])
CXXFLAGS="$CXXFLAGS $GINAC_CFLAGS"
AC_SUBST(GINAC_LIBS)

PKG_CHECK_MODULES([LIBXML2], [libxml-2.0])
CXXFLAGS="$CXXFLAGS $LIBXML2_CFLAGS"
AC_SUBST(LIBXML2_LIBS)

PKG_CHECK_MODULES([APR1], [apr-1])
CXXFLAGS="$CXXFLAGS $APR1_CFLAGS"
AC_SUBST(APR1_LIBS)

# Check for MTL
#AC_ARG_WITH([mtl],
#          AC_HELP_STRING([--with-mtl=PFX], [Use the MTL installation in PFX]),
#	  [
#	    if test -f "$withval/include/mtl/mtl.h"; then
#              CXXFLAGS="$CXXFLAGS -I$withval/include"
#	    else
#	      AC_MSG_ERROR("$withval does not appear to be a prefix for a MTL installation.")
#             fi
#          ]
#          )

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([pow sqrt])

#AC_MSG_CHECKING([for MTL])
#AC_TRY_COMPILE([@%:@include <mtl/mtl.h>], [], AC_MSG_RESULT([yes]), AC_MSG_RESULT([no]); AC_MSG_ERROR([Failed to find MTL.]))

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([navier_stokes/Makefile])
AC_CONFIG_FILES([heat/Makefile])
AC_CONFIG_FILES([mass_matrix_generator_2d/Makefile])
AC_CONFIG_FILES([steady_stokes_dolfin/Makefile])
AC_CONFIG_FILES([matmult_generator/Makefile])
AC_CONFIG_FILES([lib/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_OUTPUT
