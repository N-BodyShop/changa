#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([ChaNGa], [2.0], [changa-users@u.washington.edu])
AC_CONFIG_SRCDIR([Reductions.h])
AC_CONFIG_HEADER([config.h
		  ../utility/structures/config.h])

AC_ARG_VAR([CHARMC], "Charm compiler")
# Checks for programs.
AC_PATH_PROG([CHARMC], [charmc], [no], [$PATH:$CHARM_DIR/bin:../charm/bin:../../charm/bin])

if test "$CHARMC" = "no" ; then
   AC_MSG_WARN([Cannnot locate charmc])
   AC_MSG_WARN([Try setting your CHARM_DIR environment variable to point at the charm source directory])
   exit
fi

CHARM_PATH=${CHARMC%/bin/charmc}
CONV_CONFIG=${CHARM_PATH}/tmp/conv-config.sh
CHARMINC=${CHARM_PATH}/include
AC_SUBST([CHARM_PATH])
. ${CONV_CONFIG}
CXX=${CMK_SEQ_CXX}
CC=${CMK_SEQ_CC}
LD=${CMK_SEQ_LD}
LDXX=${CMK_SEQ_LDXX}
AC_SUBST([AR], [${CMK_AR}])
CPPFLAGS="$CPPFLAGS -I.."
AC_CHECK_FILE(../libs, [LDFLAGS="$LDFLAGS -L../libs"])

AC_PROG_CXX($CHARMC)
# AC_PROG_CC($CHARMC)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([sys/bittypes.h endian.h limits.h malloc.h stddef.h stdint.h stdlib.h string.h sys/param.h unistd.h values.h])
AC_CHECK_HEADER(netinet/in.h,,[MISSING_NETINET=yes])
if test "${MISSING_NETINET}" = "yes" ; then
   CPPFLAGS="$CPPFLAGS -I/usr/include"
fi

AC_CHECK_HEADER(rpc/xdr.h,,[MISSING_XDR=yes])
AC_SEARCH_LIBS([xdr_float], [xdr],,[MISSING_XDR=yes])
if test "${MISSING_XDR}" = "yes" ; then
   echo "XDR missing, attempting to use xdr with package"
   if test -d ../xdr ; then
   	ln -s $PWD/../xdr ../rpc
   else
	echo "Try downloading the xdr package and place it in the parent directory."
	echo "And do a make in ../xdr"
	exit
   fi
   CPPFLAGS="$CPPFLAGS -I.."
   LDFLAGS="$LDFLAGS -L../xdr"
fi

# Cooling:
AC_ARG_ENABLE([cooling],
	[AS_HELP_STRING([--enable-cooling], [enable gas cooling (planet,cosmo)])],
	[cooling=$enableval], [cooling=no])
OBJECTS_COOLING=""
if test "$cooling" = "no" ; then
   echo "Adiabatic SPH only selected"
   FLAG_COOLING="-DCOOLING_NONE"
elif test "$cooling" = "planet" ; then
   echo "Planet cooling selected"
   FLAG_COOLING="-DCOOLING_PLANET"
   OBJECTS_COOLING="cooling_planet.o stiff.o"
   LIBS="$LIBS -lgsl -lgslcblas"
elif test "$cooling" = "cosmo" ; then
   echo "Cosmo cooling selected"
   FLAG_COOLING="-DCOOLING_COSMO"
   OBJECTS_COOLING="cooling_cosmo.o stiff.o"
   LIBS="$LIBS -lgsl -lgslcblas"
else
	echo "Unknown cooling selected, using Adiabatic"
	FLAG_COOLING="-DCOOLING_NONE"
	cooling="no"
fi

if test "$cooling" != "no" ; then
	AC_CHECK_HEADER(gsl/gsl_odeiv.h,,[MISSING_GSL=yes])
	AC_SEARCH_LIBS([gsl_odeiv_step_alloc], [gsl],,[MISSING_GSL=yes])
fi
if test "${MISSING_GSL}" = "yes" ; then
   echo "GSL missing, attempting to search in parent directory"
   for gsl in $PWD/../gsl $PWD/../gsl-1.14 $PWD/../gsl-1.13 ; do
	   if test -d $gsl ; then
		CPPFLAGS="$CPPFLAGS -I$gsl/include"
		CFLAGS="$CFLAGS -I$gsl/include"
		LDFLAGS="$LDFLAGS -L$gsl/lib"
		MISSING_GSL="no"
		LIBS="$LIBS -lgslcblas"
	   fi
   done
   if test "${MISSING_GSL}" = "yes" ; then
	echo GSL (GNU Scientific Libary) is missing for cooling code
	exit
   fi
fi

AC_SUBST([FLAG_COOLING])
AC_SUBST([OBJECTS_COOLING])
 
# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_TYPES([u_int64_t, u_quad_t, int64_t, quad_t, quad])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt strchr])
AC_CHECK_LIB([xdr], [xdr_hyper])
AC_CHECK_FUNCS([xdr_hyper])

AC_CONFIG_FILES([Makefile
		 ../utility/structures/Makefile
		])
AC_OUTPUT
