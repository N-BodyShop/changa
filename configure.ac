#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([ChaNGa], [3.2], [changa-users@u.washington.edu])
AC_DEFINE_UNQUOTED([NBODY_PACKAGE_VERSION], [$PACKAGE_VERSION])
AC_CONFIG_SRCDIR([Reductions.h])
AC_CONFIG_HEADER([config.h ../utility/structures/config.h])

AC_ARG_VAR([CHARMC], "Charm compiler")
# Checks for programs.
AC_PATH_PROG([CHARMC], [charmc], [no], [$CHARM_DIR/bin:$PATH:../charm/bin:../../charm/bin])

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

AC_CHECK_HEADERS([rpc/types.h])
AC_CHECK_HEADER([rpc/xdr.h],,[MISSING_XDR=yes],
[#ifdef HAVE_RPC_TYPES_H
#include <rpc/types.h>
#endif
])
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

# SSE2 vector optimization for gravity
AC_ARG_ENABLE([sse2],
	[AS_HELP_STRING([--enable-sse2], [enable sse2 gravity vectorization])],
	[sse2=$enableval], [sse2=no])
if test "$sse2" = "yes" ; then
    echo "SSE2 selected"
    FLAG_SSE="-DCMK_USE_SSE2"
else
    FLAG_SSE=""
fi
AC_SUBST([FLAG_SSE])

# AVX vector optimization for gravity
AC_ARG_ENABLE([avx],
	[AS_HELP_STRING([--enable-avx], [enable avx gravity vectorization])],
	[avx=$enableval], [avx=no])
if test "$avx" = "yes" ; then
    echo "AVX selected"
    FLAG_AVX="-DCMK_USE_AVX"
else
    FLAG_AVX=""
fi
AC_SUBST([FLAG_AVX])

# physical softening in comoving coordinates:
AC_ARG_ENABLE([changesoft],
	[AS_HELP_STRING([--enable-changesoft], [enable physical softening])],
	[changesoft=$enableval], [changesoft=yes])
if test "$changesoft" = "yes" ; then
    echo "CHANGESOFT selected"
    FLAG_CHANGESOFT="-DCHANGESOFT"
else
    FLAG_CHANGESOFT=""
fi
AC_SUBST([FLAG_CHANGESOFT])

# 128 bit keys:
AC_ARG_ENABLE([bigkeys],
	[AS_HELP_STRING([--enable-bigkeys], [enable 128 bit hash keys])],
	[bigkeys=$enableval], [bigkeys=no])
if test "$bigkeys" = "yes" ; then
    echo "BIGKEYS selected"
    FLAG_BIGKEYS="-DBIGKEYS"
else
    FLAG_BIGKEYS=""
fi
AC_SUBST([FLAG_BIGKEYS])

# DTADJUST
AC_ARG_ENABLE([dtadjust],
	[AS_HELP_STRING([--enable-dtadjust], [enable emergency timestep adjust])],
	[dtadjust=$enableval], [dtadjust=no])
if test "$dtadjust" = "yes" ; then
    echo "DTADJUST selected"
    FLAG_DTADJUST="-DDTADJUST"
else
    FLAG_DTADJUST=""
fi

AC_SUBST([FLAG_DTADJUST])

# Wendland C4 Kernel:
AC_ARG_ENABLE([wendland],
	[AS_HELP_STRING([--enable-wendland], [enable Wendland C4 Kernel])],
	[wendland=$enableval], [wendland=no])
if test "$wendland" = "yes" ; then
    echo "WENDLAND selected"
    FLAG_KERNEL="-DWENDLAND"
else
    FLAG_KERNEL=""
fi
AC_SUBST([FLAG_KERNEL])

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
elif test "$cooling" = "cosmo" ; then
   echo "Cosmo cooling selected"
   FLAG_COOLING="-DCOOLING_COSMO"
   OBJECTS_COOLING="cooling_cosmo.o stiff.o"
elif test "$cooling" = "grackle" ; then
   echo "Grackle cooling selected"
   FLAG_COOLING="-DCOOLING_GRACKLE"
   OBJECTS_COOLING="cooling_grackle.o"
   ac_cv_libhdf5="yes"
   AC_CHECK_LIB([hdf5], [H5Fcreate], [], [ac_cv_libhdf5=no]) 
   if test $ac_cv_libhdf5 = "no" ; then
       AC_MSG_WARN([HDF5 library needed for Grackle, but not found]);
       AC_MSG_WARN([Try setting LDFLAGS environment with -L<hd5libdir>]);
       exit
   fi
   ac_cv_libgrackle="yes"
   AC_CHECK_HEADER(grackle.h, [], [],
[[
#define CONFIG_BFLOAT_8
]])
   if test "$ac_cv_header_grackle_h" == no
   then
       AC_MSG_WARN([grackle header needed for Grackle, but not found])
       AC_MSG_WARN([Try setting CPPFLAGS environment with -I <grackledir>])
       exit
   fi
   AC_CHECK_LIB([grackle], [initialize_chemistry_data], [], [ac_cv_libgrackle=no], [-lhdf5]) 
   if test $ac_cv_libgrackle = "no" ; then
       AC_MSG_WARN([grackle library needed for Grackle, but not found]);
       AC_MSG_WARN([Try setting LDFLAGS environment with -L<gracklelibdir>]);
       AC_MSG_WARN([Also remember to set LD_LIBRARY_PATH when running, if grackle and hdf5 are dynamic libraries]);
       exit
   fi
else
	echo "Unknown cooling selected, using Adiabatic"
	FLAG_COOLING="-DCOOLING_NONE"
	cooling="no"
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
