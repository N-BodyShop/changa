#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([ChaNGa], [3.2], [changa-users@u.washington.edu])
AC_DEFINE_UNQUOTED([NBODY_PACKAGE_VERSION], [$PACKAGE_VERSION])
AC_CONFIG_SRCDIR([Reductions.h])
AC_CONFIG_HEADER([config.h $srcdir/../utility/structures/config.h])

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

# ----------------------------------------------------------------
# CUDA support
#
# Here we define and process a `--with-cuda[=PREFIX]` option.
#
# If the option is not given, or is given without an argument or with an
# argument of "yes", we search for CUDA in the following locations:
#
#   /usr/local/cuda			standard CUDA install directory
#
#   /usr/lib/nvidia-cuda-toolkit 	install directory for the system-
#					installed CUDA toolkit package on Debian
#                                       and Debian-based systems
#
# There are some arguments to AC_MSG_* with ugly formatting here; the formatting
# is intentional (it produces nice, readable output in the generated `configure`
# script).
ENABLE_CUDA="auto"
AC_ARG_VAR([CUDA_DIR], [CUDA toolkit install directory])
AC_ARG_WITH([cuda],[AS_HELP_STRING([--with-cuda@<:@=PATH@:>@],
        [Use CUDA (installed in PATH) @<:@unspecified=auto@:>@])],
        [CUDA_DIR=$withval],[CUDA_DIR=""])

# Interpret "--with-cuda" (i.e. no argument) as "enable CUDA support, please."
if test "$CUDA_DIR" = "yes"; then
    ENABLE_CUDA="yes"
    CUDA_DIR=""
fi

if test "$ENABLE_CUDA" != "no"; then
    if test -z "$CUDA_DIR"; then
        CUDA_DIR_AUTOMATIC=yes
        AC_MSG_CHECKING([for CUDA toolkit directory])
        CUDA_CANDIDATE_DIRS="/usr/local/cuda /usr/lib/nvidia-cuda-toolkit"

        for dir in $CUDA_CANDIDATE_DIRS; do
            if test -d "$dir"; then
                CUDA_DIR="$dir"
                break
            fi
        done
        if test -z "$CUDA_DIR"; then
            if test "$ENABLE_CUDA" = "yes"; then
                # User explicitly requested CUDA support; being unable to find a toolkit
                # directory is an error.
                AC_MSG_RESULT(not found)
                AC_MSG_ERROR(\
[CUDA support explicitly enabled but the toolkit was not found by automatic path
    search.  Specify CUDA toolkit location with `--with-cuda=PATH`.])

            else
                AC_MSG_RESULT([not found.  CUDA support will be disabled.])
                ENABLE_CUDA=no
            fi
        else
            AC_MSG_RESULT($CUDA_DIR)
        fi
    fi
fi

if test "$ENABLE_CUDA" != "no"; then
    AC_PATH_PROG([NVCC_PATH],nvcc,[],[$CUDA_DIR/bin])
    if test -z "$NVCC_PATH" ; then
        if test "$ENABLE_CUDA" = "yes"; then
            if test "$CUDA_DIR_AUTOMATIC" = "yes"; then
                AC_MSG_ERROR(\
[`bin/nvcc` not found in automatically-selected CUDA toolkit directory.  Specify
    CUDA toolkit location with `--with-cuda=PATH`.])

            else
                AC_MSG_ERROR([`bin/nvcc` not found in specified CUDA install prefix.])
            fi
        else                    # ENABLE_CUDA is "auto"
            AC_MSG_NOTICE(\
[`bin/nvcc` not found in automatically-selected CUDA toolkit directory; CUDA
    support will be disabled.  If this is not what you want, try specifying a
    valid CUDA toolkit install prefix using `--with-cuda=PATH`.])
            ENABLE_CUDA="no"
        fi
    fi
fi

# Set up `ENABLE_CUDA` for use by Make.
if test "$ENABLE_CUDA" = "no"; then
    ENABLE_CUDA=""
else
    ENABLE_CUDA="yes"
fi


AC_SUBST([CUDA_DIR])
AC_SUBST([NVCC_PATH])
AC_SUBST([ENABLE_CUDA])

AC_CONFIG_FILES([Makefile
		 $srcdir/../utility/structures/Makefile
		])
AC_OUTPUT
