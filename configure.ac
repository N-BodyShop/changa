#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([ChaNGa], [3.3], [changa-users@u.washington.edu])
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

# Use single-precision for gravity calculations
AC_ARG_ENABLE([float],
	[AS_HELP_STRING([--enable-float], [use single-precision for gravity calculations])],
	[float=$enableval], [float=no])
if test "$float" = "yes" ; then
    echo "float selected"
    FLAG_FLOAT="-DCOSMO_FLOAT"
else
    FLAG_FLOAT=""
fi
AC_SUBST([FLAG_FLOAT])

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
	[dtadjust=$enableval], [dtadjust=yes])
if test "$dtadjust" = "yes" ; then
    echo "DTADJUST selected"
    FLAG_DTADJUST="-DDTADJUST"
else
    FLAG_DTADJUST=""
fi

AC_SUBST([FLAG_DTADJUST])

# Select kernel
# Set defaults
KERNEL="m4kernel"
FLAG_M4KERNEL="-DM4KERNEL=0"
FLAG_WENDLAND="-DWENDLAND=0"
FLAG_M6KERNEL="-DM6KERNEL=0"
# Wendland C4 Kernel:
AC_ARG_ENABLE([wendland],
	[AS_HELP_STRING([--enable-wendland], [enable Wendland C4 Kernel])],
	[wendland=$enableval], [wendland=no])
if test "$wendland" = "yes" ; then
    KERNEL="wendland"
fi
# M6 quintic spline kernel
AC_ARG_ENABLE([m6kernel],
	[AS_HELP_STRING([--enable-m6kernel], [enable M6 cubic spline kernel])],
	[m6kernel=$enableval], [m6kernel=no])
if test "$m6kernel" = "yes" ; then
    KERNEL="m6kernel"
fi
# Now set up kernel flags
if test "$KERNEL" = "m4kernel" ; then
    echo "M4KERNEL selected"
    FLAG_M4KERNEL="-DM4KERNEL=1"
elif test "$KERNEL" = "wendland" ; then
    echo "WENDLAND selected"
    FLAG_WENDLAND="-DWENDLAND=1"
elif test "$KERNEL" = "m6kernel" ; then
    echo "M6KERNEL selected"
    FLAG_M6KERNEL="-DM6KERNEL=1"
fi
KERNEL_FLAGS="$FLAG_M4KERNEL $FLAG_WENDLAND $FLAG_M6KERNEL"
echo "$KERNEL_FLAGS"
AC_SUBST([KERNEL_FLAGS])

# Velocity damping
AC_ARG_ENABLE([damping],
	[AS_HELP_STRING([--enable-damping], [enable velocity damping for glasses])],
	[damping=$enableval], [damping=no])
if test "$damping" = "yes" ; then
    echo "DAMPING selected"
    FLAG_DAMPING="-DDAMPING"
else
    FLAG_DAMPING=""
fi
AC_SUBST([FLAG_DAMPING])

# Cooling:
AC_ARG_ENABLE([cooling],
	[AS_HELP_STRING([--enable-cooling], [enable gas cooling (planet, cosmo, grackle)])],
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
        AC_MSG_WARN([Unknown cooling selected, using Adiabatic])
	FLAG_COOLING="-DCOOLING_NONE"
	cooling="no"
fi
AC_ARG_ENABLE([hexadecapole],
        [AS_HELP_STRING([--disable-hexadecapole], [disable use of hexadecapole expansions])],
        # Note that "$enableval" being "yes" means ENABLE hexadecapole -- the
        # help string, as one would expect, only affects how we document
        # the option.
        [case "$enableval" in
              yes | no ) hexadecapole="$enableval";;
              *) AC_MSG_ERROR("$enableval: invalid argument to `--disable-hexadecapole'") ;;
         esac],
        [hexadecapole="yes"])

HEXADECAPOLE=""
if test "$hexadecapole" = "yes"; then
    HEXADECAPOLE="-DHEXADECAPOLE"
else
    HEXADECAPOLE=""
fi

AC_SUBST([FLAG_COOLING])
AC_SUBST([OBJECTS_COOLING])
AC_SUBST([HEXADECAPOLE])

###################### CAMBRIDGE ###################
# CAMBRIDGE
AC_ARG_ENABLE([cambridge],
    [AS_HELP_STRING([--enable-cambridge], [Enable cambridge (yes, no)])],
    [cambridge=$enableval], [cambridge=no])
if test "$cambridge" = "yes" ; then
    echo "DCAMBRIDGE selected"
    FLAG_CAMBRIDGE="-DCAMBRIDGE"
else
    FLAG_CAMBRIDGE=""
fi

AC_SUBST([FLAG_CAMBRIDGE])
 
# diffusion (thermal and metal)
AC_ARG_ENABLE([diffusion],
	[AS_HELP_STRING([--enable-diffusion], [enable diffusion])],
	[diffusion=$enableval], [diffusion=yes])
if test "$diffusion" = "yes" ; then
    echo "DIFFUSION selected"
    FLAG_DIFFUSION="-DDIFFUSION"
else
    FLAG_DIFFUSION=""
fi

AC_SUBST([FLAG_DIFFUSION])

# Limit supernova feedback
# By default this also enables DIFFHARMONIC
AC_ARG_ENABLE([feedbacklimit],
    [AS_HELP_STRING([--enable-feedbacklimit], [limit diffusion of feedback energy])],
    [feedbacklimit=$enableval], [feedbacklimit=no])
if test "$feedbacklimit" = "yes" ; then
    echo "FEEDBACKDIFFLIMIT selected"
    FLAG_FEEDBACKDIFFLIMIT="-DFEEDBACKDIFFLIMIT"
    echo "DIFFUSIONHARMONIC selected"
    FLAG_DIFFHARMONIC="-DDIFFUSIONHARMONIC"
else
    FLAG_FEEDBACKDIFFLIMIT=""
    FLAG_DIFFHARMONIC=""
fi

AC_SUBST([FLAG_FEEDBACKDIFFLIMIT])
AC_SUBST([FLAG_DIFFHARMONIC])

# Calculate artificial viscosity from Cullen & Dehnen 2010
AC_ARG_ENABLE([cullenalpha],
        [AS_HELP_STRING([--enable-cullenalpha], [enable Cullen Dehnen artificial viscosity])],
        [cullenalpha=$enableval], [cullenalpha=no])
if test "$cullenalpha" = "yes" ; then
    echo "CULLENALPHA selected"
    FLAG_CULLENALPHA="-DCULLENALPHA"
else
    FLAG_CULLENALPHA=""
fi

AC_SUBST([FLAG_CULLENALPHA])

# Use alternative calculation of viscosity term similar to Morris & Monahan (1997)
AC_ARG_ENABLE([vsigvisc],
        [AS_HELP_STRING([--enable-vsigvisc], [enable alternative Monahan artificial viscosity])],
        [vsigvisc=$enableval], [vsigvisc=yes])
if test "$vsigvisc" = "yes" ; then
    echo "VSIGVISC selected"
    FLAG_VSIGVISC="-DVSIGVISC"
else
    FLAG_VSIGVISC=""
fi

AC_SUBST([FLAG_VSIGVISC])

# Richie-Thomas forces (Pressure SPH)
AC_ARG_ENABLE([rtforce],
	[AS_HELP_STRING([--enable-rtforce], [enable Richie-Thomas forces])],
	[rtforce=$enableval], [rtforce=yes])
if test "$rtforce" = "yes" ; then
    echo "RTFORCE selected"
    FLAG_RTFORCE="-DRTFORCE"
else
    FLAG_RTFORCE=""
fi

AC_SUBST([FLAG_RTFORCE])

# Compiler flags
AC_ARG_ENABLE([opts],
	[AS_HELP_STRING([--enable-opts="options"], [Set compiler flags (override defaults)])],
	[opts=$enableval], [])
if test "$opts" = "yes" ; then
    # The user called without arguments.  Assume this is an empty string
    opts=""
fi
# Test if 'opts' is defined
if test "${opts+x}" != ""  ; then
    OPTS="$opts"
else
    # Default options
    OPTS="-O3 -lpthread"
fi
echo "compiler options: $OPTS"
AC_SUBST([OPTS])

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

####################### CUDA #######################
m4_include([cuda.ac])
##################### END CUDA #####################


AC_CONFIG_FILES([Makefile
		 $srcdir/../utility/structures/Makefile
		])
AC_OUTPUT



