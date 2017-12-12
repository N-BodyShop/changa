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

MINIMUM_CHARM_VERSION=60800
if test $($CHARMC -V | awk '{print $3}') -lt $MINIMUM_CHARM_VERSION; then
	AC_MSG_ERROR([[At least Charm++ $MINIMUM_CHARM_VERSION is required]])
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

# Compiler-specific options
m4_include([compiler.ac])

# Require C++11
COMPILER_FLAG_REQUIRE([-std=c++11],[ChaNGa requires full C++11 support])

# Check if optimized debug builds are supported
COMPILER_FLAG_CHECK([-Og], [OPT_DEBUG], [-Og], [-O0])

# Set compiler type
COMPILER_MACRO_CHECK([__GNUC__],[IS_GNU])
COMPILER_MACRO_CHECK([__ICC__],[IS_ICC])

# @synopsis ARG_ENABLE(name, description, flag_var, val_if_enabled, default_val)
# @summary Wrapper around AC_ARG_ENABLE that does extra checking and automatic exporting
AC_DEFUN([ARG_ENABLE],
[dnl
	AC_ARG_ENABLE([$1], [AS_HELP_STRING([--enable-$1], [$2])],
	[case "$enableval" in
		  yes | no ) val=$enableval;;
		  *) AC_MSG_ERROR([invalid argument for '--enable-$1': $enableval]);;
	 esac],
	[val=$5])
	if test x$val = xyes; then $3=$4; else $3=""; fi
	AC_SUBST([$3])
])

# Use single-precision for gravity calculations
ARG_ENABLE([float], [use single-precision gravity calculations], [FLAG_FLOAT], [-DCOSMO_FLOAT], [no])

# Select target CPU architecture (this determines SIMD type)
# --- Deprecated flags ---
ARG_ENABLE([sse2], [Deprecated. Use --enable-arch=sse], [flag_sse_deprecated], [sse], [no])
ARG_ENABLE([avx], [Deprecated. Use --enable-arch=avx], [flag_avx_deprecated], [avx], [no])
# -----------------------
AC_ARG_ENABLE([arch], [AS_HELP_STRING([--enable-arch], [Set compiler target architecture])],
			  [arch=$enableval], [arch=none])
# --enable-arch overrules when one of the deprecated flags is also given
if test x$arch = xnone -a x$flag_sse_deprecated != x; then arch=$flag_sse_deprecated; fi
if test x$arch = xnone -a x$flag_avx_deprecated != x; then arch=$flag_avx_deprecated; fi
case $arch in
	no|none	) FLAG_ARCH="" ;;
 	yes|sse	) FLAG_ARCH=-DCMK_USE_SSE2 ;;
	avx		) FLAG_ARCH=-DCMK_USE_AVX ;;
	*		) AC_MSG_ERROR("invalid argument for '--enable-arch': $arch");;
esac
AC_SUBST([FLAG_ARCH])

# physical softening in comoving coordinates
ARG_ENABLE([changesoft], [enable physical softening], [FLAG_CHANGESOFT], [-DCHANGESOFT], [yes])

# 128 bit keys:
ARG_ENABLE([bigkeys], [enable 128 bit hash keys], [FLAG_BIGKEYS], [-DBIGKEYS], [no])

# DTADJUST
ARG_ENABLE([dtadjust], [enable emergency timestep adjust], [FLAG_DTADJUST], [-DDTADJUST], [yes])

# Select SPH interpolation kernel
# --- Deprecated flags ---
ARG_ENABLE([wendland], [Deprecated. Use --enable-sph-kernel=wendland], [flag_sph_wendland_deprecated], [wendland], [no])
ARG_ENABLE([m6kernel], [Deprecated. Use --enable-sph-kernel=m6kernel], [flag_sph_m6kernel_deprecated], [m6kernel], [no])
#------------------------
AC_ARG_ENABLE([sph-kernel], [AS_HELP_STRING([--enable-sph-kernel], [SPH interpolation kernel])],
				[sph_kernel=$enableval], [sph_kernel=none])
if test x$sph_kernel = xnone -a x$flag_sph_m6kernel_deprecated != x; then sph_kernel=$flag_sph_m6kernel_deprecated; fi
if test x$sph_kernel = xnone -a x$flag_sph_wendland_deprecated != x; then sph_kernel=$flag_sph_wendland_deprecated; fi
case $sph_kernel in
	yes|none|"m4kernel"	) FLAG_SPH_KERNEL="-DM4KERNEL=1 -DWENDLAND=0 -DM6KERNEL=0";;
	wendland	 		) FLAG_SPH_KERNEL="-DM4KERNEL=0 -DWENDLAND=1 -DM6KERNEL=0";;
	m6kernel	 		) FLAG_SPH_KERNEL="-DM4KERNEL=0 -DWENDLAND=0 -DM6KERNEL=1";;
	*			 		) AC_MSG_ERROR([invalid argument for '--enable-sph-kernel': $sph_kernel]);;
esac
AC_SUBST([FLAG_SPH_KERNEL])

# Velocity damping
ARG_ENABLE([damping], [enable velocity damping in gasses], [FLAG_DAMPING], [-DDAMPING], [no])

# Hexadecapole expansion for gravity
ARG_ENABLE([hexadecapole], [enable hexadecapole expansions], [HEXADECAPOLE], [-DHEXADECAPOLE], [yes])
 
# diffusion (thermal and metal)
ARG_ENABLE([diffusion], [enable thermal and metal diffusion], [FLAG_DIFFUSION], [-DDIFFUSION], [yes])

# Limit supernova feedback
# By default this also enables DIFFHARMONIC
ARG_ENABLE([feedbacklimit], [limit diffusion of feedback energy], [FLAG_FEEDBACKDIFFLIMIT], [-DFEEDBACKDIFFLIMIT], [no])
if test x$FLAG_FEEDBACKDIFFLIMIT = x"-DFEEDBACKDIFFLIMIT"; then
	FLAG_DIFFHARMONIC="-DDIFFUSIONHARMONIC";
fi
AC_SUBST([FLAG_DIFFHARMONIC])

# Calculate artificial viscosity from Cullen & Dehnen 2010
ARG_ENABLE([cullenalpha], [enable Cullen Dehnen artificial viscosity], [FLAG_CULLENALPHA], [-DCULLENALPHA], [no])

# Use alternative calculation of viscosity term similar to Morris & Monahan (1997)
ARG_ENABLE([vsigvisc], [enable alternative Monahan artificial viscosity], [FLAG_VSIGVISC], [-DVSIGVISC], [yes])

# Richie-Thomas forces (Pressure SPH)
ARG_ENABLE([rtforce], [enable Richie-Thomas forces], [FLAG_RTFORCE], [-DRTFORCE], [yes])

# Default load balancer
ARG_ENABLE([default_lb], [Default load balancer], [DEFAULT_LB], [MultistepLB_notopo], [yes])

# Enable projections
ARG_ENABLE([projections], [enable projections], [PROJECTIONS], [yes], [no])
if test x$PROJECTIONS = xyes -a x$CMK_TRACE_ENABLED = x0; then
	AC_MSG_ERROR([You must enable tracing in Charm++ before using it in ChaNGa])
fi

# Cooling
AC_ARG_ENABLE([cooling],
	[AS_HELP_STRING([--enable-cooling], [enable gas cooling (planet, cosmo, grackle)])],
	[cooling=$enableval], [cooling=no])
OBJECTS_COOLING=""
if test "$cooling" = "no" ; then
   FLAG_COOLING="-DCOOLING_NONE"
elif test "$cooling" = "planet" ; then
   FLAG_COOLING="-DCOOLING_PLANET"
   OBJECTS_COOLING="cooling_planet.o stiff.o"
elif test "$cooling" = "cosmo" ; then
   FLAG_COOLING="-DCOOLING_COSMO"
   OBJECTS_COOLING="cooling_cosmo.o stiff.o"
elif test "$cooling" = "grackle" ; then
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
fi

# GPU local tree walk
AC_ARG_ENABLE([gpu_local_tree_walk],
  [AS_HELP_STRING([--enable-gpu_local_tree_walk], [Enable gpu_local_tree_walk (yes, no)])],
  [gpu_local_tree_walk=$enableval], [gpu_local_tree_walk=no])
GPU_LOCAL_TREE_WALK=""
if test "$gpu_local_tree_walk" = "yes" ; then
  echo "DGPU_LOCAL_TREE_WALK selected"
  GPU_LOCAL_TREE_WALK="-DGPU_LOCAL_TREE_WALK"
else
  GPU_LOCAL_TREE_WALK=""
fi

AC_SUBST([GPU_LOCAL_TREE_WALK])

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

echo "***************************************************************"
echo "Configuration Summary for AC_PACKAGE_NAME AC_PACKAGE_VERSION:"
echo
echo "    "Charm path "      " $CHARM_PATH
echo "    "Charm compiler"   " $CHARMC
echo
echo "    "Gravity Flags "   " $HEXADECAPOLE $FLAG_FLOAT $FLAG_SSE$FLAG_AVX $FLAG_CHANGESOFT $FLAG_DTADJUST
echo "    "SPH flags "       " $FLAG_SPH_KERNEL $FLAG_DAMPING $FLAG_COOLING $FLAG_DIFFUSION $FLAG_FEEDBACKDIFFLIMIT $FLAG_DIFFHARMONIC $FLAG_CULLENALPHA $FLAG_VSIGVISC $FLAG_RTFORCE
echo "    "Misc Flags "      " $FLAG_BIGKEYS $PROJECTIONS
echo "    "Load balancer "   " $DEFAULT_LB
echo "    "CPU architecture" " $FLAG_ARCH

if test x$ENABLE_CUDA = xyes; then
echo
echo "    "CUDA directory" " $CUDA_DIR
echo "    "CUDA SDK "      " $CUDA_SDK_DIR
echo "    "CUDA nvcc "     " $NVCC_PATH
fi
