#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([ChaNGa], [3.3], [changa-users@u.washington.edu])
AC_DEFINE_UNQUOTED([NBODY_PACKAGE_VERSION], [$PACKAGE_VERSION])
AC_CONFIG_SRCDIR([Reductions.h])

# Save the user-provided flags
# NOTE: This must be done _before_ the Charm scripts are imported
#       as they provide their own values for these flags.
AC_SUBST([USER_CPPFLAGS], $CPPFLAGS)
AC_SUBST([USER_CFLAGS], $CFLAGS)
AC_SUBST([USER_CXXFLAGS], $CXXFLAGS)
AC_SUBST([USER_LDFLAGS], $LDFLAGS)

if test "x$STRUCT_DIR" = "x"; then
	STRUCT_DIR=$srcdir/../utility/structures
fi
AC_SUBST([STRUCT_DIR])

AC_CONFIG_HEADER([config.h $STRUCT_DIR/config.h])

AC_ARG_VAR([CHARMC], "Charm compiler")
# Checks for programs.
AC_PATH_PROG([CHARMC], [charmc], [no], [$CHARM_DIR/bin:$PATH:../charm/bin:../../charm/bin])

if test "$CHARMC" = "no" ; then
   AC_MSG_WARN([Cannnot locate charmc])
   AC_MSG_WARN([Try setting your CHARM_DIR environment variable to point at the charm source directory])
   exit
fi

# Make sure git is available
AC_CHECK_PROG(HAVE_GIT,[git],[yes],[no])
test "$HAVE_GIT" == "no" && AC_MSG_ERROR([git not found])

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
AC_SUBST([CMK_SMP], [${CMK_SMP}])

AC_PROG_CXX($CHARMC)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([sys/bittypes.h endian.h limits.h malloc.h stddef.h stdint.h stdlib.h string.h sys/param.h unistd.h values.h])
AC_CHECK_HEADER(netinet/in.h,,[MISSING_NETINET=yes])
if test "${MISSING_NETINET}" = "yes" ; then
   CPPFLAGS="$CPPFLAGS -I/usr/include"
fi

AC_CHECK_HEADERS([rpc/types.h])
AC_CHECK_HEADER([tirpc/netconfig.h],[HAVE_TIRPC=yes])
AC_CHECK_HEADER([rpc/xdr.h],,[MISSING_XDR=yes],
[#ifdef HAVE_RPC_TYPES_H
#include <rpc/types.h>
#endif
])
if test "${MISSING_XDR}" = "yes" -a "${HAVE_TIRPC}" = yes; then
   echo "XDR missing but TIRPC installed"
   CPPFLAGS="$CPPFLAGS -I/usr/include/tirpc"
   AC_SEARCH_LIBS([xdr_float], [xdr tirpc], [MISSING_XDR=no], [MISSING_XDR=yes])
else
   AC_SEARCH_LIBS([xdr_float], [xdr],,[MISSING_XDR=yes])
fi
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

# Set up the cha_commitid file
VOLD=$(test -e VERSION && cat VERSION)
VNEW=$(cd $srcdir && git describe --long --always)
if test -n "$VNEW" -a "$VOLD" != "$VNEW"; then
    echo "const char * const Cha_CommitID = \"$VNEW\";" > cha_commitid.c
fi

# Compiler-specific options
m4_include([compiler.ac])

# Require C++11
COMPILER_FLAG_REQUIRE([-std=c++11],[ChaNGa requires full C++11 support])

# Require dependency-generation flags
COMPILER_FLAG_REQUIRE([-MMD -MP -MF .dummy],[ChaNGa requires the dependency-generation flags '-MMD -MP -MF'])
rm -f .dummy

# Check if optimized debug builds are supported
COMPILER_FLAG_CHECK([-Og], [OPT_DEBUG], [-Og], [-O0])

# Set compiler type
COMPILER_MACRO_CHECK([__GNUC__],[IS_GNU]) # also true for Clang
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

# Select target CPU architecture (this determines SIMD type)
# --- Deprecated flags ---
ARG_ENABLE([sse2], [DEPRECATED. Use --enable-arch=sse2], [flag_sse_deprecated], [sse], [no])
ARG_ENABLE([avx], [DEPRECATED. Use --enable-arch=avx], [flag_avx_deprecated], [avx], [no])
# -----------------------
AC_ARG_ENABLE([arch], [AS_HELP_STRING([--enable-arch], [set compiler target architecture])],
			  [arch=$enableval], [arch=none])
# --enable-arch overrules when one of the deprecated flags is also given
if test x$arch = xnone -a x$flag_sse_deprecated != x; then arch=$flag_sse_deprecated; fi
if test x$arch = xnone -a x$flag_avx_deprecated != x; then arch=$flag_avx_deprecated; fi
case $arch in
	no|none ) FLAG_ARCH="" ;;
 	sse2    ) FLAG_ARCH=-DCMK_USE_SSE2 ;;
	avx     ) FLAG_ARCH=-DCMK_USE_AVX ;;
	*       ) AC_MSG_ERROR("invalid argument for '--enable-arch': $arch");;
esac
AC_SUBST([FLAG_ARCH])

# -----------------------------------------------------------------------------
# ---- Gravity Options --------------------------------------------------------
# -----------------------------------------------------------------------------

# Use single-precision for gravity calculations
ARG_ENABLE([float], [use single-precision gravity calculations], [FLAG_FLOAT], [-DCOSMO_FLOAT], [no])

# Hexadecapole expansion for gravity
ARG_ENABLE([hexadecapole], [hexadecapole expansions in gravity], [HEXADECAPOLE], [-DHEXADECAPOLE], [yes])

# physical softening in comoving coordinates
ARG_ENABLE([changesoft], [physical softening], [FLAG_CHANGESOFT], [-DCHANGESOFT], [yes])

# 128 bit keys:
ARG_ENABLE([bigkeys], [128 bit hash keys], [FLAG_BIGKEYS], [-DBIGKEYS], [no])

# DTADJUST
ARG_ENABLE([dtadjust], [emergency timestep adjust], [FLAG_DTADJUST], [-DDTADJUST], [yes])

# Interaction list version
AC_ARG_ENABLE([interlist],
	[AS_HELP_STRING([--enable-interlist], [interaction list version])],
	[interlist=$enableval], [interlist=2])
case $interlist in
	1|2) FLAG_INTERLIST="-DINTERLIST_VER=$interlist";;
	*  ) AC_MSG_ERROR([invalid argument for '--enable-interlist': $interlist]);;
esac
AC_SUBST([FLAG_INTERLIST])

# Tree building
AC_ARG_ENABLE([tree-build],
	[AS_HELP_STRING([--enable-tree-build], [tree building (merge-remote,split-phase)])],
	[tree_build=$enableval], [tree_build=merge-remote])
case $tree_build in
	# merges remote requests before sending; local trees built after requests sent
	merge-remote ) FLAG_TREE_BUILD="-DMERGE_REMOTE_REQUESTS";;
	# no merging of remote requests on PE; local trees built after requests sent
	split-phase  ) FLAG_TREE_BUILD="-DSPLIT_PHASE_TREE_BUILD";;
	*            ) AC_MSG_ERROR([invalid argument for '--enable-tree-build': $tree_build]);;
esac
AC_SUBST([FLAG_TREE_BUILD])

# -----------------------------------------------------------------------------
# ---- Gas Options ------------------------------------------------------------
# -----------------------------------------------------------------------------

# NSMOOTHINNER smoothing adjustment
ARG_ENABLE([nsmoothinner], [nsmoothinner feedback], [FLAG_NSMOOTHINNER], [-DNSMOOTHINNER], [no])

# Gas particle splitting
ARG_ENABLE([splitgas], [gas particle splitting], [FLAG_SPLITGAS], [-DSPLITGAS], [no])

# Select SPH interpolation kernel
# --- Deprecated flags ---
ARG_ENABLE([wendland], [DEPRECATED. Use --enable-sph-kernel=wendland], [flag_sph_wendland_deprecated], [wendland], [none])
ARG_ENABLE([m6kernel], [DEPRECATED. Use --enable-sph-kernel=m6], [flag_sph_m6kernel_deprecated], [m6kernel], [none])
#------------------------
AC_ARG_ENABLE([sph_kernel], [AS_HELP_STRING([--enable-sph-kernel], [SPH interpolation kernel (m4,m6,wendland)])],
				[sph_kernel=$enableval], [sph_kernel=sph_none])
if test x$sph_kernel = xsph_none -a x$flag_sph_m6kernel_deprecated != x; then sph_kernel=$flag_sph_m6kernel_deprecated; fi
if test x$sph_kernel = xsph_none -a x$flag_sph_wendland_deprecated != x; then sph_kernel=$flag_sph_wendland_deprecated; fi
case $sph_kernel in
	sph_none|"m4") FLAG_SPH_KERNEL="-DM4KERNEL=1 -DWENDLAND=0 -DM6KERNEL=0";;
	m6           ) FLAG_SPH_KERNEL="-DM4KERNEL=0 -DWENDLAND=0 -DM6KERNEL=1";;
	wendland     ) FLAG_SPH_KERNEL="-DM4KERNEL=0 -DWENDLAND=1 -DM6KERNEL=0";;
	*            ) AC_MSG_ERROR([invalid argument for '--enable-sph-kernel': $sph_kernel]);;
esac
AC_SUBST([FLAG_SPH_KERNEL])

# Cooling
AC_ARG_ENABLE([cooling],
	[AS_HELP_STRING([--enable-cooling], [gas cooling (none,planet,cosmo,boley,metal,H2,grackle)])],
	[cooling=$enableval], [cooling=none])
case $cooling in
	no|none	) FLAG_COOLING="-DCOOLING_NONE";;
	planet	) FLAG_COOLING="-DCOOLING_PLANET";;
	cosmo	) FLAG_COOLING="-DCOOLING_COSMO";;
	boley   ) FLAG_COOLING="-DCOOLING_BOLEY";;
	metal   ) FLAG_COOLING="-DCOOLING_METAL";;
	H2      ) FLAG_COOLING="-DCOOLING_MOLECULARH";;
	grackle	) FLAG_COOLING="-DCOOLING_GRACKLE";;
	*		) AC_MSG_ERROR([invalid argument for '--enable-cooling': $cooling]);;
esac
AC_SUBST([FLAG_COOLING])

if test "$cooling" = "grackle" ; then
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

# Velocity damping
ARG_ENABLE([damping], [velocity damping in glasses], [FLAG_DAMPING], [-DDAMPING], [no])

# diffusion (thermal and metal)
ARG_ENABLE([diffusion], [thermal and metal diffusion], [FLAG_DIFFUSION], [-DDIFFUSION], [yes])

# Limit supernova feedback
# By default this also enables DIFFHARMONIC
ARG_ENABLE([feedbacklimit], [limit diffusion of feedback energy], [FLAG_FEEDBACKDIFFLIMIT], [-DFEEDBACKDIFFLIMIT], [no])
if test x$FLAG_FEEDBACKDIFFLIMIT = x"-DFEEDBACKDIFFLIMIT"; then
	FLAG_DIFFHARMONIC="-DDIFFUSIONHARMONIC";
fi
AC_SUBST([FLAG_DIFFHARMONIC])

# SUPERBUBBLE FEEDBACK
ARG_ENABLE([superbubble], [superbubble feedback], [FLAG_SUPERBUBBLE], [-DSUPERBUBBLE], [no])
if test "$superbubble" = "yes" -a "$FLAG_DIFFUSION" = "no" ; then
    AC_MSG_ERROR([SUPERBUBBLE requires diffusion.  Enable it with --enable-diffusion]);
fi
if test "$superbubble" = "yes" -a "$FLAG_COOLING" = "-DCOOLING_NONE" ; then
    AC_MSG_ERROR([SUPERBUBBLE requires cooling.  Enable it with --enable-cooling]);
fi

# Calculate artificial viscosity from Cullen & Dehnen 2010
ARG_ENABLE([cullenalpha], [Cullen Dehnen artificial viscosity], [FLAG_CULLENALPHA], [-DCULLENALPHA], [no])

# Use alternative calculation of viscosity term similar to Morris & Monahan (1997)
ARG_ENABLE([vsigvisc], [alternative Monahan artificial viscosity], [FLAG_VSIGVISC], [-DVSIGVISC], [yes])

# Richie-Thomas forces (Pressure SPH)
ARG_ENABLE([rtforce], [Richie-Thomas forces], [FLAG_RTFORCE], [-DRTFORCE], [yes])

# SIDM interactions
ARG_ENABLE([sidminter], [SIDM interactions], [FLAG_SIDMINTERACT], [-DSIDMINTERACT], [no])

# -----------------------------------------------------------------------------
# ---- Charm++ Options --------------------------------------------------------
# -----------------------------------------------------------------------------

# Default load balancer
ARG_ENABLE([default-lb], [default load balancer], [DEFAULT_LB], [MultistepLB_notopo], [yes])

# Enable projections
ARG_ENABLE([projections], [Charm++ projections], [PROJECTIONS], [yes], [no])
if test x$PROJECTIONS = xyes -a x$CMK_TRACE_ENABLED = x0; then
	AC_MSG_ERROR([You must enable tracing in Charm++ before using it in ChaNGa])
fi

# -----------------------------------------------------------------------------
# ---- Compiler Options -------------------------------------------------------
# -----------------------------------------------------------------------------

AC_ARG_ENABLE([opts], [AS_HELP_STRING([--enable-opts], [DEPRECATED - Do not use])], [], [])

# Sanitizers
AC_ARG_ENABLE([sanitizer],
	[AS_HELP_STRING([--enable-sanitizer], [LLVM sanitizers])],
	[sanitizer=$enableval], [sanitizer=no])
case $sanitizer in
	no|none) FLAG_SANITIZER=no;;
	yes    ) FLAG_SANITIZER=address;;  # default to the AddressSanitizer
	*      ) FLAG_SANITIZER=$sanitizer;;
esac
if test ! x$FLAG_SANITIZER = xno; then
	COMPILER_FLAG_REQUIRE([-fsanitize=$FLAG_SANITIZER],[$CXX does not support the $FLAG_SANITIZER sanitizer])
	COMPILER_FLAG_REQUIRE([-fno-omit-frame-pointer],[-fno-omit-frame-pointer is required when using the sanitizers])
fi
AC_SUBST([FLAG_SANITIZER])

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

AC_CONFIG_FILES([Makefile $STRUCT_DIR/Makefile])
AC_OUTPUT

echo "***************************************************************"
echo "Configuration Summary for AC_PACKAGE_NAME AC_PACKAGE_VERSION:"
echo
echo "    "Charm path "      " $CHARM_PATH
echo "    "Charm compiler"   " $CHARMC
echo "    "C++ compiler"     " $CXX
echo
echo "    "Gravity Flags "   " $HEXADECAPOLE $FLAG_FLOAT $FLAG_SSE$FLAG_AVX $FLAG_CHANGESOFT $FLAG_DTADJUST $FLAG_TREE_BUILD
echo "    "SPH flags "       " $FLAG_SPH_KERNEL $FLAG_DAMPING $FLAG_COOLING $FLAG_DIFFUSION $FLAG_FEEDBACKDIFFLIMIT $FLAG_DIFFHARMONIC $FLAG_CULLENALPHA $FLAG_VSIGVISC $FLAG_RTFORCE
echo "    "Misc Flags "      " projections=$PROJECTIONS $FLAG_BIGKEYS sanitizer=$FLAG_SANITIZER
echo "    "Load balancer "   " $DEFAULT_LB
echo "    "CPU architecture" " $FLAG_ARCH

if test x$ENABLE_CUDA = xyes; then
	echo
	echo "    "CUDA directory" " $CUDA_DIR
	echo "    "CUDA nvcc "     " $NVCC_PATH
	
	if test ! x$FLAG_SANITIZER = xno; then
		echo
		echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		echo !
		echo ! You are using a sanitizer with CUDA. In order to avoid spurious errors,
		echo ! you will need to export the following before running ChaNGa
		echo !
		echo !"     " ASAN_OPTIONS=protect_shadow_gap=0:replace_intrin=0:detect_leaks=0
		echo !
		echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	fi
fi

echo
