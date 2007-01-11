#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([ChaNGa], [0.1], [changa-users@u.washington.edu])
AC_CONFIG_SRCDIR([Reductions.h])
AC_CONFIG_HEADER([config.h
		  ../structures/config.h])

AC_ARG_VAR([CHARMC], "Charm compiler")
# Checks for programs.
AC_PATH_PROG([CHARMC], [charmc], [no], [$PATH:$CHARM_DIR/bin:../charm/bin:../../charm/bin])

if test "$CHARMC" = "no" ; then
   AC_MSG_WARN([Cannnot locate charmc])
   AC_MSG_WARN([Try setting your CHARM_DIR environment variable to point at the charm source directory])
   exit
fi

AC_PROG_CXX($CHARMC)
# AC_PROG_CC($CHARMC)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([sys/bittypes.h endian.h limits.h malloc.h stddef.h stdint.h stdlib.h string.h sys/param.h unistd.h values.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_TYPES([u_int64_t, u_quad_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt strchr xdr_hyper])

AC_CONFIG_FILES([Makefile
		 ../structures/Makefile
		])
AC_OUTPUT
