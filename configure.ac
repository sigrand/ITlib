AC_INIT([ITlib], [0.01.0], [vadim@sigrand.com], [ITlib], [http://www.sigrand.com])

AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.10 no-define])
AC_CONFIG_HEADERS([config.h])
#AC_CONFIG_FILES([Makefile])

#AC_ISC_POSIX
AC_PROG_CC
AM_PROG_CC_C_O
#AC_HEADER_STDC
AM_PROG_LIBTOOL
AC_PROG_RANLIB
AC_PROG_LIBTOOL

AC_ARG_ENABLE(debug,[--enable-debug  Enable debug notices],[debug=$enableval], [debug=no])

#GCC_FLAGS="-Wall"

if test "x${debug}" = "xyes"; then
  CFLAGS="-Wall -O0 -g"
else
  CFLAGS="-O2"
fi
AC_SUBST(CFLAGS)

AC_OUTPUT(Makefile  src/Makefile src/utils/Makefile examples/Makefile )
