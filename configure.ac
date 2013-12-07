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

dnl === WITHLIB_OPTION([opt_pfx], [outvar_pfx])
dnl ===   Defines --with-<opt_pfx>{include,lib}dir options which set
dnl ===   the variables <outvar_pfx>_{INCLUDES,LIBS}.
AC_DEFUN([WITHLIB_OPTION],
  [AC_ARG_WITH([$1includedir],
               AS_HELP_STRING([--with-$1includedir=DIR],
                              [use $2 includes from DIR]),
               $2_INCLUDES="-I$withval")
   AC_ARG_WITH([$1libdir],
               AS_HELP_STRING([--with-$1libdir=DIR],
                              [use $2 libraries from DIR]),
               [$2_LIBS="-L$withval"])])


dnl === check for PNG support ===

PNG_INCLUDES=""
PNG_LIBS=""
AC_PATH_PROGS(LIBPNG_CONFIG,
              [libpng-config libpng15-config libpng14-config libpng12-config])
if test -n "$LIBPNG_CONFIG"; then
  PNG_INCLUDES=`$LIBPNG_CONFIG --cflags`
  PNG_PREFIX=`$LIBPNG_CONFIG --prefix`
  if test "${PNG_PREFIX}/lib" != "/usr/lib" ; then
    PNG_LIBS="-L${PNG_PREFIX}/lib"
  fi
  PNG_LIBS="$PNG_LIBS `$LIBPNG_CONFIG --libs`"
fi

WITHLIB_OPTION([png], [PNG])

SAVED_CPPFLAGS=$CPPFLAGS
SAVED_LIBS=$LIBS
CPPFLAGS="$PNG_INCLUDES $CPPFLAGS"
LIBS="$PNG_LIBS $LIBS"

AC_CHECK_HEADER(png.h,
  AC_SEARCH_LIBS(png_get_libpng_ver, [png],
                 [test "$ac_cv_search_png_get_libpng_ver" = "none required" \
                    || PNG_LIBS="$PNG_LIBS $ac_cv_search_png_get_libpng_ver"
                  PNG_INCLUDES="$PNG_INCLUDES -DWEBP_HAVE_PNG"
                  AC_DEFINE(WEBP_HAVE_PNG, [1],
                            [Set to 1 if PNG library is installed])
                  png_support=yes
                 ],
                 [AC_MSG_WARN(Optional png library not found)
                  PNG_LIBS=""
                  PNG_INCLUDES=""
                 ],
                 [$MATH_LIBS]),
  [AC_MSG_WARN(png library not available - no png.h)
   PNG_LIBS=""
   PNG_INCLUDES=""
  ],
)
AC_SUBST(PNG_LIBS)
AC_SUBST(PNG_INCLUDES)

CPPFLAGS=$SAVED_CPPFLAGS
LIBS=$SAVED_LIBS

AC_OUTPUT(Makefile  src/Makefile src/stereo/Makefile src/utils/Makefile src/transforms/Makefile src/filters/Makefile src/hdr/Makefile src/segmentation/Makefile examples/Makefile src/motion/Makefile)
