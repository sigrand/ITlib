AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CPPFLAGS += -I$(top_srcdir)/src/utils

bin_PROGRAMS = itlib

itlib_SOURCES = itlib.c
#cwebp_SOURCES += jpegdec.c jpegdec.h
#cwebp_SOURCES += pngdec.c pngdec.h

itlib_CPPFLAGS  = $(AM_CPPFLAGS) $(PNG_INCLUDES)
#cwebp_CPPFLAGS += $(JPEG_INCLUDES) $(PNG_INCLUDES)
itlib_LDADD = ../src/libit.la  $(PNG_LIBS)
#$(JPEG_LIBS) $(PNG_LIBS) $(TIFF_LIBS)
