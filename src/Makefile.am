SUBDIRS = utils transforms filters hdr segmentation stereo
#segment hdr filters transforms libit 

AM_CPPFLAGS = -I$(top_srcdir)/src
lib_LTLIBRARIES = libit.la

#common_HEADERS =
#common_HEADERS += libit/types.h
#commondir = $(includedir)/libit

libit_la_SOURCES =
#libitinclude_HEADERS =
#libitinclude_HEADERS += webp/encode.h

#noinst_HEADERS =
#noinst_HEADERS += webp/format_constants.h

libit_la_LIBADD =
libit_la_LIBADD += utils/libitutils.la
libit_la_LIBADD += transforms/libittransforms.la
libit_la_LIBADD += filters/libitfilters.la
libit_la_LIBADD += hdr/libithdr.la
libit_la_LIBADD += segmentation/libitseg.la
libit_la_LIBADD += stereo/libitstereo.la

