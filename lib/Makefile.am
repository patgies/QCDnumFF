lib_LTLIBRARIES = libQCDNUM.la

libQCDNUM_la_SOURCES = 

libQCDNUM_la_LDFLAGS = $(VERSIONFLAGS)

libQCDNUM_la_LIBADD = $(FCLIBS) \
	$(top_builddir)/mbutil/libmbutil.la \
	$(top_builddir)/wstore/libwstore.la \
    $(top_builddir)/splint/libsplint.la \
	$(top_builddir)/qcdnum/libqcdnum.la \
	$(top_builddir)/zmstf/libzmstf.la \
	$(top_builddir)/hqstf/libhqstf.la \
	$(top_builddir)/CxxWrap/libQCDNUMCCWrap.la

## Hack to get linking done by CXX
libQCDNUM_la_LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
