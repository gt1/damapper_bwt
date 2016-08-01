ACLOCAL_AMFLAGS=-I m4

noinst_HEADERS = map.h prs_uint64_t_pair.h prs.ic prs.ih

MANPAGES = 

man_MANS = ${MANPAGES}

EXTRA_DIST = ${MANPAGES}

bin_PROGRAMS = damapper_bwt

damapper_bwt_SOURCES = damapper_bwt.cpp map.c prs_uint64_t_pair.c
damapper_bwt_LDADD = ${LIBMAUS2LIBS} ${DALIGNERLIBS}
damapper_bwt_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${DALIGNERLDFLAGS} ${AM_LDFLAGS}
damapper_bwt_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${DALIGNERCPPFLAGS} -DUSE_PRS
