AC_INIT(damapper_bwt,0.0.0,[tischler@mpi-cbg.de],[damapper_bwt],[http://www.mpi-cbg.de])
AC_CANONICAL_SYSTEM
AC_PROG_LIBTOOL

AM_INIT_AUTOMAKE([std-options subdir-objects parallel-tests])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_SYSTEM

AC_ARG_ENABLE(optimization,
        AS_HELP_STRING([--enable-optimization],[use compiler optimization (default yes)]),
        [optimization=${enableval}],[optimization=yes])
AC_ARG_ENABLE(debug,
        AS_HELP_STRING([--enable-debug],[use compiler debug flags (default no)]),
        [debug=${enableval}],[debug=no])
AC_ARG_ENABLE(profile,
        AS_HELP_STRING([--enable-profile],[use compiler profiling flags (default no)]),
        [profile=${enableval}],[profile=no])
AC_ARG_ENABLE(fast,
        AS_HELP_STRING([--enable-fast],[disable evaluation of assertions (default no)]),
        [fast=${enableval}],[fast=no])

# check for C compiler
CFLAGS_SAVE="${CFLAGS}"
CFLAGS=""
AC_PROG_CC
CFLAGS="${CFLAGS_SAVE}"

CFLAGS_SAVE="${CFLAGS}"
CFLAGS=""
CXXFLAGS_SAVE="${CXXFLAGS}"
CXXFLAGS=""
AM_PROG_CC_C_O
CFLAGS="${CFLAGS_SAVE}"
CXXFLAGS="${CXXFLAGS_SAVE}"

# check for C++ compiler
CXXFLAGS_SAVE="${CXXFLAGS}"
CXXFLAGS=""
AC_PROG_CXX
CXXFLAGS="${CXXFLAGS_SAVE}"

if test "${debug}" = "yes" ; then
	CXXDEBUGFLAGS="-g -O0 -rdynamic"
	CDEBUGFLAGS="-g -O0 -rdynamic"
        CXXFLAGS="${CXXFLAGS} ${CXXDEBUGFLAGS}"
        CFLAGS="${CFLAGS} ${CDEBUGFLAGS}"
	openmp="no"

	AC_MSG_NOTICE([Adding debug options ${CDEBUGFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])
	AC_MSG_NOTICE([Adding debug options ${CXXDEBUGFLAGS} to C compiler ${CXX} flags, full list now ${CXXFLAGS}])
else
	if test "${profile}" = "yes" ; then
	        CXXPROFFLAGS="-g -pg -rdynamic"
	        CPROFFLAGS="-g -pg -rdynamic"

	        CXXFLAGS="${CXXPROFFLAGS} -g -pg -rdynamic"
	        CFLAGS="${CPROFFLAGS} -g -pg -rdynamic"

		AC_MSG_NOTICE([Adding profiling options ${CPROFFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])
		AC_MSG_NOTICE([Adding profiling options ${CXXPROFFLAGS} to C compiler ${CXX} flags, full list now ${CXXFLAGS}])
	else
		if test "${optimization}" = "yes" ; then
			COPTFLAGS=
			case ${CC} in
				gcc)
					COPTFLAGS="-O3 -rdynamic"
					;;
				*-mingw32msvc-gcc)
					COPTFLAGS="-O3 -rdynamic"
					;;
				cl.exe)
					COPTFLAGS="-O2 -Ob2 -Ot -Oy"
					;;
			esac
			CFLAGS="${CFLAGS} ${COPTFLAGS}"
			AC_MSG_NOTICE([Adding optimisation options ${COPTFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])

			CXXOPTFLAGS=
			case ${CXX} in
				g++)
					CXXOPTFLAGS="-O3 -rdynamic"
					;;
				*-mingw32msvc-g++)
					CXXOPTFLAGS="-O3 -rdynamic"
					;;
				cl.exe)
					CXXOPTFLAGS="-O2 -Ob2 -Ot -Oy"
					;;
			esac
			CXXFLAGS="${CXXFLAGS} ${CXXOPTFLAGS}"

			AC_MSG_NOTICE([Adding optimisation options ${CXXOPTFLAGS} to C++ compiler ${CXX} flags, full list now ${CXXFLAGS}])
		fi
	fi
fi


# check for pkg-config tool
PKG_PROG_PKG_CONFIG

if test "${fast}" = "yes" ; then
	CPPFLAGS="${CPPFLAGS} -DNDEBUG"
fi

# check for standard headers
AC_LANG_PUSH([C++])
AC_HEADER_STDC
AC_LANG_POP

if test "${CXX}" = "g++" ; then
	CPPFLAGS="${CPPFLAGS} -W -Wall"
fi

AC_ARG_WITH([libmaus2],
            [AS_HELP_STRING([--with-libmaus2@<:@=PATH@:>@], [path to installed libmaus2 library @<:@default=@:>@])],
            [with_libmaus2=${withval}],
            [with_libmaus2=])

if test ! -z "${with_libmaus2}" ; then
	PKGCONFIGPATHSAVE="${PKG_CONFIG_PATH}"
	if test -z "${PKG_CONFIG_PATH}" ; then
		export PKG_CONFIG_PATH="${with_libmaus2}/lib/pkgconfig"
	else
		export PKG_CONFIG_PATH="${with_libmaus2}/lib/pkgconfig:${PKG_CONFIG_PATH}"
	fi
fi

PKG_CHECK_MODULES([libmaus2],[libmaus2 >= 2.0.239])

if test ! -z "${with_libmaus2}" ; then
	if test ! -z "${PKGCONFIGPATHSAVE}" ; then
		export PKG_CONFIG_PATH="${PKGCONFIGPATHSAVE}"
	fi
fi

LIBMAUS2CPPFLAGS="${libmaus2_CFLAGS} ${libmaus2seqchksumsfactory_CFLAGS} ${libmaus2digests_CFLAGS}"
LIBMAUS2LIBS="${libmaus2_LIBS} ${libmaus2seqchksumsfactory_LIBS} ${libmaus2digests_LIBS}"

CPPFLAGS_SAVE="${CPPFLAGS}"
LDFLAGS_SAVE="${LDFLAGS}"
LIBS_SAVE="${LIBS}"

CPPFLAGS="${CPPFLAGS} ${libmaus2_CFLAGS}"
LIBS="${LIBS} ${libmaus2_LIBS}"

AC_LANG_PUSH([C++])
AC_MSG_CHECKING([whether we can compile a program using libmaus2])
	AC_LINK_IFELSE([AC_LANG_SOURCE([
#include <libmaus2/util/ArgInfo.hpp>

				int main(int argc, char * argv[[]]) {
					::libmaus2::util::ArgInfo const arginfo(argc,argv);
					return 0;
			}])],
			have_libmaus2=yes,
			have_libmaus2=no
		)
AC_LANG_POP
AC_MSG_RESULT($have_libmaus2)

if test "${have_libmaus2}" = "no" ; then
	AC_MSG_ERROR([Required libmaus2 is not available.])
fi

AC_ARG_WITH([daligner],
	[AS_HELP_STRING([--with-daligner@<:@=PATH@:>@], [path to installed daligner align library @<:@default=no@:>@])],
	[with_daligner=${withval}],
	[with_daligner=no])
DABWT_HAVE_DALIGNER=
DALIGNERCPPFLAGS=
DALIGNERLDFLAGS=
DALIGNERLIBS=


if test "${with_daligner}" != "no" ; then
	DALIGNERDIR=
	if test "${with_daligner}" != "yes" ; then
		DALIGNERDIR=${with_daligner}
	fi
	
	CPPFLAGS_SAVE=${CPPFLAGS}
	CXXFLAGS_SAVE=${CXXFLAGS}
	LIBS_SAVE=${LIBS}
	LDFLAGS_SAVE=${LDFLAGS}
	
	if test -z "${DALIGNERDIR}" ; then
		DALIGNERCPPFLAGS=
		DALIGNERLDFLAGS=
		DALIGNERLIBS="-lalign"
		
		CPPFLAGS="${CPPFLAGS} ${DALIGNERCPPFLAGS}"
		LDFLAGS="${LDFLAGS} ${DALIGNERLDFLAGS}"
		LIBS="${LIBS} ${DALIGNERLIBS}"

                AC_LANG_PUSH([C++])
                AC_CHECK_HEADER(align.h, [daligner_align_h=yes], [daligner_align_h=no])
                AC_LANG_POP

		if test "${daligner_align_h}" = "yes" ; then
			AC_MSG_CHECKING([for align.h])
			AC_LANG_PUSH([C++])
			AC_TRY_LINK([#include <align.h>
				],[
				return 0;
			],[daligner_align_h=yes],[daligner_align_h=no])
			AC_MSG_RESULT([${daligner_align_h}])	
                fi
	else
		DALIGNERCPPFLAGS="-I${DALIGNERDIR}/include"
		DALIGNERLDFLAGS="-L${DALIGNERDIR}/lib -Wl,-rpath=${DALIGNERDIR}/lib"
		DALIGNERLIBS="-lalign"
		
		CPPFLAGS="${CPPFLAGS} ${DALIGNERCPPFLAGS}"
		LDFLAGS="${LDFLAGS} ${DALIGNERLDFLAGS}"
		LIBS="${LIBS} ${DALIGNERLIBS}"
		
		AC_LANG_PUSH([C++])
		AC_CHECK_HEADER(align.h, [daligner_align_h=yes], [daligner_align_h=no])
		AC_LANG_POP

		if test "${daligner_align_h}" = "yes" ; then
			AC_MSG_CHECKING([for align.h])
			AC_LANG_PUSH([C++])
			AC_TRY_LINK([#include <align.h>
				],[
				return 0;
			],[daligner_align_h=yes],[daligner_align_h=no])
			AC_MSG_RESULT([${daligner_align_h}])	
                fi

		if test ${daligner_align_h} = "no" ; then
			AS_UNSET(ac_cv_header_align_h)
			
			CPPFLAGS=${CPPFLAGS_SAVE}
			CXXFLAGS=${CXXFLAGS_SAVE}
			LIBS=${LIBS_SAVE}
			LDFLAGS=${LDFLAGS_SAVE}

			DALIGNERCPPFLAGS="-I${DALIGNERDIR}/include"
			DALIGNERLDFLAGS="-L${DALIGNERDIR}/lib"
			DALIGNERLIBS="-lalign"
		
			CPPFLAGS="${CPPFLAGS} ${DALIGNERCPPFLAGS}"
			LDFLAGS="${LDFLAGS} ${DALIGNERLDFLAGS}"
			LIBS="${LIBS} ${DALIGNERLIBS}"
		
			AC_LANG_PUSH([C++])
			AC_CHECK_HEADER(align.h, [daligner_align_h=yes], [daligner_align_h=no])
			AC_LANG_POP
			if test "${daligner_align_h}" = "yes" ; then
				AC_MSG_CHECKING([for align.h])
				AC_LANG_PUSH([C++])
				AC_TRY_LINK([#include <align.h>
					],[
					return 0;
				],[daligner_align_h=yes],[daligner_align_h=no])
				AC_MSG_RESULT([${daligner_align_h}])	
			fi
		fi

		if test ${daligner_align_h} = "no" ; then
			AS_UNSET(ac_cv_header_align_h)
			
			CPPFLAGS=${CPPFLAGS_SAVE}
			CXXFLAGS=${CXXFLAGS_SAVE}
			LIBS=${LIBS_SAVE}
			LDFLAGS=${LDFLAGS_SAVE}

			DALIGNERCPPFLAGS="-I${DALIGNERDIR}"
			DALIGNERLDFLAGS="-L${DALIGNERDIR}"
			DALIGNERLIBS="-lalign"
		
			CPPFLAGS="${CPPFLAGS} ${DALIGNERCPPFLAGS}"
			LDFLAGS="${LDFLAGS} ${DALIGNERLDFLAGS}"
			LIBS="${LIBS} ${DALIGNERLIBS}"
		
			AC_LANG_PUSH([C++])
			AC_CHECK_HEADER(align.h, [daligner_align_h=yes], [daligner_align_h=no])
			AC_LANG_POP
			if test "${daligner_align_h}" = "yes" ; then
				AC_MSG_CHECKING([for align.h])
				AC_LANG_PUSH([C++])
				AC_TRY_LINK([#include <align.h>
					],[
					return 0;
				],[daligner_align_h=yes],[daligner_align_h=no])
				AC_MSG_RESULT([${daligner_align_h}])	
			fi
		fi
	fi
	
	
	if test "${daligner_align_h}" = "yes" ; then	
AC_MSG_CHECKING([whether we can compile a program using the daligner library])
AC_LANG_PUSH([C++])
	AC_TRY_LINK([
#include <align.h>
#include <cstdlib>
#include <cstring>

struct DalignerWorkData
{
	Work_Data * data;
	
	DalignerWorkData() : data(NULL)
	{
		data = New_Work_Data();

		if ( ! data )
		{
			abort();
		}
	}
	
	~DalignerWorkData()
	{
		Free_Work_Data(data);
	}
};
],[
	char t[] = "ACGTACGT";
	char q[] = "ACGTTACGT";
	size_t const tn = strlen(t);
	size_t const qn = strlen(q);
	Path path;
	path.trace = NULL;
	path.tlen = 0;
	path.diffs = 40;
	path.abpos = 0;
	path.bbpos = 0;
	path.aepos = tn;
	path.bepos = qn;

	Alignment algn;
	algn.path = &path;
	algn.flags = 0;
	algn.aseq = t;
	algn.bseq = q;
	algn.alen = tn;
	algn.blen = qn;
	
	DalignerWorkData DWD;

	Compute_Trace_ALL(&algn,DWD.data);

	float freq[] = { 0.25, 0.25, 0.25, 0.25 };
	Align_Spec * spec = New_Align_Spec(0.75 /* correlation */, 100 /* tspace */, &freq[0]);
	Work_Data * workdata = New_Work_Data();
	Alignment align;

	Path * la_path = Local_Alignment(&align,workdata,spec,0,0,0,-1,-1);

	Free_Work_Data(workdata);
	Free_Align_Spec(spec);

	return 0;
],[daligner=yes],[daligner=no])
AC_MSG_RESULT([${daligner}])
		if test "${daligner}" = "yes" ; then
			DABWT_HAVE_DALIGNER="#define DABWT_HAVE_DALIGNER"
		fi
	fi

	CPPFLAGS=${CPPFLAGS_SAVE}
	CXXFLAGS=${CXXFLAGS_SAVE}
	LIBS=${LIBS_SAVE}
	LDFLAGS=${LDFLAGS_SAVE}
fi

if test -z "${DABWT_HAVE_DALIGNER}" ; then
	AC_MSG_ERROR([Required daligner library is not available.])
fi

AC_MSG_NOTICE([Using flags ${CFLAGS} for C compiler ${CC}])
AC_MSG_NOTICE([Using flags ${CXXFLAGS} for C++ compiler ${CXX}])

PACKAGE_NAME=${PACKAGE}
PACKAGE_VERSION=${VERSION}

AC_SUBST([DALIGNERCPPFLAGS])
AC_SUBST([DALIGNERLDFLAGS])
AC_SUBST([DALIGNERLIBS])

AC_SUBST([PACKAGE_NAME])
AC_SUBST([PACKAGE_VERSION])
# libmaus2
AC_SUBST([LIBMAUS2CPPFLAGS])
AC_SUBST([LIBMAUS2LIBS])
# 
AC_OUTPUT(Makefile src/Makefile)
