dnl
dnl Copyright 2010 by Virginia Polytechnic Institute and State
dnl University. All rights reserved. Virginia Polytechnic Institute and
dnl State University (Virginia Tech) owns the software and its
dnl associated documentation.
dnl

AC_PREREQ(2.63)

AC_INIT([opencl-apps], [0])
AC_CONFIG_AUX_DIR([m4])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.10])

PAC_PUSH_ALL_FLAGS
LT_INIT
PAC_POP_ALL_FLAGS

PAC_PROG_CC
PAC_PROG_CXX

AC_CONFIG_HEADER(config.h)

# Strict flags should come first as it might change the way the
# remaining tests behave
PAC_ARG_STRICT

AM_PROG_MKDIR_P
AC_PROG_INSTALL
AC_C_INLINE
AC_C_CONST
AC_C_RESTRICT
PAC_ARG_CACHING
PAC_C_MACRO_VA_ARGS
PAC_C_GNU_ATTRIBUTE
PAC_CC_FUNCTION_NAME_SYMBOL

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


### The OpenCL SDK allows us to build some additional applications

# NOTE: We prioritize the OpenCL SDK compared to the regular OpenCL in
# order to avoid picking the system OpenCL. This is because the OpenCL
# library seems to have a dependency on some scripts that ship with
# the SDK, which is not installed on the system. So if the person who
# installed the system OpenCL library deleted his/her SDK source, we
# will be in trouble.
AC_ARG_WITH([opencl-sdk],
	AC_HELP_STRING([--with-opencl-sdk],[path to the OpenCL SDK]),
	[OPENCL_SDK_PREFIX=$withval],[OPENCL_SDK_PREFIX=])

if test ! -z "$OPENCL_SDK_PREFIX" ; then
   # FIXME: check for when the 32-bit version is required
   PAC_APPEND_FLAG([-I${OPENCL_SDK_PREFIX}/include],[CPPFLAGS])
   PAC_APPEND_FLAG([-I${OPENCL_SDK_PREFIX}/samples/opencl/SDKUtil/include],[CPPFLAGS])

   PAC_APPEND_FLAG([-L${OPENCL_SDK_PREFIX}/lib/x86_64],[LDFLAGS])
   PAC_APPEND_FLAG([-L${OPENCL_SDK_PREFIX}/samples/opencl/SDKUtil/build/debug/x86_64],[LDFLAGS])

   # Use rpath to force the use of these libraries
   PAC_APPEND_FLAG([-Wl,-rpath=${OPENCL_SDK_PREFIX}/lib/x86_64],[LDFLAGS])
   PAC_APPEND_FLAG(
	[-Wl,-rpath=${OPENCL_SDK_PREFIX}/samples/opencl/SDKUtil/build/debug/x86_64],
	[LDFLAGS])
fi
AC_SUBST(OPENCL_SDK_PREFIX)


### We need OpenCL
AC_ARG_WITH([opencl], AC_HELP_STRING([--with-opencl],[path to OpenCL]),
	[OPENCL_PREFIX=$withval],[OPENCL_PREFIX=])
if test ! -z "$OPENCL_PREFIX" ; then
   PAC_APPEND_FLAG([-L${OPENCL_PREFIX}],[LDFLAGS])
   PAC_APPEND_FLAG([-Wl,-rpath=${OPENCL_PREFIX}],[LDFLAGS])
fi


### Check for the libraries that we need
AC_CHECK_LIB(OpenCL,clGetPlatformIDs,[PAC_PREPEND_FLAG([-lOpenCL],[LIBS])],
	AC_MSG_ERROR([OpenCL library not found]))

AC_LANG_PUSH([C++])
PAC_PUSH_FLAG([LIBS])
PAC_APPEND_FLAG([-lSDKUtil],[LIBS])
AC_MSG_CHECKING([for SDKSample::initialize in -lSDKUtil])
AC_TRY_LINK([
	#include <SDKCommon.hpp>
	#include <SDKApplication.hpp>
	#include <SDKCommandArgs.hpp>
	#include <SDKFile.hpp>],
	[streamsdk::Option x;],found_sdkutil=yes,found_sdkutil=no)
AC_MSG_RESULT([${found_sdkutil}])
if test "${found_sdkutil}" = "no" ; then
   PAC_POP_FLAG([LIBS])
fi
AC_LANG_POP([C++])


### Check what applications the user wants us to build
#AC_ARG_WITH([apps],
#	AC_HELP_STRING([--with-apps],[Applications to build (gem,matmult,mattran,swat)]),
#	[app_list=$withval],[app_list=gem,matmult,mattran,swat])
AC_ARG_WITH([apps],
	AC_HELP_STRING([--with-apps],[Applications to build (nw,srad,lud,kmeans,samplecl,crc,cfd,bsort)]),
	[app_list=$withval],[app_list=nw,srad,lud,kmeans,samplecl,crc,cfd,bsort])
app_names="`echo $app_list | sed -e 's/:/ /g' -e 's/,/ /g'`"

#build_gem=false
#build_matmult=false
#build_mattran=false
#build_swat=false
build_nw=false
build_srad=false
build_lud=false
build_kmeans=false
build_samplecl=false
build_crc=false
build_cfd=false
build_bsort=false

selected_apps=""
for app in ${app_names}; do
    case "$app" in
#	 gem)
#		# We don't need to test anything here, as we would
#		# have already aborted if OpenCL wasn't present
#		build_gem=true
#		selected_apps="$selected_apps gem"
#		;;
#	 matmult)
#		if test "${found_sdkutil}" = "yes"; then
#		   build_matmult=true
#		   selected_apps="$selected_apps matmult"
#		fi
#		;;
#	 mattran)
#		if test "${found_sdkutil}" = "yes"; then
#		   build_mattran=true
#		   selected_apps="$selected_apps mattran"
#		fi
#		;;
#	 swat)
#		if test "${found_sdkutil}" = "yes"; then
#		   build_swat=true
#		   selected_apps="$selected_apps swat"
#		fi
#		;;
     cfd)
        build_cfd=true
        selected_apps="$selected_apps cfd"
        ;;
     crc)
        build_crc=true
        selected_apps="$selected_apps crc"
        ;;
     bsort)
        build_bsort=true
        selected_apps="$selected_apps bsort"
        ;;
     nw)
        build_nw=true
        selected_apps="$selected_apps nw"
        ;;
     srad)
        build_srad=true
        selected_apps="$selected_apps srad"
        ;;
     lud)
        build_lud=true
        selected_apps="$selected_apps lud"
        ;;
     kmeans)
        build_kmeans=true
        selected_apps="$selected_apps kmeans"
        ;;
     samplecl)
        build_samplecl=true
        selected_apps="$selected_apps samplecl"
        ;;
	 *)
		AC_MSG_ERROR([unrecognized application $app])
		;;
    esac
done

AC_MSG_RESULT([building applications... ${selected_apps}])
#AM_CONDITIONAL([build_gem], [$build_gem])
#AM_CONDITIONAL([build_matmult], [$build_matmult])
#AM_CONDITIONAL([build_mattran], [$build_mattran])
#AM_CONDITIONAL([build_swat], [$build_swat])
AM_CONDITIONAL([build_nw], [$build_nw])
AM_CONDITIONAL([build_srad], [$build_srad])
AM_CONDITIONAL([build_lud], [$build_lud])
AM_CONDITIONAL([build_kmeans], [$build_kmeans])
AM_CONDITIONAL([build_samplecl], [$build_samplecl])
AM_CONDITIONAL([build_crc], [$build_crc])
AM_CONDITIONAL([build_cfd], [$build_cfd])
AM_CONDITIONAL([build_bsort], [$build_bsort])

# Final output
AC_OUTPUT(Makefile)
