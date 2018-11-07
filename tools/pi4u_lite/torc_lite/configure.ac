dnl
dnl  Torc-Lite runtime library
dnl  == Copyright (C) 2011 the Torc-Lite Team
dnl  == Department of Computer Science, University of Ioannina
dnl
dnl This file is part of Torc-Lite.
dnl
dnl Torc-Lite is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl Torc-Lite is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with OMPi; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
dnl
dnl This file must be processed by autoconf

dnl Configure.in file for Torc-Lite

AC_INIT(Torc-Lite, 0.1a,,Torc-Lite)
AC_CONFIG_AUX_DIR(auxdir)
AC_CANONICAL_SYSTEM

AC_CONFIG_HEADER(include/torc_config.h)
AM_INIT_AUTOMAKE

AC_PROG_CC(mpicc)
AC_PROG_F77(mpif90)

AC_LANG_C
AC_PROG_CPP

AC_PROG_MAKE_SET
AC_PROG_RANLIB
AM_PROG_AS
# torc_lite + some fortran support
AC_F77_WRAPPERS	

#=============================================================================
# Adjust CPPFLAGS.
#=============================================================================
CPPFLAGS=${CPPFLAGS}' -D_REENTRANT  '
AC_MSG_NOTICE([CPPFLAGS=${CPPFLAGS}])

#=============================================================================
# Adjust CPPFLAGS.
#=============================================================================
CPPFLAGS=${CPPFLAGS}' -D_REENTRANT  '
AC_MSG_NOTICE([CPPFLAGS=${CPPFLAGS}])


#=============================================================================
# Check pthread.h
#=============================================================================
AC_CHECK_HEADERS(pthread.h)

#=============================================================================
# Setup ps_config.h
#=============================================================================
AH_TOP(
[
#ifndef _PS_CONFIG
#define _PS_CONFIG
])

AH_BOTTOM([#endif])

#=============================================================================
# Check if -lrt is needed for we sched_yield()
#=============================================================================
needRT=""
AC_CHECK_FUNCS(sched_yield, [ break ],     
	[ AC_CHECK_LIB(rt, sched_yield, [ needRT=yes ],
	[ AC_DEFINE_UNQUOTED(sched_yield(), 0, not found) ])
	])

if test ".$needRT" = ".yes"; then
	LDFLAGS="${LDFLAGS} -lrt"
fi

#=============================================================================
# Set synchronization method
#=============================================================================
syncopt=""
mutexspin=""
AC_ARG_WITH(sync,
	AC_HELP_STRING([--with-sync=method], [posix synchronization mechanism (mutex, mutex_try, spin, spin_try)]),
	[
	case $withval in
		mutex| mutex_try )	syncopt=$withval 
					mutexspin="mutex" ;;
		spin | spin_try)	syncopt=$withval 
					mutexspin="spin" ;;
		* ) AC_MSG_ERROR([invalid synchronization mechanism -- allowed:mutex,mutex_try,spin,spin_try]) ;;
	esac
	]
)


if test ".$syncopt" = "."; then
	AC_MSG_NOTICE([posix synchronization method not defined - using default (mutex_try)])
	syncopt=mutex_try
	mutexspin=mutex
fi


if test ".$mutexspin" = ".mutex"; then
	# Check for mutex locks
	AC_MSG_CHECKING([for usable pthread_mutex_t])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
	[[	#include <pthread.h>
		pthread_mutex_t lock;
	]],
	[[	pthread_mutex_init(&lock,NULL); 
		pthread_mutex_lock(&lock);
		pthread_mutex_unlock(&lock); 
		pthread_mutex_trylock(&lock); 
	]])],,
	AC_MSG_FAILURE([pthread_mutex_t not available]))
	AC_MSG_RESULT([yes])

elif test ".$mutexspin" = ".spin"; then
	# Check for spin locks
	# Usually they are available with the _XOPEN_SOURCE macro
	AC_MSG_CHECKING([for usable pthread_spinlock_t])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
	[[
	#define _XOPEN_SOURCE 600
	
	#include <pthread.h>
	pthread_spinlock_t lock;
	]],
	[[
	pthread_spin_init(&lock,0); 
	pthread_spin_lock(&lock); 
	pthread_spin_unlock(&lock); 
	pthread_spin_trylock(&lock);
	]])],,
	AC_MSG_FAILURE([pthread_spinlock_t not available]))
	AC_MSG_RESULT([yes])
else
	AC_MSG_FAILURE([no appropriate synchronization method found])
fi

# 
# 

if test ".$syncopt" = ".mutex"; then
	(cd include; echo "#define POSIX_MUTEX_LOCK" > ps_config.h)
elif test ".$syncopt" = ".mutex_try"; then
	(cd include; echo "#define POSIX_MUTEX_TRYLOCK" > ps_config.h)
elif test ".$syncopt" = ".spin"; then
	(cd include; echo  "#define POSIX_SPIN_LOCK" > ps_config.h)
elif test ".$syncopt" = ".spin_try"; then
	(cd include; echo  "#define POSIX_SPIN_TRYLOCK" > ps_config.h)
else
	AC_MSG_FAILURE([no appropriate synchronization method found])
fi


#==============================================================================
# Set max number of nodes (mpi processes)
#==============================================================================

NNODES="64"
AC_ARG_WITH(maxnodes,
	AC_HELP_STRING([--with-maxnodes=num], [maximum number of nodes (default: 64)]),
	[
	NNODES=$withval
	]
)
AC_SUBST(NNODES)

#==============================================================================
# Set max number of virtual processors
#==============================================================================
NVPS="64"
AC_ARG_WITH(maxvps,
	AC_HELP_STRING([--with-maxvps=num], [maximum number of virtual processors (default: 64)]),
	[
	NVPS=$withval
	]
)
AC_SUBST(NVPS)

#==============================================================================
# Enable debugging
#==============================================================================
AC_ARG_ENABLE(debug,
	AC_HELP_STRING([--enable-debug], [enable debugging (default:disabled)]),
	if eval "test x$enable_debug = xyes"; then
		DEBUGFLAG="-g"
	else
		DEBUGFLAG=""
	fi,
	DEBUGFLAG="")

AC_SUBST(DEBUG_FLAG)

#=============================================================================
# Check default cache line size
#=============================================================================
defcls=""

AC_MSG_CHECKING([for default cache line size])
AC_RUN_IFELSE([
AC_LANG_PROGRAM([
/* Borrowed from OMPi (lib/sysdeps.h) */
#include <stdio.h>

/*
 * Processors
 */
#if defined(pentium4) || defined(__pentium4)
  #define __cpu_x86 686
#elif defined(i686) || defined(__i686) || \
    defined(pentiumpro) || defined(__pentiumpro) || \
    defined(pentium2) || defined(__pentium2) || \
    defined(pentium3) || defined(__pentium3)
  #define __cpu_x86 686
  #define __cpu_686_class p6       /* Only used for CACHE_LINE size */
#elif defined(i586) || defined(__i586) || defined(pentium) || defined(__pentium)
  #define __cpu_x86 586
#elif defined(i486) || defined(__i486)
  #define __cpu_x86 486
#elif defined(i386) || defined(__i386)
  #define __cpu_x86 386
#endif

#if defined(sparcv9) || defined(__sparcv9)
  #define __cpu_sparc v9
#elif defined(sparc) || defined(__sparc)
  #define __cpu_sparc v8
#endif

#if defined(mips) || defined(__mips)
  #define __cpu_mips 1
#endif

/*
 * CACHE_LINE size (for L2 cache).
 * A generally good value here is 128 bytes, so it is the default, to be
 * on the safe side.
 * However, different processors have different sizes, so below
 * we try to match their values.
 */
#if defined(__cpu_x86)
  /* Up to and including P6 cpus (Pentium Pro, Pentium II, Pentium III),
   * the L1 & L2 line size is 32 bytes. All newer CPUs have cache lines
   * of 64 bytes. Netbutst CPUs (Pentium 4, Pentium D, Pentium Extreme and
   * some Xeons), in addition have a "sector" size of 128 bytes.
   * #if __cpu_686_class == p6
   *   #define CACHE_LINE 32
   * #else
   *   #define CACHE_LINE 64
   * #endif
   * All known Athlons have 64-byte cache lines.
   */
  #define CACHE_LINE 64
#elif defined(__cpu_sparc)
  /* UltraSparcs (v9) use L2 cache lines of 64 bytes, up to IIIi.
   * IIIcu has line sizes 64 - 512 bytes (!) but anyways, they have
   * the 64 bytes as "a unit of fill & coherency".
   * IV ahs 64 - 128 bytes
   * IV+ and T1 came back to 64-byte cache lines.
   *
   * So, 64 bytes should be ok.
   */
  #define CACHE_LINE 64
#elif defined(__cpu_mips)
   /* The specifications say that the R4000 has a variable size for
    * its cache lines, between 4 and 32 words. MIPS 10000 (MIPS IV)
    * has 16-32 words lines, and it is a 64-bit CPU, which means
    * that it has lines of size 128-256 bytes.
    * The following may need to be modified to 256(?).
    */
  #define CACHE_LINE 128
#else
  #define CACHE_LINE 128       /* A default value */
#endif


int main(int argc, char *argv[]) {
	FILE *fp;
	
	if ((fp = fopen("conftestval", "w")) == NULL)
		exit(1);
	fprintf(fp, "%d\n", CACHE_LINE);
	fclose(fp);

	/* die successfully */
	exit(0);
}
])],
[defcls=`cat conftestval`], [defcls="128"], [defcls="128"])

AC_MSG_RESULT([$defcls])

#=============================================================================
# Set cache line size
#=============================================================================
clsize=""
AC_ARG_WITH(cachelinesize,
	AC_HELP_STRING([--with-cachelinesize=value], [cache line size]),
	[
		clsize=$withval
		AC_MSG_NOTICE([cache line size is set to $clsize])
	]
)

if test ".$clsize" = "."; then
	clsize=$defcls
	AC_MSG_NOTICE([using default cache line size ($clsize bytes) ])
fi

AC_DEFINE_UNQUOTED(CACHE_LINE_SIZE, $clsize, [explicit cache line size])


(cd include; echo "#define CACHE_LINE_SIZE $clsize" >> ps_config.h)

#=============================================================================
# Create Torc-Lite_cflags, Torc-Lite_libs
#=============================================================================
instdir=""
if test ".$prefix" = ".NONE"; then
instdir="/usr/local"
else
instdir=$prefix
fi

if test ".$mutexspin" = ".spin"; then
(cd scripts; echo "echo -I${instdir}/include -D_REENTRANT -D_XOPEN_SOURCE=600" > torc_cflags; chmod +x torc_cflags)
else
(cd scripts; echo "echo -I${instdir}/include -D_REENTRANT" > torc_cflags; chmod +x torc_cflags)
fi
(cd scripts; echo "echo -L${instdir}/lib -ltorc -lpthread" > torc_libs; chmod +x torc_libs)

#==============================================================================
# Torc-Lite 
# ==============================================================================

# MPI installation path
AC_ARG_WITH(mpi, 
	[  --with-mpi=<dir>        Location of the MPI installation],
	if test x$withval = xyes; then
		AC_MSG_ERROR(--with-mpi must be given a pathname)
	else
		CFLAGS="${CFLAGS} -I${withval}/include"
		FFLAGS="${FFLAGS} -I${withval}/include"
		LDFLAGS="${LDFLAGS} -L${withval}/lib"
		MPILIB="-lmpi"
	fi)

# MPI library name
AC_ARG_WITH(mpilib, 
	[  --with-mpilib=<mpilib>        Name of MPI library],
	if test x$withval = xyes; then
		AC_MSG_WARN(--with-mpilib must be given a library. Using default: mpich)
		MPILIB="-lmpich"
	else
		MPILIB="-l${withval}"
	fi)

# MPI include directory
AC_ARG_WITH(mpiincdir, 
	[  --with-mpiincdir=<mpiincdir>        Directory of MPI header files],
	if test x$withval = xyes; then
		AC_MSG_ERROR(--with-mpiincdir must be given a directory)
	else
		CFLAGS="${CFLAGS} -I${withval}"
		FFLAGS="${FFLAGS} -I${withval}"
	fi)

# MPI library directory
AC_ARG_WITH(mpilibdir, 
	[  --with-mpilibdir=<mpilibdir>        Directory of MPI library],
	if test x$withval = xyes; then
		AC_MSG_ERROR(--with-mpilibdir must be given a directory)
	else
		LDFLAGS="${LDFLAGS} -L${withval}"
	fi)

# Statistics
#AC_ARG_ENABLE(stats,
#	AC_HELP_STRING([--enable-stats], [enable statistics (default:disabled)]),
#	if eval "test x$enable_stats = xyes"; then
		AC_DEFINE(TORC_STATS,,[Runtime statistics])
#	fi,)

# Non blocking communication (useful for mpi libraries that do not provide thread safety) 
#AC_ARG_ENABLE(nonblocking,
#	AC_HELP_STRING([--enable-nonblocking], [enable nonblocking server thread (default:disabled)]),
#	if eval "test x$enable_nonblocking = xyes"; then
#		AC_DEFINE(NONBLOCKING,,[Non-blocking server thread])
#	fi,)

###########################################################################

#
#
#
dnl AC_CHECK_DECL(__GNUC__, [ echo yes ], [ echo no])

AC_ARG_ENABLE(debug,
	AC_HELP_STRING([--enable-debug], [enable debugging (default:disabled)]),
	if eval "test x$enable_debug = xyes"; then
		DEBUGFLAG="-g"
	else
		DEBUGFLAG=""
	fi,
	DEBUGFLAG="")


LDFLAGS="${LDFLAGS} ${MPILIB}"
LDFLAGS="${LDFLAGS} -lpthread"

echo "Making with LDFLAGS=$LDFLAGS"

AC_SUBST(MPICC)
AC_SUBST(MPILIB)

# Check for MPI
AC_MSG_CHECKING([for MPI/C support])
mpic_enabled=""
AC_LANG_PUSH(C)
AC_LINK_IFELSE([
AC_LANG_PROGRAM([
	#include <stdio.h>
	#include <mpi.h>
	int main(int argc, char **argv) {
		printf("hello\n");    
		MPI_Init(&argc, &argv);
		MPI_Finalize();
		return 0;
	}
]
)
],
[mpic_enabled=yes], [mpic_enabled=no])
AC_LANG_POP(C)

#echo $mpic_enabled

if test ".$mpic_enabled" = ".yes"; then
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
#	AC_MSG_FAILURE([MPI/C not enabled])
fi                 

AC_OUTPUT(Makefile src/Makefile scripts/Makefile demo/Makefile)
