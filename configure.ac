#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([fcap],[1.0],[daniel.ortiz@ub.edu])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_HEADERS([config.h])
AX_PREFIX_CONFIG_H([./src/fcap_config.h])
AM_PATH_PYTHON([3])

dnl Initialize Libtool
LT_CONFIG_LTDL_DIR([libltdl])
LT_INIT([dlopen])
LTDL_INIT([recursive])

# Save previous compilation flags
CFLAGS_AUX=${CFLAGS}
CXXFLAGS_AUX=${CXXFLAGS}

# Checks for Python modules.
AX_PYTHON_MODULE(matplotlib,fatal,python)
AX_PYTHON_MODULE(sklearn,fatal,python)
AX_PYTHON_MODULE(scipy,fatal,python)
AX_PYTHON_MODULE(numpy,fatal,python)
AX_PYTHON_MODULE(statsmodels,fatal,python)
AX_PYTHON_MODULE(json,fatal,python)

# Checks for programs.
AC_PATH_PROGS(PYTHON, python)
AC_PATH_PROGS(RSCRIPT, Rscript)
AC_PATH_PROGS(SPLIT, split)
AC_PATH_PROGS(SED, sed)
AC_PATH_PROGS(GREP, grep)
AC_PATH_PROGS(UNIQ, uniq)
AC_PATH_PROGS(HEAD, head)
AC_PATH_PROGS(TAIL, tail)
AC_PATH_PROGS(FIND, find)
AC_PATH_PROGS(LN, ln)
AC_PATH_PROGS(WC, wc)
AC_PATH_PROGS(BASENAME, basename)
AC_PATH_PROGS(DIRNAME, dirname)
AC_PATH_PROGS(MKTEMP, mktemp)
AC_PATH_PROGS(GZIP, gzip)
AC_PATH_PROGS(SSH, ssh)
AC_PATH_PROGS(QSUB, qsub)
AC_PATH_PROGS(QSTAT, qstat)
AC_PATH_PROGS(QRLS, qrls)

# BASH may already be set in the shell, if the admin then changes the
# the /bin/sh symlink to a non-bash shell, then the installation will
# break
unset BASH
AC_PATH_PROGS(BASH, bash, "no", [$PATH:/bin:/usr/bin:/usr/sbin])
if test x$BASH = "xno"; then
    AC_MSG_ERROR([bash is required])
fi

# check that qsub works
if test ! -z "$QSUB"; then
 echo -n "checking whether qsub works... "
 echo "" | ${QSUB} >/dev/null 2>&1 
 if test $? -eq 0; then
  QSUB_WORKS="yes"
  echo "yes"
 else
  QSUB_WORKS="no"
  echo "no"
 fi
else
 QSUB_WORKS="no"
fi
AC_SUBST(QSUB_WORKS)

# check qsub parameter (-terse)
if test "${QSUB_WORKS}" = "yes"; then
 echo -n "checking whether qsub accepts -terse... "
 echo "" | ${QSUB} -terse >/dev/null 2>&1
 if test $? -eq 0; then
  QSUB_TERSE_OPT="-terse"
  echo "yes"
 else
  QSUB_TERSE_OPT=""
  echo "no"
 fi 
else
 QSUB_TERSE_OPT=""
fi
AC_SUBST(QSUB_TERSE_OPT)

# check qstat parameter (-j)
if test "${QSUB_WORKS}" = "yes"; then
 echo -n "checking whether qstat accepts -j... "
 jid=0
 ${QSTAT} -j ${jid} 2>&1 | $GREP "invalid option" >/dev/null 2>&1
 if test $? -eq 0; then
  QSTAT_J_OPT=""
  echo "no"
 else
  QSTAT_J_OPT="-j"
  echo "yes"
 fi 
else
 QSTAT_J_OPT=""
fi
AC_SUBST(QSTAT_J_OPT)


# Search awk

## GAWK
AC_PATH_PROG(GAWK, gawk, no)
if test "$GAWK" = no; then

   ## NAWK
   AC_PATH_PROG(NAWK, nawk, no)
   if test "$NAWK" = no; then
        
        ## MAWK
	AC_PATH_PROG(MAWK, mawk, no)
        if test "$MAWK" = no; then
	      AWK=awk
        else
              # mawk found
              AWK=$MAWK
        fi
   else
        # nawk found
        AWK=$NAWK
   fi
else
        # gawk found
        AWK=$GAWK
fi

AC_SUBST(AWK)

AC_SUBST(USE_SOURCE_OPT)

# check awk delete array
AC_MSG_CHECKING(whether awk supports delete array )
USE_FAST_DELETE="no"
echo "" | $AWK '{ a[[1]]=1; delete a; if (a[[1]]!=0) exit (1) }' && USE_FAST_DELETE="yes"
if test ${USE_FAST_DELETE} = "yes"; then
   AC_MSG_RESULT(yes)   
else
   AC_MSG_RESULT(no) 
fi

AC_SUBST(USE_FAST_DELETE)

# Search sort and check parameters
AC_PATH_PROGS(SORT,sort)

echo -n "checking whether sort accepts -S <buffsize>... "
echo "" | ${SORT} -S 131072 >/dev/null 2>&1
if test $? -eq 0; then
   sortpars="-S 131072"
   echo "yes"
else
   sortpars=""
   echo "no"
fi 
AC_SUBST(sortpars)

echo -n "checking whether sort accepts -T <dir>... "
echo "" | ${SORT} -T /tmp >/dev/null 2>&1
if test $? -eq 0; then
   sortT="yes"
   echo "yes"
else
   sortT="no"
   echo "no"
fi 
AC_SUBST(sortT)

echo -n "checking whether sort accepts -m... "
echo "" | ${SORT} -m >/dev/null 2>&1
if test $? -eq 0; then
   sortm="yes"
   echo "yes"
else
   sortm="no"
   echo "no"
fi 
AC_SUBST(sortm)


# Check for compilers
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O

CFLAGS="${CFLAGS_AUX}"
CXXFLAGS="${CXXFLAGS_AUX}"

# Checks for basic flags
AC_LANG(C++)
for f in -Wno-deprecated
do                        
    AC_MSG_CHECKING(whether ${CXX} supports ${f})
    saved_cxxflags="${CXXFLAGS}"
    CXXFLAGS="${CXXFLAGS} ${f}"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }],
                      AC_MSG_RESULT(yes),
                      AC_MSG_RESULT(no)
                      CXXFLAGS="${saved_cxxflags}")],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
done

dnl Check if we have enable debug support.
AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE(debug,[  --enable-debug   turn on debugging], enable_debug=yes,
              enable_debug=no)

dnl Yes, debug is enabled
if test "$enable_debug" = "yes"; then

AC_MSG_RESULT(yes)

AC_DEFINE(DEBUG,,[Define if the debugging mode is enabled])

# Checks for flags
AC_LANG(C)
for f in -g -Wall  
do                        
    AC_MSG_CHECKING(whether ${CC} supports ${f})
    saved_cflags="${CFLAGS}"
    CFLAGS="${CFLAGS} ${f}"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }],
                      AC_MSG_RESULT(yes),
                      AC_MSG_RESULT(no)
                      CFLAGS="${saved_cflags}")],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
done

AC_LANG(C++)
for f in -g -Wall
do                        
    AC_MSG_CHECKING(whether ${CXX} supports ${f})
    saved_cxxflags="${CXXFLAGS}"
    CXXFLAGS="${CXXFLAGS} ${f}"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }],
                      AC_MSG_RESULT(yes),
                      AC_MSG_RESULT(no)
                      CXXFLAGS="${saved_cxxflags}")],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
done

else

AC_MSG_RESULT(no)

# Checks for flags
AC_LANG(C)
for f in -Wall -O2 
do                        
    AC_MSG_CHECKING(whether ${CC} supports ${f})
    saved_cflags="${CFLAGS}"
    CFLAGS="${CFLAGS} ${f}"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }],
                      AC_MSG_RESULT(yes),
                      AC_MSG_RESULT(no)
                      CFLAGS="${saved_cflags}")],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
done

AC_LANG(C++)
for f in -Wall -O2
do                        
    AC_MSG_CHECKING(whether ${CXX} supports ${f})
    saved_cxxflags="${CXXFLAGS}"
    CXXFLAGS="${CXXFLAGS} ${f}"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }],
                      AC_MSG_RESULT(yes),
                      AC_MSG_RESULT(no)
                      CXXFLAGS="${saved_cxxflags}")],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
done

fi

for top_builddir in . .. ../.. $ac_auxdir $ac_auxdir/..; do
  test -f $top_builddir/configure && break
done

# Checks for libraries.
# AC_PROG_RANLIB

# Adds an extra argument to configure (profiling)
AC_MSG_CHECKING(whether to enable profiling)
AC_ARG_ENABLE(developer,[  --enable-developer   enable developer features],enable_developer=yes,
              enable_developer=no)
dnl Yes, profiling is enabled
if test "$enable_developer" = "yes"; then
AC_MSG_RESULT(yes)

if test "$enable_debug" = "yes"; then
 flags="-pg"
else
 flags="-g -pg"
fi

AC_LANG(C)
for f in ${flags}
do                        
    AC_MSG_CHECKING(whether ${CC} supports ${f})
    saved_cflags="${CFLAGS}"
    CFLAGS="${CFLAGS} ${f}"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }],
                      AC_MSG_RESULT(yes),
                      AC_MSG_RESULT(no)
                      CFLAGS="${saved_cflags}")],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
done

AC_LANG(C++)
for f in ${flags}
do                        
    AC_MSG_CHECKING(whether ${CXX} supports ${f})
    saved_cxxflags="${CXXFLAGS}"
    CXXFLAGS="${CXXFLAGS} ${f}"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }],
                      AC_MSG_RESULT(yes),
                      AC_MSG_RESULT(no)
                      CXXFLAGS="${saved_cxxflags}")],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
done

else
AC_MSG_RESULT(no)

fi

# Checks if the compiler supports ISO C++ standard library
AC_CXX_HAVE_STD
if test "$ac_cv_cxx_have_std" = no; then
  AC_MSG_ERROR(C++ standard library not present!)
fi

# Checks if the compiler supports the Standard Template Library
AC_CXX_HAVE_STL
if test "$ac_cv_cxx_have_stl" = no; then
  AC_MSG_ERROR(C++ standard template library not present!)
fi

# Checks for lm library
AC_CHECK_LIB([m], [log])

# Checks for glpk library
LGLPK=""
libaux=${LIBS}
AC_CHECK_LIB([glpk], [glp_free])
if test $ac_cv_lib_glpk_glp_free = yes; then
   LGLPK="-lglpk"
   AC_DEFINE(HAVE_GLPK,,[Define if the glpk library is present])
fi
AM_CONDITIONAL([HAVE_GLPK], [test x${ac_cv_lib_glpk_glp_free} = xyes])
AC_SUBST(LGLPK)
LIBS=${libaux}

# Checks for gmp library
LGMP=""
libaux=${LIBS}
AC_CHECK_LIB([gmp], [__gmpz_init])
if test $ac_cv_lib_gmp___gmpz_init = yes; then
   LGMP="-lgmp"
   AC_DEFINE(HAVE_GMP,,[Define if the gmp library is present])
fi
AC_SUBST(LGMP)
LIBS=${libaux}

# Checks for zlib library
ZLIB=""
libaux=${LIBS}
AC_CHECK_LIB([z], [inflate])
if test $ac_cv_lib_z_inflate = yes; then
   ZLIB="-lz"
   AC_DEFINE(HAVE_ZLIB,,[Define if the zlib library is present])
fi
AC_SUBST(ZLIB)
LIBS=${libaux}

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_CHECK_TYPE(ssize_t, int)
AC_TYPE_OFF_T

# Checks for large file support
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Checks for library functions.
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday pow getdelim])

 # Some systems do not supply getline()
AC_MSG_CHECKING([if getline() is supported])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
    [[
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
    ]],[[
    char* buf = 0;
    size_t len = 0;
    FILE* f=NULL;
    getline(&buf,&len,f);
    ]]
  )],
  [AC_DEFINE(HAVE_GETLINE,[],[Define if getline() is supported])
   AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)])

# Checks for doxygen
AC_PATH_PROG(DOXYGEN,doxygen)

# Process options given in shell variables

# Adds an extra argument to configure (enable glpk-based programs)
AC_MSG_CHECKING(whether to enable glpk-based programs)
AC_ARG_ENABLE(glpk-progs,[  --enable-glpk-progs   enable glpk-based programs],enable_glpk_progs=yes,enable_glpk_progs=no)
dnl Yes, glpk-based programs enabled
if test "$enable_glpk_progs" = "yes"; then
 AC_MSG_RESULT(yes)
else
 AC_MSG_RESULT(no)
fi
AM_CONDITIONAL([ENABLE_GLPK_BASED_PROGS], [test x${enable_glpk_progs} = xyes])

# Change path of distribution include directory
includedir=${prefix}/include/fba

#AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([Makefile doc/Makefile m4/Makefile src/Makefile src/common/Makefile src/fba/Makefile utils/Makefile examples/Makefile])
AC_OUTPUT
