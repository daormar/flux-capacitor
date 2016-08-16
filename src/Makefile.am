# general settings

LGLPK=@LGLPK@
LGMP=@LGMP@
ZLIB=@ZLIB@

AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS= common fba

if ENABLE_GLPK_BASED_PROGS
if HAVE_GLPK
GLPK_BASED_PROGS=solve_lp_probl
endif
endif

bin_PROGRAMS = $(GLPK_BASED_PROGS)
lib_LTLIBRARIES = libfcap_common.la # libcap.la

# flags

AM_CXXFLAGS= -I$(srcdir)/common -I$(srcdir)/fba

LDADD = $(LGLPK) $(LGMP) $(ZLIB)

# sources

common_h= fcap_config.h common/ErrorDefs.h common/options.h		\
common/mem_alloc_utils.h common/awkInputStream.h common/getline.h	\
common/getdelim.h common/ctimer.h
common_defs= common/options.cc common/mem_alloc_utils.cc	\
common/awkInputStream.cc common/getline.c common/getdelim.c	\
common/ctimer.c

# programs
solve_lp_probl_SOURCES = fba/solve_lp_probl.cc fba/solve_lp_probl_pars.h
solve_lp_probl_LDFLAGS = libfba_common.la # libfba.la

##########

# include headers
include_HEADERS = $(common_h) 

##########

# libraries
libfcap_common_la_SOURCES= $(common_h) $(common_defs)
libfcap_common_la_LIBADD= $(LDADD)

# libfcap_la_SOURCES= $(common_h) $(common_defs)
# libfcap_la_LIBADD= $(LDADD)
