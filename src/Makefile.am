# general settings

LGLPK=@LGLPK@
LGMP=@LGMP@
ZLIB=@ZLIB@

AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS= common fba

bin_PROGRAMS = solve_lp_probl
lib_LTLIBRARIES = libfba_common.la # libfba.la

# flags

AM_CXXFLAGS= -I$(srcdir)/common -I$(srcdir)/fba

LDADD = $(LGLPK) $(LGMP) $(ZLIB) #-lsw_models -lnlp_common -L$(SW_MODELS_HOME)/lib

# sources

common_h= fba_config.h common/ErrorDefs.h common/options.h		\
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
libfba_common_la_SOURCES= $(common_h) $(common_defs)
libfba_common_la_LIBADD= $(LDADD)

# libfba_la_SOURCES= $(common_h) $(common_defs)
# libfba_la_LIBADD= $(LDADD)
