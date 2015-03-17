# general settings

LGMP=@LGMP@
ZLIB=@ZLIB@

bin_PROGRAMS = solve_lp_probl
lib_LTLIBRARIES = libfba.la

# flags

# The following variables must be defined externally if the required
# libraries have been installed in local directories:
#
# NLP_COMMON_HOME

AM_CXXFLAGS= #-I$(SW_MODELS_HOME)/include -I$(NLP_COMMON_HOME)/include

LDADD = #-lsw_models -lnlp_common -L$(SW_MODELS_HOME)/lib

# sources

common_h= fba_config.h
common_defs=

# edit_dist_h= BaseEditDist.h _editDist.h EditDistForStr.h	\
# EditDistForVec.h EditDistForVecString.h
# edit_dist_defs= EditDistForStr.cc EditDistForVecString.cc

# programs
solve_lp_probl_SOURCES = solve_lp_probl.cc
solve_lp_probl_LDFLAGS = libfba.la

##########

# include headers
include_HEADERS = $(common_h) 

##########

# libraries
fba_la_SOURCES= $(common_h) $(common_defs)

fba_la_LIBADD= $(LDADD)
