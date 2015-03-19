/*
fba package
Copyright (C) 2013 Daniel Ortiz-Mart\'inez
 
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
*/
 
/**
 * @file getdelim.h
 * @brief Defines the getdelim function
 */

#ifndef _getdelim_h
#define _getdelim_h

#if HAVE_CONFIG_H
#  include <fba_config.h>
#else
#  define FBA_TIME_WITH_SYS_TIME 1
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#if FBA_HAVE_UNISTD_H
# include <unistd.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif
#ifndef FBA_HAVE_GETDELIM
# define GETDELIM_BUFFER 128
  ssize_t getdelim(char **lineptr, size_t *n, int delimiter, FILE *stream);
#endif
#ifdef __cplusplus  
}
#endif

#endif
