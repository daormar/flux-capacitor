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
 
#ifndef _solve_lp_probl_pars
#define _solve_lp_probl_pars

//--------------- Include files --------------------------------------

#if HAVE_CONFIG_H
#  include <fba_config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <options.h>

using namespace std;

//--------------- Structs --------------------------------------------

struct solve_lp_probl_pars
{
  bool f_given;
  std::string f_str;
  bool o_given;
  std::string o_str;
  bool v_given;
  bool v1_given;

  solve_lp_probl_pars()
    {
      default_values();
    }

  void default_values(void)
    {
      f_given=false;
      o_given=false;
      v_given=false;
      v1_given=false;      
    }
};

#endif
