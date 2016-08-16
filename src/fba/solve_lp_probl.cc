/*
fcap package
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
 * @file solve_lp_probl.cc
 * @brief Utility to solve linear programming problems.
 */

//--------------- Include files ---------------------------------------

#include <options.h>
#include <ctimer.h>
#include <set>
#include <string>
#include <ErrorDefs.h>
#include "solve_lp_probl_pars.h"
#include <glpk.h>

//--------------- Constants -------------------------------------------


//--------------- Function Declarations -------------------------------

int processParameters(solve_lp_probl_pars pars);
int handleParameters(int argc,
                     char *argv[],
                     solve_lp_probl_pars& pars);
int takeParameters(int argc,
                   const vector<std::string>& argv_stl,
                   solve_lp_probl_pars& pars);
int checkParameters(solve_lp_probl_pars& pars);
void printParameters(solve_lp_probl_pars pars);
void printUsage(void);
void version(void);

//--------------- Constants -------------------------------------------


//--------------- Global variables ------------------------------------


//--------------- Function Definitions --------------------------------


//--------------- main function
int main(int argc,char *argv[])
{
  solve_lp_probl_pars pars;
    
  if(handleParameters(argc,argv,pars)==ERROR)
  {
    return ERROR;
  }
  else
  {
    return processParameters(pars);
  }
}

//--------------- processParameters function
int processParameters(solve_lp_probl_pars pars)
{
      // Create glp_prob
  glp_prob *glpProbPtr=glp_create_prob();

      // Read lp problem from file in CPLEX format
  glp_cpxcp *parmPtr=NULL;
  
  int ret=glp_read_lp(glpProbPtr,parmPtr,pars.f_str.c_str());
  if(ret)
  {
        // If ret not equal to zero, release memory and terminate
    cerr<<"Error while reading problem in CPLEX format\n";
    glp_delete_prob(glpProbPtr);
    return ERROR;
  }

      // Solve problem
  ret=glp_simplex(glpProbPtr,NULL);
  if(ret)
  {
    cerr<<"Warning: problem not sucessfully solved\n";
  }

      // Print solution
  ret=glp_print_sol(glpProbPtr,pars.o_str.c_str());
  
      // Delete glp_prob
  glp_delete_prob(glpProbPtr);

  return OK;
}

//--------------- handleParameters function
int handleParameters(int argc,
                     char *argv[],
                     solve_lp_probl_pars& pars)
{
  if(argc==1 || readOption(argc,argv,"--version")!=-1)
  {
    version();
    return ERROR;
  }
  if(readOption(argc,argv,"--help")!=-1)
  {
    printUsage();
    return ERROR;   
  }

  vector<std::string> argv_stl=argv2argv_stl(argc,argv);
  if(takeParameters(argc,argv_stl,pars)==ERROR)
  {
    return ERROR;
  }
  else
  {
    if(checkParameters(pars)==OK)
    {
      printParameters(pars);
      return OK;
    }
    else
    {
      return ERROR;
    }
  }
}

//--------------- takeparameters function
int takeParameters(int argc,
                   const vector<std::string>& argv_stl,
                   solve_lp_probl_pars& pars)
{
  int i=1;
  unsigned int matched;
  
  while(i<argc)
  {
    matched=0;
    
        // -w parameter
    if(argv_stl[i]=="-f" && !matched)
    {
      pars.f_given=true;
      if(i==argc-1)
      {
        cerr<<"Error: no value for -f parameter."<<endl;
        return ERROR;
      }
      else
      {
        pars.f_str=argv_stl[i+1];
        ++matched;
        ++i;
      }
    }

        // -o parameter
    if(argv_stl[i]=="-o" && !matched)
    {
      pars.o_given=true;
      if(i==argc-1)
      {
        cerr<<"Error: no value for -o parameter."<<endl;
        return ERROR;
      }
      else
      {
        pars.o_str=argv_stl[i+1];
        ++matched;
        ++i;
      }
    }

        // -v parameter
    if(argv_stl[i]=="-v" && !matched)
    {
      pars.v_given=true;
      ++matched;
    }

        // -v1 parameter
    if(argv_stl[i]=="-v1" && !matched)
    {
      pars.v1_given=true;
      ++matched;
    }

        // Check if current parameter is not valid
    if(matched==0)
    {
      cerr<<"Error: parameter "<<argv_stl[i]<<" not valid."<<endl;
      return ERROR;
    }
    ++i;
  }
  return OK;
}

//--------------- checkParameters function
int checkParameters(solve_lp_probl_pars& pars)
{
  if(!pars.f_given)
  {
    cerr<<"Error: -f parameter not given!"<<endl;
    return ERROR;
  }

  if(!pars.o_given)
  {
    cerr<<"Error: -o parameter not given!"<<endl;
    return ERROR;
  }

  return OK;
}

//--------------- printParameters function
void printParameters(solve_lp_probl_pars pars)
{
  cerr<<"File with lp problem in CPLEX format: "<<pars.f_str<<endl;
  cerr<<"Output file name: "<<pars.o_str<<endl;
}

//--------------- printUsage function
void printUsage(void)
{
  cerr<<"Usage: solve_lp_probl      -f <string> -o <string>\n";
  cerr<<"                           [-v|-v1] [--help] [--version]\n\n";
  cerr<<"-f <string>                File lp problem in CPLEX format.\n";
  cerr<<"-o <string>                Output file name.\n";
  cerr<<"-v | -v1                   Verbose modes.\n";
  cerr<<"--help                     Display this help and exit.\n";
  cerr<<"--version                  Output version information and exit.\n";
}

//--------------- version function
void version(void)
{
  cerr<<"solve_lp_probl is part of the fcap package "<<endl;
  cerr<<"fcap version "<<FCAP_VERSION<<endl;
//  cerr<<"fcap is GNU software written by Daniel Ortiz"<<endl;
}
