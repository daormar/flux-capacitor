"""
Flux Capacitor package
Copyright (C) 2015-2018 Daniel Ortiz-Mart\'inez
 
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
"""
 
# *- python -*

import sys, getopt, numpy
from libsbml import *

##################################################
def extract_info(model_file,outprefix):
    # Read SBML file
    reader = SBMLReader()
    document = reader.readSBMLFromFile(model_file)
#    numerr=document.getNumErrors()
    model = document.getModel()

    ### Write files

    # Write gene ids
    gene_ids_file=outprefix+"_gene_ids.csv"
    write_gene_ids(model,gene_ids_file)
    
##################################################
def write_gene_ids(model,gene_ids_file):
    file=open(react_file_name, 'w')
    
##################################################
def print_help():
    print("extract_sbml_mod_info -m <string> -o <string> [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-m <string> :    path to file with metabolic model in SBML format", file=sys.stderr)
    print("-o <string> :    prefix for output files", file=sys.stderr)
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

##################################################
def main(argv):
    # take parameters
    m_given=False
    m_val = ""
    o_given=False
    o_val=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hm:o:",["help","model=","outpref="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                print_help()
                sys.exit()
            elif opt in ("-m", "--model"):
                m_val = arg
                m_given=True
            elif opt in ("-o", "--outpref"):
                o_val = arg
                o_given=True

    # print parameters
    if(m_given==True):
        print("m is %s" % (m_val), file=sys.stderr)
    else:
        print("Error: -m option not given", file=sys.stderr)
        sys.exit(2)

    if(o_given==True):
        print("o is %s" % (o_val), file=sys.stderr)
    else:
        print("Error: -o option not given", file=sys.stderr)
        sys.exit(2)

    # Process parameters
    extract_info(m_val,o_val)
        
if __name__ == "__main__":
    main(sys.argv)
