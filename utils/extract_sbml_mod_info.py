# Author: Daniel Ortiz Mart\'inez
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
    print >> sys.stderr, "extract_sbml_mod_info -m <string> -o <string> [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-m <string> :    path to file with metabolic model in SBML format"
    print >> sys.stderr, "-o <string> :    prefix for output files"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

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
        print >> sys.stderr, "m is %s" % (m_val)
    else:
        print >> sys.stderr, "Error: -m option not given"
        sys.exit(2)

    if(o_given==True):
        print >> sys.stderr, "o is %s" % (o_val)
    else:
        print >> sys.stderr, "Error: -o option not given"
        sys.exit(2)

    # Process parameters
    extract_info(m_val,o_val)
        
if __name__ == "__main__":
    main(sys.argv)
