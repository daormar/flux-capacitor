# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, numpy

##################################################
class params_info:
    def __init__(self):
        self.phenotypefile=''

##################################################
def load_params_info(filename):
    pinfo=params_info()
    file = open(filename, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(" ")
        if(fields[0]=="-p"):
            pinfo.phenotypefile=fields[3]
        
    return pinfo

##################################################
class phenotype_info:
    def __init__(self):
        self.sample_to_pheno={}

##################################################
def load_phenotype_info(filename):
    phenoinfo=phenotype_info()
    file = open(filename, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        if(lineno!=1):
            phenoinfo.sample_to_pheno[fields[1]]=fields[2]
            
        lineno=lineno+1

    return phenoinfo

##################################################
def load_reaction_ids(react_file):
    reacts=set()
    file = open(react_file, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        reacts.add(fields[1])

    return reacts
    
##################################################
def load_fluxes_for_samples(fluxes_files_dir,phenoinfo):
    sample_fluxes={}
    for sample in phenoinfo.sample_to_pheno.keys():
        fluxes_file=fluxes_files_dir+"/fluxes_"+sample+".csv"
        sample_fluxes[sample]=load_fluxes_for_sample(fluxes_file)
    return sample_fluxes

##################################################
def load_fluxes_for_sample(fluxes_file):
    sflux={}
    file = open(fluxes_file, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        sflux[fields[1]]=fields[2]
        
    return sflux
    
##################################################
def print_gathered_fluxes(reacts,sample_fluxes):

    # Print first row
    samples_csv=",".join(sample_fluxes.keys())
    print ","+samples_csv

    # Iterate over reactions
    for react in reacts:
        # Print sample fluxes per reaction
        fluxes=""
        for sample in sample_fluxes.keys():
            fluxes=fluxes+","+sample_fluxes[sample][react]
        print react+fluxes

##################################################
def print_help():
    print >> sys.stderr, "gather_sample_fluxes -d <string> [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-d <string> :    output directory of auto_fba experiment"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
def main(argv):
    # take parameters
    d_given=False
    d_val = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hd:",["help","fbadir="])
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
            elif opt in ("-d", "--fbadir"):
                d_val = arg
                d_given=True

    # print parameters
    if(d_given==True):
        print >> sys.stderr, "d is %s" % (d_val)
    else:
        print >> sys.stderr, "Error: -d option not given"
        sys.exit(2)

    # load parameter file
    parfile=d_val+"/params.txt"
    parinfo=load_params_info(parfile)

    # load phenotype file
    phenoinfo=load_phenotype_info(parinfo.phenotypefile)

    # load reaction ids
    react_file=d_val+"/minfo/model_reaction_ids.csv"
    reacts=load_reaction_ids(react_file)
    
    # load fluxes for each sample
    fluxes_files_dir=d_val+"/sol"
    sample_fluxes=load_fluxes_for_samples(fluxes_files_dir,phenoinfo)

    # print gathered fluxes
    print_gathered_fluxes(reacts,sample_fluxes)
    
if __name__ == "__main__":
    main(sys.argv)
