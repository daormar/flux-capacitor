# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt

# define sbml_info class
class sbml_info:
    def __init__(self):
        genemap={}
        metabmap={}
        reactlist=[]
        rlowbndlist=[]
        ruppbndlist=[]
        stoicheqdict={}
    # def __repr__(self):
    #     return str(self.seqid)+" "+str(self.genomepos)+" "+str(self.lineno)+" "+str(self.charno)

# define sparse_st_mat_elem class
class sparse_st_mat_elem:
    def __init__(self):
        v=0
        coef=0

# extract_sbml_info
def extract_sbml_info(sbmlf):
    sbmli=sbml_info()
    sbmli.genemap=read_gene_map(sbmlf+"_gene_ids.csv")
    sbmli.metabmap=read_metab_map(sbmlf+"_metabolite_ids.csv")
    sbmli.reactlist=read_react_list(sbmlf+"_reaction_ids.csv")
    sbmli.rlowbndlist=read_rlowbnd_list(sbmlf+"_reaction_lowbnds.csv")
    sbmli.ruppbndlist=read_ruppbnd_list(sbmlf+"_reaction_uppbnds.csv")
    sbmli.stoicheqdict=read_sparse_st_matrix(sbmlf+"_sparse_st_matrix.csv")
    return sbmli

# read_gene_map
def read_gene_map(filename):
    genemap={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        genemap[int(fields[0])]=fields[1]
    return genemap

# read_metab_map
def read_metab_map(filename):
    metabmap={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        metabmap[int(fields[0])]=fields[1]
    return metabmap

# read_react_list
def read_react_list(filename):
    reactlist=[]
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        reactlist.append(fields[0])
    return reactlist

# read_rlowbnd_list
def read_rlowbnd_list(filename):
    rlowbndlist=[]
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        rlowbndlist.append(float(fields[0]))
    return rlowbndlist

# read_ruppbnd_list
def read_ruppbnd_list(filename):
    ruppbndlist=[]
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        ruppbndlist.append(float(fields[0]))
    return ruppbndlist

# read_sparse_st_matrix
def read_sparse_st_matrix(filename):
    # initialize variables
    stoicheqdict={}

    # open file
    file = open(filename, 'r')
    
    # read file line by line
    lineno=0
    for line in file:
        if(lineno>1):
            line=line.strip("\n")
            fields=line.split(" ")
#            print fields
            if fields[0] in stoicheqdict:
                elem=sparse_st_mat_elem()
                elem.v=fields[1]
                elem.coef=fields[2]
                stoicheqdict[fields[0]].append(elem) 
            else:
                elem=sparse_st_mat_elem()
                elem.v=fields[1]
                elem.coef=fields[2]
                stoicheqdict[fields[0]]=[]
                stoicheqdict[fields[0]].append(elem)
        lineno=lineno+1

    # return result
    return stoicheqdict

# load_abspres_info
def load_abspres_info(abspresf):
    abspres_info={}
    file = open(abspresf, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        if(fields[0]!="NA"):
            abspres_info[fields[0]]=fields[1]
#            print fields[0],fields[1]
    return abspres_info

# print_cplex_problem
def print_cplex_problem(sbmli,abspres_info):
    # Print objective function
    print "TBD"

# main
def main(argv):
    # take parameters
    s_given=False
    a_given=False
    sbmlf = ""
    abspresf= ""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:a:",["sbmlf=","abspresf="])
    except getopt.GetoptError:
        print >> sys.stderr, "create_cplex_file -s <string> -a <string>"
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print "create_cplex_file -s <string> -a <string>"
            print ""
            print "-s <string> :    prefix of SBML info files"
            print "-a <string> :    file with absent/present genes data"
            sys.exit()
        elif opt in ("-s", "--sbmlf"):
            sbmlf = arg
            s_given=True
        elif opt in ("-a", "--abspresf"):
            abspresf = arg
            a_given=True

    # print parameters
    if(s_given==True):
        print >> sys.stderr, "s is %s" % (sbmlf)
    else:
        print >> sys.stderr, "Error: -s option not given"
        sys.exit(2)

    if(a_given==True):
        print >> sys.stderr, "a is %s" % (abspresf)
    else:
        print >> sys.stderr, "Error: -a option not given"
        sys.exit(2)

    # load sbml info
    sbmli=extract_sbml_info(sbmlf)

    # load absent/present genes info
    abspres_info=load_abspres_info(abspresf)

    # print problem in cplex format
    print_cplex_problem(sbmli,abspres_info)

if __name__ == "__main__":
    main(sys.argv)
