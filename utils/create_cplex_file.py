# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt

# define NA_ class
class NA_(object):
    instance = None # Singleton (so `val is NA` will work)
    def __new__(self):
        if NA_.instance is None:
            NA_.instance = super(NA_, self).__new__(self)
        return NA_.instance
    def __str__(self): return "NA"
    def __repr__(self): return "NA_()"
    def __and__(self, other):
        if self is other or other:
            return self
        else:
            return other
    __rand__ = __and__
    def __or__(self, other):
        if self is other or other:
            return other
        else:
            return self
    __ror__ = __or__
    def __xor__(self, other):
        return self
    __rxor__ = __xor__
    def __eq__(self, other):
        return self is other
    __req__ = __eq__
    def __nonzero__(self):
        raise TypeError("bool(NA) is undefined.")

# define sbml_info class
class sbml_info:
    def __init__(self):
        genemap={}
        metabmap={}
        reactlist=[]
        gprrlist=[]
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
#    sbmli.reactlist=read_react_list(sbmlf+"_reaction_ids.csv")
    sbmli.gprrlist=read_gprr_list(sbmlf+"_gpr_rules.csv")
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
        idx=fields[1].rfind(".")
#        print fields[1],idx,fields[1][0:idx]
#        genemap[int(fields[0])]=fields[1]
        genemap[int(fields[0])]=fields[1][0:idx]
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

# # read_react_list
# def read_react_list(filename):
#     reactlist=[]
#     file = open(filename, 'r')
#     # read file line by line
#     for line in file:
#         line=line.strip("\n")
#         fields=line.split(",")
#         reactlist.append(fields[0])
#     return reactlist

def read_gprr_list(filename):
    gprrlist=[]
    gprrlist.append(None)
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        gprrlist.append(fields[0])
    return gprrlist

# read_rlowbnd_list
def read_rlowbnd_list(filename):
    rlowbndlist=[]
    rlowbndlist.append(None)
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
    ruppbndlist.append(None)
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

# load_idmap_info
def load_idmap_info(idmapf):
    idmap_info={}
    file = open(idmapf, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        if(fields[0]!="NA"):
            idmap_info[fields[0]]=fields[1]
#            print fields[0],fields[1]
    return idmap_info

# obtain_hlreact_set
def obtain_hlreact_set(sbmli,abspres_info,idmap_info):
    # Create x vector
    x=[]
    x.append(-1)
    for i in sbmli.genemap:
        if(sbmli.genemap[i] in abspres_info):
#            print i,sbmli.genemap[i],abspres_info[sbmli.genemap[i]]
            if(abspres_info[sbmli.genemap[i]]=="NA"):
                x.append(NA_())
            else:
                x.append(int(abspres_info[sbmli.genemap[i]]))
        else:
            x.append(NA_())
#            print i,sbmli.genemap[i],x[i]

    # Initialize result variable
    result=[]

    # Iterate over gpr reactions
    for i in xrange(len(sbmli.gprrlist)):
        if(sbmli.gprrlist[i]!=""):
            try:
                tmp=eval(sbmli.gprrlist[i])
                if(tmp==NA_()):
                    tmp=0.5
            except TypeError:
                print >> sys.stderr,"obtain_hlreact_set(), TypeError:",sbmli.gprrlist[i]
                tmp=0.5
            except SyntaxError:
                print >> sys.stderr,"obtain_hlreact_set(), SyntaxError:",sbmli.gprrlist[i]
                tmp=0.5
        else:
            tmp=0.5
#        print "****",sbmli.gprrlist[i],tmp
        result.append(tmp)

    # Return result
    return result

# print_obj_func
def print_obj_func(hlreact_set):
    print "Maximize"
    for i in range(1,len(hlreact_set)):
        if(i<len(hlreact_set)-1):
            if(hlreact_set[i]==1):
                st="+ yp%d + ym%d" % (i,i)
                print st,
            elif(hlreact_set[i]==0):
                st="+ yp%d" % (i)
                print st,
        else:
            if(hlreact_set[i]==1):
                st="+ yp%d + ym%d" % (i,i)
                print st
            elif(hlreact_set[i]==0):
                st="+ yp%d" % (i)
                print st
            elif(hlreact_set[i]==0.5):
                print ""

# print_flux_boundaries
def print_flux_boundaries(sbmli):
    # # Print lowerbounds
    # for i in range(1,len(sbmli.rlowbndlist)):
    #     varname="v%d" % (i)
    #     print varname,">",sbmli.rlowbndlist[i]

    # # Print upperbounds
    # for i in range(1,len(sbmli.ruppbndlist)):
    #     varname="v%d" % (i)
    #     print varname,"<",sbmli.ruppbndlist[i]    

    # Print flux upper and lower bounds
    for i in range(1,len(sbmli.rlowbndlist)):
        varname="v%d" % (i)
        print sbmli.rlowbndlist[i],"<",varname,"<",sbmli.ruppbndlist[i]

# print_bin_vars
def print_bin_vars(hlreact_set):
    print "Binary"
    for i in range(1,len(hlreact_set)):
        if(i<len(hlreact_set)-1):
            if(hlreact_set[i]==1):
                st="yp%d ym%d" % (i,i)
                print st,
            elif(hlreact_set[i]==0):
                st="yp%d" % (i)
                print st,
        else:
            if(hlreact_set[i]==1):
                st="yp%d ym%d" % (i,i)
                print st
            elif(hlreact_set[i]==0):
                st="yp%d" % (i)
                print st
            elif(hlreact_set[i]==0.5):
                print ""

# print_cplex_problem
def print_cplex_problem(sbmli,hlreact_set):
    
    # Print objective function
    print_obj_func(hlreact_set)

    ## Print constraints
    print "Subject To"

    # Print steady state constraints
    print "TBD"

    # Print flux boundaries
    print_flux_boundaries(sbmli)

    # Print ids of binary variables
    print_bin_vars(hlreact_set)

# main
def main(argv):
    # take parameters
    s_given=False
    a_given=False
    sbmlf = ""
    abspresf= ""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:a:m:",["sbmlf=","abspresf=","idmap="])
    except getopt.GetoptError:
        print >> sys.stderr, "create_cplex_file -s <string> -a <string> -m <string>"
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print "create_cplex_file -s <string> -a <string> -m <string>"
            print ""
            print "-s <string> :    prefix of SBML info files"
            print "-a <string> :    file with absent/present genes data"
            print "-m <string> :    file with mapping between probeset ids and entrez ids" 
            print ""
            sys.exit()
        elif opt in ("-s", "--sbmlf"):
            sbmlf = arg
            s_given=True
        elif opt in ("-a", "--abspresf"):
            abspresf = arg
            a_given=True
        elif opt in ("-m", "--idmap"):
            idmap = arg
            m_given=True

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

    if(m_given==True):
        print >> sys.stderr, "m is %s" % (idmap)
    else:
        print >> sys.stderr, "Error: -m option not given"
        sys.exit(2)

    # load sbml info
    sbmli=extract_sbml_info(sbmlf)

    # load absent/present genes info
    abspres_info=load_abspres_info(abspresf)

    # load mapping between probset ids and entrez ids
    idmap_info=load_idmap_info(idmap)

    # Obtain highly/lowly expressed reactions
    hlreact_set=obtain_hlreact_set(sbmli,abspres_info,idmap_info)

    # print problem in cplex format
    print_cplex_problem(sbmli,hlreact_set)

if __name__ == "__main__":
    main(sys.argv)
