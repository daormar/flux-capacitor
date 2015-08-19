# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt

##################################################
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

##################################################
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

##################################################
class sparse_st_mat_elem:
    def __init__(self):
        v=0
        coef=0

##################################################
def extract_sbml_info(sbmlf):
    sbmli=sbml_info()
    sbmli.genemap=read_gene_map(sbmlf+"_gene_ids.csv")
    sbmli.metabmap=read_metab_map(sbmlf+"_metabolite_ids.csv")
    sbmli.gprrlist=read_gprr_list(sbmlf+"_gpr_rules.csv")
    sbmli.rlowbndlist=read_rlowbnd_list(sbmlf+"_reaction_lowbnds.csv")
    sbmli.ruppbndlist=read_ruppbnd_list(sbmlf+"_reaction_uppbnds.csv")
    sbmli.stoicheqdict=read_sparse_st_matrix(sbmlf+"_sparse_st_matrix.csv")
    return sbmli

##################################################
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

##################################################
def read_metab_map(filename):
    metabmap={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        metabmap[int(fields[0])]=fields[1]
    return metabmap

##################################################
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

##################################################
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

##################################################
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

##################################################
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
            key=int(fields[0])
            if key in stoicheqdict:
                elem=sparse_st_mat_elem()
                elem.v=int(fields[1])
                elem.coef=float(fields[2])
                stoicheqdict[key].append(elem) 
            else:
                elem=sparse_st_mat_elem()
                elem.v=int(fields[1])
                elem.coef=float(fields[2])
                stoicheqdict[key]=[]
                stoicheqdict[key].append(elem)
        lineno=lineno+1

    # return result
    return stoicheqdict

##################################################
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

##################################################
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

##################################################
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
    result.append(None)

    # Iterate over gpr reactions
    for i in xrange(1,len(sbmli.gprrlist)):
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

##################################################
def gen_vname(i):
    return "v%05d" % (i)

##################################################
def gen_yplus_name(i):
    return "yp%05d" % (i)

##################################################
def gen_yminus_name(i):
    return "ym%05d" % (i)

##################################################
def print_obj_func(hlreact_set):

    # Print header
    print "Maximize"

    # Print objective function
    for i in range(1,len(hlreact_set)):
        if(i<len(hlreact_set)-1):
            if(hlreact_set[i]==1):
                st="+ " + gen_yplus_name(i) + " + " + gen_yminus_name(i)
                print st,
            elif(hlreact_set[i]==0):
                st="+ " + gen_yplus_name(i)
                print st,
        else:
            if(hlreact_set[i]==1):
                st="+ " + gen_yplus_name(i) + " + " + gen_yminus_name(i)
                print st
            elif(hlreact_set[i]==0):
                st="+ " + gen_yplus_name(i)
                print st
            elif(hlreact_set[i]==0.5):
                print ""


    # Print footer
    print ""

##################################################
def print_flux_boundaries(sbmli,hlreact_set):

    # Print header
    print "Bounds"

    # Print flux upper and lower bounds
    for i in range(1,len(sbmli.rlowbndlist)):
        varname=gen_vname(i)
        print sbmli.rlowbndlist[i],"<=",varname,"<=",sbmli.ruppbndlist[i]

    # Print footer
    print ""

##################################################
def print_bin_vars(hlreact_set):

    # Print header
    print "Binary"

    # Print binary variables
    for i in range(1,len(hlreact_set)):
        if(i<len(hlreact_set)-1):
            if(hlreact_set[i]==1):
                st=gen_yplus_name(i) + " " + gen_yminus_name(i)
                print st,
            elif(hlreact_set[i]==0):
                st=gen_yplus_name(i)
                print st,
        else:
            if(hlreact_set[i]==1):
                st=gen_yplus_name(i) + " " + gen_yminus_name(i)
                print st
            elif(hlreact_set[i]==0):
                st=gen_yplus_name(i)
                print st
            elif(hlreact_set[i]==0.5):
                print ""

    # Print footer
    print ""

##################################################
def print_constraints(sbmli,hlreact_set):
    
    # Print header
    print "Subject To"

    # Iterate over metabolites
    for k in sbmli.metabmap:
        # Obtain metabname and modify it to avoid problems with solvers
        # such as CPLEX
        metabname=sbmli.metabmap[k]
        metabname=metabname.replace("[","_")
        metabname=metabname.replace("]","_")
        # Print constraint
        print "_"+metabname+":",
        for i in range(len(sbmli.stoicheqdict[k])):
            vname=gen_vname(sbmli.stoicheqdict[k][i].v)
            if(sbmli.stoicheqdict[k][i].coef > 0.0):
                print "+",sbmli.stoicheqdict[k][i].coef,vname,
            else:
                print "-",-sbmli.stoicheqdict[k][i].coef,vname,
        print "= 0"

    # Init epsilon
    epsilon=1

    # Print lower bounds for R_H
    for i in range(1,len(sbmli.rlowbndlist)):
        if(hlreact_set[i]==1):
            vname=gen_vname(i)
            ypname=gen_yplus_name(i)
            coef=sbmli.rlowbndlist[i]-epsilon
            if(coef>=0):
                print vname,"+",coef,ypname,">=",sbmli.rlowbndlist[i]
            else:
                print vname,coef,ypname,">=",sbmli.rlowbndlist[i]

    # Print upper bounds for R_H
    for i in range(1,len(sbmli.ruppbndlist)):
        if(hlreact_set[i]==1):
            vname=gen_vname(i)
            ymname=gen_yminus_name(i)
            coef=sbmli.ruppbndlist[i]+epsilon
            if(coef>=0):
                print vname,"+",coef,ymname,"<=",sbmli.ruppbndlist[i]
            else:
                print vname,coef,ymname,"<=",sbmli.ruppbndlist[i]

    # Print upper bounds for R_L
    for i in range(1,len(sbmli.ruppbndlist)):
        if(hlreact_set[i]==0):
            vname=gen_vname(i)
            ypname=gen_yplus_name(i)
            # print vname,">=",sbmli.rlowbndlist[i],"-",ypname
            # print vname,"<=",sbmli.ruppbndlist[i],"-",ypname
            print vname,"+",ypname,">=",sbmli.rlowbndlist[i]
            print vname,"+",ypname,"<=",sbmli.ruppbndlist[i]

    # Print footer
    print "End"

##################################################
def print_lp_problem(sbmli,hlreact_set):
    
    # Print objective function
    print_obj_func(hlreact_set)

    # Print constraints
    print_constraints(sbmli,hlreact_set)

    # Print flux boundaries
    print_flux_boundaries(sbmli,hlreact_set)

    # Print ids of binary variables
    print_bin_vars(hlreact_set)

##################################################
def print_help():
    print >> sys.stderr, "create_lp_file -s <string> -a <string> -m <string> -c <int> [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-s <string> :    prefix of SBML info files"
    print >> sys.stderr, "-a <string> :    file with absent/present genes data"
    print >> sys.stderr, "-m <string> :    file with mapping between probeset ids and entrez ids" 
    print >> sys.stderr, "-c <int>    :    fba criterion used to generate the lp file. The criterion"
    print >> sys.stderr, "                 can be selected from the following list,"    
    print >> sys.stderr, "                 0 -> Shlomi et al. 2008"    
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
def create_lp_file_shlomi(sbmlf,abspresf,idmapf):
    # load sbml info
    sbmli=extract_sbml_info(sbmlf)

    # load absent/present genes info
    abspres_info=load_abspres_info(abspresf)

    # load mapping between probeset ids and entrez ids
    idmap_info=load_idmap_info(idmapf)

    # Obtain highly/lowly expressed reactions
    hlreact_set=obtain_hlreact_set(sbmli,abspres_info,idmap_info)

    # print problem in lp format
    print_lp_problem(sbmli,hlreact_set)

##################################################
def main(argv):
    # take parameters
    s_given=False
    a_given=False
    sbmlf = ""
    abspresf= ""
    c_given=False
    crit=0
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:a:m:c:",["sbmlf=","abspresf=","idmapf=","crit="])
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
            elif opt in ("-s", "--sbmlf"):
                sbmlf = arg
                s_given=True
            elif opt in ("-a", "--abspresf"):
                abspresf = arg
                a_given=True
            elif opt in ("-m", "--idmapf"):
                idmapf = arg
                m_given=True
            elif opt in ("-c", "--crit"):
                crit = int(arg)
                c_given=True

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
        print >> sys.stderr, "m is %s" % (idmapf)
    else:
        print >> sys.stderr, "Error: -m option not given"
        sys.exit(2)

    print >> sys.stderr, "c is %s" % (crit)

    # create lp file according to selected criterion
    if(crit==0):
        create_lp_file_shlomi(sbmlf,abspresf,idmapf)
        
if __name__ == "__main__":
    main(sys.argv)
