# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys

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
        self.genemap={}
        self.metabmap={}
        self.reactmap={}
        self.gprrmap={}
        self.rlowbndmap={}
        self.ruppbndmap={}
        self.stoicheqdict={}
        self.objfun=[]

##################################################
class sparse_st_mat_elem:
    def __init__(self):
        self.v=0
        self.coef=0

##################################################
def extract_sbml_info(sbmlf):
    sbmli=sbml_info()
    sbmli.genemap=read_gene_map(sbmlf+"_gene_ids.csv")
    sbmli.metabmap=read_metab_map(sbmlf+"_metabolite_ids.csv")
    sbmli.reactmap=read_react_map(sbmlf+"_reaction_ids.csv")
    sbmli.gprrmap=read_gprr_map(sbmlf+"_gpr_rules.csv")
    sbmli.rlowbndmap=read_rlowbnd_map(sbmlf+"_reaction_lowbnds.csv")
    sbmli.ruppbndmap=read_ruppbnd_map(sbmlf+"_reaction_uppbnds.csv")
    sbmli.stoicheqdict=read_sparse_st_matrix(sbmlf+"_sparse_st_matrix.csv")
    sbmli.objfun=read_objfun(sbmlf+"_obj_fun.csv")
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
def read_react_map(filename):
    reactmap={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        reactmap[int(fields[0])]=fields[1]
    return reactmap

##################################################
def read_gprr_map(filename):
    gprrmap={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        gprrmap[int(fields[0])]=fields[1]
    return gprrmap

##################################################
def read_rlowbnd_map(filename):
    rlowbndmap={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        rlowbndmap[int(fields[0])]=float(fields[1])
    return rlowbndmap

##################################################
def read_ruppbnd_map(filename):
    ruppbndmap={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        ruppbndmap[int(fields[0])]=float(fields[1])
    return ruppbndmap

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
def read_objfun(filename):
    objfun=[]
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        objfun.append(int(fields[0]))
    return objfun

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
    return idmap_info

##################################################
def obtain_hlreact_set(sbmli,abspres_info):
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

    # Iterate over gpr rules
    for key,val in sbmli.gprrmap.items():
        if(val!=""):
            try:
                tmp=eval(val)
                if(tmp==NA_()):
                    tmp=0.5
            except TypeError:
                print >> sys.stderr,"obtain_hlreact_set(), TypeError:",val
                tmp=0.5
            except SyntaxError:
                print >> sys.stderr,"obtain_hlreact_set(), SyntaxError:",val
                tmp=0.5
        else:
            tmp=0.5
#        print "****",val,tmp
        result.append(tmp)

    # Return result
    return result

##################################################
def gen_vname(i):
    return "v%06d" % (i)

##################################################
def gen_yplus_h_name(i):
    return "yph%06d" % (i)

##################################################
def gen_yplus_l_name(i):
    return "ypl%06d" % (i)

##################################################
def gen_yminus_name(i):
    return "ymh%06d" % (i)

##################################################
def clean_string(s):

    # Clean string to ensure it is appropriate to be included in dot
    # files
    result=s.replace("[","_")
    result=result.replace("]","_")
    result=result.replace("(","_")
    result=result.replace(")","_")

    # Return result
    return result
