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
        genemap={}
        metabmap={}
        reactmap={}
        gprrlist=[]
        rlowbndlist=[]
        ruppbndlist=[]
        stoicheqdict={}
        objfun=[]
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
    sbmli.reactmap=read_react_map(sbmlf+"_reaction_ids.csv")
    sbmli.gprrlist=read_gprr_list(sbmlf+"_gpr_rules.csv")
    sbmli.rlowbndlist=read_rlowbnd_list(sbmlf+"_reaction_lowbnds.csv")
    sbmli.ruppbndlist=read_ruppbnd_list(sbmlf+"_reaction_uppbnds.csv")
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

    # Iterate over gpr rules
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
def gen_yplus_h_name(i):
    return "yph%05d" % (i)

##################################################
def gen_yplus_l_name(i):
    return "ypl%05d" % (i)

##################################################
def gen_yminus_name(i):
    return "ymh%05d" % (i)

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
