# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_help():
    print >> sys.stderr, "plot_metab_network -s <string> -d <string> -f <string> [-e <string>] [-t <int>]"
    print >> sys.stderr, "                   [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-s <string> :    prefix of SBML info files"
    print >> sys.stderr, "-d <string> :    file with data for reaction id's (fluxes, p-values, etc.)"
    print >> sys.stderr, "-f <string> :    file containing the reaction id's to be included in the plot"
    print >> sys.stderr, "-e <string> :    file containing the metabolite id's considered to be external"
    print >> sys.stderr, "                 (if not provided, all metabolites are considered internal)"
    print >> sys.stderr, "-t <int>    :    plot type"
    print >> sys.stderr, "                 0 -> reactions + reaction senses + metabolites +"
    print >> sys.stderr, "                      stoichiometric coefs. (default option)"
    print >> sys.stderr, "                 1 -> reactions + reaction senses + metabolites"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
def read_reactdata(dataf):
    
    # initialize result
    reactdata={}

    # read file line by line
    file = open(dataf, 'r')
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        reactid=int(fields[0])
        value=float(fields[1])
        reactdata[reactid]=value
        
    # return result
    return reactdata

##################################################
def load_rids_filt_file(filterf):

    # initialize result
    rids={}

    # read file line by line
    file = open(filterf, 'r')
    for line in file:
        line=line.strip("\n")
        fields=line.split(" ")
        idx=int(fields[0])
        rids[idx]=True

    # return result
    return rids

##################################################
def load_extern_metab_file(extermf):

    # initialize result
    extern_metab_set=set()

    # read file line by line
    file = open(filterf, 'r')
    for line in file:
        line=line.strip("\n")
        fields=line.split(" ")
        extern_metab_set.insert(int(fields[0]))

    # return result
    return extern_metab_set

##################################################
def assign_color(flux):

    # Set value of result variable
    if(flux>-1 and flux<1):
        result="black"
    elif(flux>=1):
        result="red"
    elif(flux<=-1):
        result="blue"

    # Return result
    return result


##################################################
def print_header():

    print "digraph word_graph {"
    print "rankdir=LR;"
    print "overlap=false;"
    print "splines=true;"
#    print "splines=ortho;"
    print "K=1;"

##################################################
def box_reaction_representation(sbmli,included_rids):

    print "node [shape = box]; ",
    for rid in sbmli.rlowbndmap:
        if(rid in included_rids):
            vname=fba.gen_vname(rid)
            reactname=sbmli.reactmap[rid]
            print "_"+fba.clean_string(reactname),";",
    print ""

##################################################
def circle_metab_representation():

    print "node [shape = circle];"

##################################################
def gen_node_id_and_name_for_metab(extern_metab_set,mid,metabname,reactname):

    if(mid in extern_metab_set):
        return metabname+"_"+reactname,metabname
    else:
        return metabname,metabname

##################################################
def print_arc_zero(sbmli,extern_metab_set,reactdata,vcoef,mid,rid):
    # Initialize variables
    metabname=sbmli.metabmap[mid]
    color=assign_color(reactdata[rid])
    reactname=sbmli.reactmap[rid]
    clreactname=fba.clean_string(reactname)
    clmetabname=fba.clean_string(metabname)
    nodeid,nodename=gen_node_id_and_name_for_metab(extern_metab_set,mid,clmetabname,clreactname)
    metabnode_string="{"+"_"+nodeid +" [label=\""+nodename+"\"]}"
    reactnode_string="{"+"_"+clreactname +" [label=\""+reactname+"\"]}"

    # Print arc
    if(vcoef>=0):
        print reactnode_string,"->", metabnode_string,"[ label = \"",vcoef,"\", color =",color,"];"
    else:
        print metabnode_string,"->",reactnode_string,"[ label = \"",vcoef,"\", color =",color,"];"

##################################################
def print_arc_one(sbmli,extern_metab_set,reactdata,vcoef,mid,rid):
    # Initialize variables
    metabname=sbmli.metabmap[mid]
    color=assign_color(reactdata[rid])
    reactname=sbmli.reactmap[rid]
    clreactname=fba.clean_string(reactname)
    clmetabname=fba.clean_string(metabname)
    nodeid,nodename=gen_node_id_and_name_for_metab(extern_metab_set,mid,clmetabname,clreactname)
    metabnode_string="{"+"_"+nodeid +" [label=\""+nodename+"\"]}"
    reactnode_string="{"+"_"+clreactname +" [label=\""+reactname+"\"]}"

    # Print arc
    if(vcoef>=0):
        print reactnode_string,"->",metabnode_string,"[ color =",color,"];"
    else:
        print metabnode_string,"->",reactnode_string,"[ color =",color,"];"

##################################################
def process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation):
    # Iterate over metabolites
    for mid in sbmli.metabmap:
        for i in range(len(sbmli.stoicheqdict[mid])):
            rid=sbmli.stoicheqdict[mid][i].v
            if(rid in included_rids):
                vcoef=sbmli.stoicheqdict[mid][i].coef
                if(arc_representation==0):
                    print_arc_zero(sbmli,extern_metab_set,reactdata,vcoef,mid,rid)
                elif(arc_representation==1):
                    print_arc_one(sbmli,extern_metab_set,reactdata,vcoef,mid,rid)

##################################################
def print_footer():
    print "}"

##################################################
def print_metab_network_type_zero(sbmli,extern_metab_set,reactdata,included_rids):

    # Print header
    print_header()

    ## Set representation for the different nodes

    # Set representation for reactions
    box_reaction_representation(sbmli,included_rids)

    # Set representation for metabolites
    circle_metab_representation()

    ## Process stochiometric relations
    arc_representation=0
    process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation)

    # Print footer
    print_footer()

##################################################
def print_metab_network_type_one(sbmli,extern_metab_set,reactdata,included_rids):

    # Print header
    print_header()

    ## Set representation for the different nodes

    # Set representation for reactions
    box_reaction_representation(sbmli,included_rids)

    # Set representation for metabolites
    circle_metab_representation()

    ## Process stochiometric relations
    arc_representation=1
    process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation)

    # Print footer
    print_footer()
    
##################################################
def plot_network(sbmlf,dataf,filterf,extermf,pltype):
    # load sbml info
    sbmli=fba.extract_sbml_info(sbmlf)

    # read reaction data from sol file
    reactdata=read_reactdata(dataf)

    # load reacton id's to include in the plot if a file was given
    if(filterf!=""):
        included_rids=load_rids_filt_file(filterf)

    # load
    extern_metab_set=set()
    if(extermf!=""):
        extern_metab_set=load_extern_metab_file(extermf)

    # print metabolic network
    if(pltype==0):
        print_metab_network_type_zero(sbmli,extern_metab_set,reactdata,included_rids)
    elif(pltype==1):
        print_metab_network_type_one(sbmli,extern_metab_set,reactdata,included_rids)

##################################################
def main(argv):
    # take parameters
    s_given=False
    d_given=False
    f_given=False
    filterf=""
    e_given=False
    extermf=""
    t_given=False
    pltype=0
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:d:f:e:t:",["help","sbmlf=","dataf=","filterf=","extermf","pltype="])
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
            elif opt in ("-d", "--dataf"):
                dataf = arg
                d_given=True
            elif opt in ("-f", "--filterf"):
                filterf = arg
                f_given=True
            elif opt in ("-e", "--extermf"):
                extermf = arg
                e_given=True
            elif opt in ("-t", "--pltype"):
                pltype = int(arg)
                t_given=True

    # print parameters
    if(s_given==True):
        print >> sys.stderr, "s is %s" % (sbmlf)
    else:
        print >> sys.stderr, "Error: -s option not given"
        sys.exit(2)

    if(d_given==True):
        print >> sys.stderr, "d is %s" % (dataf)
    else:
        print >> sys.stderr, "Error: -d option not given"
        sys.exit(2)

    if(f_given==True):
        print >> sys.stderr, "f is %s" % (filterf)
    else:
        print >> sys.stderr, "Error: -f option not given"
        sys.exit(2)

    # create lp file according to selected criterion
    plot_network(sbmlf,dataf,filterf,extermf,pltype)
        
if __name__ == "__main__":
    main(sys.argv)
