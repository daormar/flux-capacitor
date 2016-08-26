# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_help():
    print >> sys.stderr, "plot_fba_network -s <string> -d <string> -f <string> [-t <int>] [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-s <string> :    prefix of SBML info files"
    print >> sys.stderr, "-d <string> :    file with data for reaction id's"
    print >> sys.stderr, "-f <string> :    file containing the reaction id's to be included in the plot"
    print >> sys.stderr, "-t <int>    :    plot type"
    print >> sys.stderr, "                 0 -> reactions + reaction values + metabolites +"
    print >> sys.stderr, "                      stoichiometric coefs. (default option)"
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

##################################################
def set_reaction_representation(sbmli,included_rids):

    print "node [shape = box]; ",
    for rid in sbmli.rlowbndmap:
        if(rid in included_rids):
            vname=fba.gen_vname(rid)
            reactname=sbmli.reactmap[rid]
            print "_"+fba.clean_string(reactname),";",
    print ""

##################################################
def set_metab_representation():

    print "node [shape = circle];"

##################################################
def print_arc(metabname,vcoef,vname,color,reactname):
    if(vcoef>=0):
        print "_"+fba.clean_string(reactname),"->","_"+fba.clean_string(metabname),"[ label = \"",vcoef,"\", color =",color," ];"
    else:
        print "_"+fba.clean_string(metabname),"->","_"+fba.clean_string(reactname),"[ label = \"",vcoef,"\", color =",color," ];"

##################################################
def process_stoich_relations(sbmli,reactdata,included_rids):
    # Iterate over metabolites
    for k in sbmli.metabmap:
        metabname=sbmli.metabmap[k]
        for i in range(len(sbmli.stoicheqdict[k])):
            rid=sbmli.stoicheqdict[k][i].v
            if(rid in included_rids):
                vcoef=sbmli.stoicheqdict[k][i].coef
                vname=fba.gen_vname(rid)
                color=assign_color(reactdata[rid])
                reactname=sbmli.reactmap[rid]
                print_arc(metabname,vcoef,vname,color,reactname)

##################################################
def print_footer():
    print "}"

##################################################
def print_fba_network(sbmli,reactdata,included_rids):

    # Print header
    print_header()

    ## Set representation for the different nodes

    # Set representation for reactions
    set_reaction_representation(sbmli,included_rids)

    # Set representation for metabolites
    set_metab_representation()

    ## Process stochiometric relations
    process_stoich_relations(sbmli,reactdata,included_rids)

    # Print footer
    print_footer()
    
##################################################
def plot_network(sbmlf,dataf,filterf):
    # load sbml info
    sbmli=fba.extract_sbml_info(sbmlf)

    # read reaction data from sol file
    reactdata=read_reactdata(dataf)

    # load reacton id's to include in the plot if a file was given
    if(filterf!=""):
        included_rids=load_rids_filt_file(filterf)

    # print fba network
    print_fba_network(sbmli,reactdata,included_rids)

##################################################
def main(argv):
    # take parameters
    s_given=False
    d_given=False
    f_given=False
    filterf=""
    t_given=False
    pltype=0
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:d:f:t:",["help","sbmlf=","dataf=","filterf=","pltype="])
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
    if(pltype==0):
        plot_network(sbmlf,dataf,filterf)
        
if __name__ == "__main__":
    main(sys.argv)
