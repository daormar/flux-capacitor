# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_help():
    print >> sys.stderr, "plot_fba_network -s <string> -p <string> [-f <string>] [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-s <string> :    prefix of SBML info files"
    print >> sys.stderr, "-p <string> :    file with lp solution generated with CPLEX"
    print >> sys.stderr, "-f <string> :    file containing the reaction id's to be included in the plot" 
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
def read_fluxes(solf):
    
    # initialize result
    fluxes={}

    # read file line by line
    file = open(solf, 'r')
    for line in file:
        line=line.strip("\n")
        fields=line.split("\"")
        if(len(fields)>5 and fields[1].startswith('v')):
            varname=fields[1]
            value=float(fields[5])
            fluxes[varname]=value
#            print varname,value
        
    # return result
    return fluxes

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
        result="green"

    # Return result
    return result

##################################################
def print_fba_network(sbmli,fluxes,included_rids):

    # Print header
    print "digraph word_graph {"
    print "rankdir=LR;"
    print "size=\"8,5\""

    ## Set representation for the different nodes

    # Set representation for reactions
    print "node [shape = box]; ",

    for vid in range(1,len(sbmli.rlowbndlist)):
        if(vid in included_rids):
            vname=fba.gen_vname(vid)
            reactname=sbmli.reactmap[vid]
            print "_"+reactname,";",
    print ""

    # Set representation for metabolites
    print "node [shape = circle];"

    ## Process stochiometric relations

    # Iterate over metabolites
    for k in sbmli.metabmap:
        metabname=sbmli.metabmap[k]
        for i in range(len(sbmli.stoicheqdict[k])):
            vid=sbmli.stoicheqdict[k][i].v
            if(vid in included_rids):
                vcoef=sbmli.stoicheqdict[k][i].coef
                vname=fba.gen_vname(vid)
                color=assign_color(fluxes[vname])
                reactname=sbmli.reactmap[vid]
                if(vcoef>=0):
                    print "_"+reactname,"->","_"+fba.clean_string(metabname),"[ label = \"",vcoef,"\", color =",color," ];"
                else:
                    print "_"+fba.clean_string(metabname),"->","_"+reactname,"[ label = \"",vcoef,"\", color =",color," ];"

    # Print footer
    print "}"

##################################################
def plot_network(sbmlf,solf,filterf):
    # load sbml info
    sbmli=fba.extract_sbml_info(sbmlf)

    # read fluxes from sol file
    fluxes=read_fluxes(solf)

    # load reacton id's to include in the plot if a file was given
    if(filterf!=""):
        included_rids=load_rids_filt_file(filterf)

    # print fba network
    print_fba_network(sbmli,fluxes,included_rids)

##################################################
def main(argv):
    # take parameters
    s_given=False
    p_given=False
    f_given=False
    filterf=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:p:f:",["sbmlf=","solf=","filterf="])
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
            elif opt in ("-p", "--solf"):
                solf = arg
                p_given=True
            elif opt in ("-f", "--filterf"):
                filterf = arg
                f_given=True

    # print parameters
    if(s_given==True):
        print >> sys.stderr, "s is %s" % (sbmlf)
    else:
        print >> sys.stderr, "Error: -s option not given"
        sys.exit(2)

    if(p_given==True):
        print >> sys.stderr, "p is %s" % (solf)
    else:
        print >> sys.stderr, "Error: -p option not given"
        sys.exit(2)

    if(f_given==True):
        print >> sys.stderr, "f is %s" % (filterf)

    # create lp file according to selected criterion
    plot_network(sbmlf,solf,filterf)
        
if __name__ == "__main__":
    main(sys.argv)
