# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_help():
    print >> sys.stderr, "obtain_accessible_reacs -s <string> -r <string> [-e <string>]"
    print >> sys.stderr, "                        [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-s <string> :    prefix of SBML info files"
    print >> sys.stderr, "-r <string> :    file containing the reaction id's used as seed for the tool"
    print >> sys.stderr, "-e <string> :    file containing the metabolite id's considered to be external"
    print >> sys.stderr, "                 (if not provided, all metabolites are considered internal)"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""


##################################################
def read_arcs(sbmlf):
    
    # initalize result
    arc_set=set()

    # obtain file name containing arc information
    sparse_mat_file=sbmlf+"_sparse_st_matrix.csv"

    # read file line by line
    lineno=1
    file = open(sparse_mat_file, 'r')
    for line in file:
        if(lineno>2):
            line=line.strip("\n")
            fields=line.split(" ")
            arc_set.add((int(fields[0]),int(fields[1])))
        lineno+=1

    # return result
    return arc_set
    
    
##################################################
def load_reac_file(reacf):

    # initialize result
    reac_set=set()

    # read file line by line
    file = open(reacf, 'r')
    for line in file:
        line=line.strip("\n")
        fields=line.split(" ")
        reac_set.add(int(fields[0]))

    # return result
    return reac_set

##################################################
def load_extern_metab_file(extermf):

    # initialize result
    extern_metab_set=set()

    # read file line by line
    file = open(extermf, 'r')
    for line in file:
        line=line.strip("\n")
        fields=line.split(" ")
        extern_metab_set.add(int(fields[0]))

    # return result
    return extern_metab_set

##################################################
def filter_arc_set(arc_set,extern_metab_set):

    filt_arc_set=set()

    for (metab,reac) in arc_set:
        if(not metab in extern_metab_set):
            filt_arc_set.add((metab,reac))
    
    return filt_arc_set
    

##################################################
def obtain_metabs_for_each_reac(filt_arc_set):

    reac_dict={}

    for (metab,reac) in filt_arc_set:
        if(reac in reac_dict):
            reac_dict[reac].append(metab)
        else:
            reac_dict[reac]=[metab]

    return reac_dict
    
##################################################
def obtain_metabs_linked_to_acc_reacs(filt_arc_set,reac_set):

    acc_reacs_metabolites=set()

    for (metab,reac) in filt_arc_set:
        if(reac in reac_set):
            acc_reacs_metabolites.add(metab)
    
    return acc_reacs_metabolites
    
##################################################
def obtain_accessible_reactions(arc_set,reac_set,extern_metab_set):

    # initialize result and related variables
    accessible_reacs=set(reac_set)
    prev_accessible_reacs=set(accessible_reacs)
    
    # filter arc set given external metabolite set
    filt_arc_set=filter_arc_set(arc_set,extern_metab_set)

    # obtain metabolites for each reaction
    reac_dict=obtain_metabs_for_each_reac(filt_arc_set)
    
    # obtain metabolites linked to accessible reactions
    acc_reacs_metabolites=obtain_metabs_linked_to_acc_reacs(filt_arc_set,reac_set)

    # start loop to obtain accessible reactions
    iterno=1
    end=False
    while(not end):
        print  >> sys.stderr,"* Iter no:",iterno," , num. acc. reacs:",len(accessible_reacs)
        # iterate over filtered reactions
        for reac in reac_dict:
            # find reactions not included in accessible reactions having
            # one metabolite in the set of metabolites linked to
            # accessible reactions
            if(not reac in accessible_reacs):
                for metab in reac_dict[reac]:
                    if(metab in acc_reacs_metabolites):
                        accessible_reacs.add(reac)
                        # add metabolites of reac to the set of metabolites
                        # linked to accessible reactions
                        for met in reac_dict[reac]:
                            acc_reacs_metabolites.add(met)
                        break
        # check stop condition
        if(len(accessible_reacs)==len(prev_accessible_reacs)):
            end=True
        else:
            prev_accessible_reacs=set(accessible_reacs)
        iterno+=1

    print >> sys.stderr,"Algorithm finished, num. acc. reacs:",len(accessible_reacs)
        
    # return result
    return accessible_reacs

##################################################
def print_reacs(accessible_reacs):

    for reac in accessible_reacs:
        print reac
    
##################################################
def main(argv):
    # take parameters
    s_given=False
    r_given=False
    reacf=""
    e_given=False
    extermf=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:r:e:",["help","sbmlf=","reacf=","extermf"])
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
            elif opt in ("-r", "--reacf"):
                reacf = arg
                r_given=True
            elif opt in ("-e", "--extermf"):
                extermf = arg
                e_given=True

    # print parameters
    if(s_given==True):
        print >> sys.stderr, "s is %s" % (sbmlf)
    else:
        print >> sys.stderr, "Error: -s option not given"
        sys.exit(2)

    if(r_given==True):
        print >> sys.stderr, "r is %s" % (reacf)
    else:
        print >> sys.stderr, "Error: -r option not given"
        sys.exit(2)

    ## process parameters

    # read arc information from network
    arc_set=read_arcs(sbmlf)
    
    # read reaction ids
    reac_set=load_reac_file(reacf)

    # load external metabolite ids
    extern_metab_set=set()
    if(extermf!=""):
        extern_metab_set=load_extern_metab_file(extermf)

    # obtain accessible reactions
    accessible_reacs=obtain_accessible_reactions(arc_set,reac_set,extern_metab_set)

    # print reactions
    print_reacs(accessible_reacs)
    
if __name__ == "__main__":
    main(sys.argv)
