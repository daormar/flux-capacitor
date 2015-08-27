# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_obj_func(hlreact_set):

    # Print header
    print "Maximize"

    # Print objective function
    for i in range(1,len(hlreact_set)):
        if(i<len(hlreact_set)-1):
            if(hlreact_set[i]==1):
                st="+ " + fba.gen_yplus_h_name(i) + " + " + fba.gen_yminus_name(i)
                print st,
            elif(hlreact_set[i]==0):
                st="+ " + fba.gen_yplus_l_name(i)
                print st,
        else:
            if(hlreact_set[i]==1):
                st="+ " + fba.gen_yplus_h_name(i) + " + " + fba.gen_yminus_name(i)
                print st
            elif(hlreact_set[i]==0):
                st="+ " + fba.gen_yplus_l_name(i)
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
        varname=fba.gen_vname(i)
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
                st=fba.gen_yplus_h_name(i) + " " + fba.gen_yminus_name(i)
                print st,
            elif(hlreact_set[i]==0):
                st=fba.gen_yplus_l_name(i)
                print st,
        else:
            if(hlreact_set[i]==1):
                st=fba.gen_yplus_h_name(i) + " " + fba.gen_yminus_name(i)
                print st
            elif(hlreact_set[i]==0):
                st=fba.gen_yplus_l_name(i)
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
        metabname=fba.clean_string(sbmli.metabmap[k])
        # Print constraint
        print "_"+metabname+":",
        for i in range(len(sbmli.stoicheqdict[k])):
            vname=fba.gen_vname(sbmli.stoicheqdict[k][i].v)
            if(sbmli.stoicheqdict[k][i].coef >= 0.0):
                print "+",sbmli.stoicheqdict[k][i].coef,vname,
            else:
                print "-",-sbmli.stoicheqdict[k][i].coef,vname,
        print "= 0"

    # Init epsilon
    epsilon=1

    # Print lower bounds for R_H
    for i in range(1,len(sbmli.rlowbndlist)):
        if(hlreact_set[i]==1):
            vname=fba.gen_vname(i)
            ypname=fba.gen_yplus_h_name(i)
            coef=sbmli.rlowbndlist[i]-epsilon
            if(coef>=0):
                print vname,"+",coef,ypname,">=",sbmli.rlowbndlist[i]
            else:
                print vname,coef,ypname,">=",sbmli.rlowbndlist[i]

    # Print upper bounds for R_H
    for i in range(1,len(sbmli.ruppbndlist)):
        if(hlreact_set[i]==1):
            vname=fba.gen_vname(i)
            ymname=fba.gen_yminus_name(i)
            coef=sbmli.ruppbndlist[i]+epsilon
            if(coef>=0):
                print vname,"+",coef,ymname,"<=",sbmli.ruppbndlist[i]
            else:
                print vname,coef,ymname,"<=",sbmli.ruppbndlist[i]

    # Print upper and lower bounds for R_L
    for i in range(1,len(sbmli.ruppbndlist)):
        if(hlreact_set[i]==0):
            vname=fba.gen_vname(i)
            ypname=fba.gen_yplus_l_name(i)
            # print vname,">=",sbmli.rlowbndlist[i],"-",ypname
            # print vname,"<=",sbmli.ruppbndlist[i],"-",ypname
            print vname,"+",ypname,">=",sbmli.rlowbndlist[i]
            print vname,"+",ypname,"<=",sbmli.ruppbndlist[i]

    # Print footer
    print ""

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

    # Print end string
    print "End"

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
    sbmli=fba.extract_sbml_info(sbmlf)

    # load absent/present genes info
    abspres_info=fba.load_abspres_info(abspresf)

    # load mapping between probeset ids and entrez ids
    idmap_info=fba.load_idmap_info(idmapf)

    # Obtain highly/lowly expressed reactions
    hlreact_set=fba.obtain_hlreact_set(sbmli,abspres_info,idmap_info)
    print >> sys.stderr,"* R_H/R_L information"
    for i in range(1,len(hlreact_set)):
        print >> sys.stderr,"%05d" % (i),hlreact_set[i]

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
