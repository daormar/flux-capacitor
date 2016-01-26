# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_biomass_obj_func(sbmlf):

    # Print header
    print "Maximize"
    
    # Print biomass objective function
    for i in range(len(sbmlf.objfun)):
        if(i<len(sbmlf.objfun)-1):
            print fba.gen_vname(sbmlf.objfun[i]),"+",
        else:
            print fba.gen_vname(sbmlf.objfun[i])

    # Print footer
    print ""

##################################################
def print_shlomi_obj_func(hlreact_set):

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
def print_flux_boundaries(sbmli):

    # Print header
    print "Bounds"

    # Print flux upper and lower bounds
    for i in sbmli.rlowbndmap:
        varname=fba.gen_vname(i)
        print sbmli.rlowbndmap[i],"<=",varname,"<=",sbmli.ruppbndmap[i]

    # Print footer
    print ""

##################################################
def print_bin_vars(hlreact_set):

    # Print header
    print "Binary"

    # Print binary variables
    for i in range(1,len(hlreact_set)):
        if(hlreact_set[i]==1):
            st=fba.gen_yplus_h_name(i) + " " + fba.gen_yminus_name(i)
            print st
        elif(hlreact_set[i]==0):
            st=fba.gen_yplus_l_name(i)
            print st
        # elif(hlreact_set[i]==0.5):
        #     print ""

    # Print footer
    print ""

##################################################
def print_st_constraints(sbmli):
    
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

    # Print footer
    print ""
    
##################################################
def print_shlomi_constraints(sbmli,hlreact_set):

    # Init epsilon
    epsilon=1

    # Print lower bounds for R_H
    for i in sbmli.rlowbndmap:
        if(hlreact_set[i]==1):
            vname=fba.gen_vname(i)
            ypname=fba.gen_yplus_h_name(i)
            coef=sbmli.rlowbndmap[i]-epsilon
            if(coef>=0):
                print vname,"+",coef,ypname,">=",sbmli.rlowbndmap[i]
            else:
                print vname,coef,ypname,">=",sbmli.rlowbndmap[i]

    # Print upper bounds for R_H
    for i in sbmli.ruppbndmap:
        if(hlreact_set[i]==1):
            vname=fba.gen_vname(i)
            ymname=fba.gen_yminus_name(i)
            coef=sbmli.ruppbndmap[i]+epsilon
            if(coef>=0):
                print vname,"+",coef,ymname,"<=",sbmli.ruppbndmap[i]
            else:
                print vname,coef,ymname,"<=",sbmli.ruppbndmap[i]

    # Print upper and lower bounds for R_L
    for i in sbmli.ruppbndmap:
        if(hlreact_set[i]==0):
            vname=fba.gen_vname(i)
            ypname=fba.gen_yplus_l_name(i)
            if(sbmli.rlowbndmap[i]>0):
                print vname,"+",sbmli.rlowbndmap[i],ypname,">=",sbmli.rlowbndmap[i]
            elif(sbmli.rlowbndmap[i]<0):
                print vname,sbmli.rlowbndmap[i],ypname,">=",sbmli.rlowbndmap[i]
            else:
                print vname,">=",sbmli.rlowbndmap[i]

            if(sbmli.ruppbndmap[i]>0):
                print vname,"+",sbmli.ruppbndmap[i],ypname,"<=",sbmli.ruppbndmap[i]
            elif(sbmli.ruppbndmap[i]<0):
                print vname,sbmli.ruppbndmap[i],ypname,"<=",sbmli.ruppbndmap[i]
            else:
                print vname,"<=",sbmli.ruppbndmap[i]

    # Print footer
    print ""

##################################################
def print_shlomi_lp_problem(sbmli,hlreact_set):
    
    # Print objective function
    print_shlomi_obj_func(hlreact_set)

    # Print constraints
    print_st_constraints(sbmli)
    print_shlomi_constraints(sbmli,hlreact_set)

    # Print flux boundaries
    print_flux_boundaries(sbmli)

    # Print ids of binary variables
    print_bin_vars(hlreact_set)

    # Print end string
    print "End"

##################################################
def print_biomass_obj_func_fva_templ():
    print "<GOAL>"
    print "<VAR>"
    print ""

##################################################
def print_shlomi_obj_func_fva_templ():
    print "<GOAL>"
    print "<VAR>"
    print ""

##################################################
def print_biomass_constraints_fva_templ(sbmli):
    # Print fba constraints
    print_st_constraints(sbmli)

    # Print fva specific constraint
    print "fba_obj_val:",

    # Print biomass objective function
    for i in range(len(sbmli.objfun)):
        if(i<len(sbmli.objfun)-1):
            print fba.gen_vname(sbmli.objfun[i]),"+",
        else:
            print fba.gen_vname(sbmli.objfun[i]),
    print ">= <RH>"
    print ""

##################################################
def print_shlomi_fva_constraint(hlreact_set):
    print "fba_obj_val:",
    for i in range(1,len(hlreact_set)):
        if(hlreact_set[i]==1):
            st="+ " + fba.gen_yplus_h_name(i) + " + " + fba.gen_yminus_name(i)
            print st,
        elif(hlreact_set[i]==0):
            st="+ " + fba.gen_yplus_l_name(i)
            print st,

    print ">= <RH>"
    print ""

##################################################
def print_shlomi_constraints_fva_templ(sbmli,hlreact_set):
    # Print fba constraints
    print_st_constraints(sbmli)
    print_shlomi_constraints(sbmli,hlreact_set)

    # Print fva specific constraint
    print_shlomi_fva_constraint(hlreact_set)

##################################################
def print_shlomi_lp_problem_fva_templ(sbmli,hlreact_set):
    
    # Print objective function
    print_shlomi_obj_func_fva_templ()

    # Print constraints
    print_shlomi_constraints_fva_templ(sbmli,hlreact_set)

    # Print flux boundaries
    print_flux_boundaries(sbmli)

    # Print ids of binary variables
    print_bin_vars(hlreact_set)

    # Print end string
    print "End"

##################################################
def print_help():
    print >> sys.stderr, "create_lp_file -s <string> -c <int> [-a <string>] [--fva]"
    print >> sys.stderr, "               [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-s <string> :  prefix of SBML info files"
    print >> sys.stderr, "-c <int>    :  fba criterion used to generate the lp file. The criterion"
    print >> sys.stderr, "               can be selected from the following list,"    
    print >> sys.stderr, "               0 -> Maximize biomass"    
    print >> sys.stderr, "               1 -> Shlomi et al. 2008"    
    print >> sys.stderr, "-a <string> :  file with absent/present genes data (required by criterion 1)"
    print >> sys.stderr, "--fva       :  generate template file for fva instead of an fba file"
    print >> sys.stderr, "--help      :  print this help message" 
    print >> sys.stderr, ""

##################################################
def create_biomass_lp_file(sbmlf):
    # load sbml info
    sbmli=fba.extract_sbml_info(sbmlf)

    # Print objective function
    print_biomass_obj_func(sbmli)

    # Print constraints
    print_st_constraints(sbmli)

    # Print flux boundaries
    print_flux_boundaries(sbmli)

    # Print end string
    print "End"

##################################################
def create_biomass_lp_file_fva_templ(sbmlf):
    # load sbml info
    sbmli=fba.extract_sbml_info(sbmlf)

    # Print objective function
    print_biomass_obj_func_fva_templ()

    # Print constraints
    print_biomass_constraints_fva_templ(sbmli)

    # Print flux boundaries
    print_flux_boundaries(sbmli)

    # Print end string
    print "End"

##################################################
def create_shlomi_lp_file(sbmlf,abspresf):
    # load sbml info
    sbmli=fba.extract_sbml_info(sbmlf)

    # load absent/present genes info
    abspres_info=fba.load_abspres_info(abspresf)

    # Obtain highly/lowly expressed reactions
    hlreact_set=fba.obtain_hlreact_set(sbmli,abspres_info)
    print >> sys.stderr,"* R_H/R_L information"
    for i in range(1,len(hlreact_set)):
        print >> sys.stderr,"%05d" % (i),hlreact_set[i]

    # print problem in lp format
    print_shlomi_lp_problem(sbmli,hlreact_set)

##################################################
def create_shlomi_lp_file_fva_templ(sbmlf,abspresf):
    # load sbml info
    sbmli=fba.extract_sbml_info(sbmlf)

    # load absent/present genes info
    abspres_info=fba.load_abspres_info(abspresf)

    # Obtain highly/lowly expressed reactions
    hlreact_set=fba.obtain_hlreact_set(sbmli,abspres_info)

    # print problem in lp format
    print_shlomi_lp_problem_fva_templ(sbmli,hlreact_set)

##################################################
def main(argv):
    # take parameters
    s_given=False
    a_given=False
    sbmlf = ""
    abspresf= ""
    c_given=False
    crit=0
    fva=False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:a:c:f",["help","sbmlf=","abspresf=","crit=","fva"])
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
            if opt in ("-f", "--fva"):
                fva=True
            elif opt in ("-s", "--sbmlf"):
                sbmlf = arg
                s_given=True
            elif opt in ("-a", "--abspresf"):
                abspresf = arg
                a_given=True
            elif opt in ("-c", "--crit"):
                crit = int(arg)
                c_given=True

    # print parameters
    if(s_given==True):
        print >> sys.stderr, "s is %s" % (sbmlf)
    else:
        print >> sys.stderr, "Error: -s option not given"
        sys.exit(2)

    if(crit==1):
        if(a_given==True):
            print >> sys.stderr, "a is %s" % (abspresf)
        else:
            print >> sys.stderr, "Error: -a option not given"
            sys.exit(2)

    print >> sys.stderr, "c is %s" % (crit)

    print >> sys.stderr, "fva flag is %s" % (fva)

    # create lp file according to selected criterion
    if(crit==0):
        if(fva==False):
            create_biomass_lp_file(sbmlf)
        else:
            create_biomass_lp_file_fva_templ(sbmlf)
    elif(crit==1):
        if(fva==False):
            create_shlomi_lp_file(sbmlf,abspresf)
        else:
            create_shlomi_lp_file_fva_templ(sbmlf,abspresf)
        
if __name__ == "__main__":
    main(sys.argv)
