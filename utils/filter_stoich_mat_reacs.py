# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_help():
    print >> sys.stderr, "filter_stoich_mat_reacs -s <string> -r <string>"
    print >> sys.stderr, "                        [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-s <string> :    file containing the stoichimetric matrix in csv format"
    print >> sys.stderr, "-r <string> :    file containing the reaction id's to be filtered"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

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
def filter_stoich_mat(stoichf,reac_set):

    # read file line by line
    lineno=1
    file = open(stoichf, 'r')
    for line in file:
        line=line.strip("\n")
        if(lineno<3):
            print line
        else:
            fields=line.split(" ")
            if(int(fields[1]) in reac_set):
                print line
        lineno+=1

##################################################
def main(argv):
    # take parameters
    s_given=False
    r_given=False
    reacf=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:r:",["help","stoichf=","reacf="])
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
            elif opt in ("-s", "--stoichf"):
                stoichf = arg
                s_given=True
            elif opt in ("-r", "--reacf"):
                reacf = arg
                r_given=True

    # print parameters
    if(s_given==True):
        print >> sys.stderr, "s is %s" % (stoichf)
    else:
        print >> sys.stderr, "Error: -s option not given"
        sys.exit(2)

    if(r_given==True):
        print >> sys.stderr, "r is %s" % (reacf)
    else:
        print >> sys.stderr, "Error: -r option not given"
        sys.exit(2)

    ## process parameters
    
    # read reaction ids
    reac_set=load_reac_file(reacf)
    
    # obtain accessible reactions
    filter_stoich_mat(stoichf,reac_set)
    
if __name__ == "__main__":
    main(sys.argv)
