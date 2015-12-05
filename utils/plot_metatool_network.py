# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_help():
    print >> sys.stderr, "plot_metatool_network -f <string> [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-f <string> :    metatool model file" 
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
class Metatool_info:
    def __init__(self):
        self.enzrev=set()
        self.enzirrev=set()
        self.metint=set()
        self.metext=set()
        self.reactin={}
        self.reactout={}

##################################################
def process_cat_line(result,fields):
    result.reactin[fields[0]]=[]
    result.reactout[fields[0]]=[]
    cat_mode="in"
    i=2
    while i<len(fields)-1:
        if(fields[i]=="="):
            cat_mode="out"
            i=i+1
        else:
            if(cat_mode=="in"):
                coef=1
                if(fields[i]=="+"):
                    i=i+1
                if(fields[i]=="-"):
                    coef=-1
                    i=i+1
                if(fields[i].isdigit()):
                    coef=coef*int(fields[i])
                    i=i+1
                metab=fields[i]
                result.reactin[fields[0]].append([coef,metab])
                i=i+1
            else:
                coef=1
                if(fields[i]=="+"):
                    i=i+1
                if(fields[i]=="-"):
                    coef=-1
                    i=i+1
                if(fields[i].isdigit()):
                    coef=coef*int(fields[i])
                    i=i+1
                metab=fields[i]
                result.reactout[fields[0]].append([coef,metab])
                i=i+1

##################################################
def load_metatool_file(metatoolf):

    # initialize result
    result= Metatool_info()

    # read file line by line
    file = open(metatoolf, 'r')
    mode=""
    for line in file:
        line=line.strip("\n")
        fields=line.split(" ")            

        if(mode=="ENZREV"):
            for i in range(len(fields)):
                result.enzrev.add(fields[i])

        elif(mode=="ENZIRREV"):
            for i in range(len(fields)):
                result.enzirrev.add(fields[i])

        elif(mode=="METINT"):
            for i in range(len(fields)):
                result.metint.add(fields[i])

        elif(mode=="METEXT"):
            for i in range(len(fields)):
                result.metext.add(fields[i])

        elif(mode=="CAT"):
            if(fields[1]==":"):
                process_cat_line(result,fields)

        # Determine mode
        if(fields[0]=="-ENZREV"):
            mode="ENZREV"
        elif(fields[0]=="-ENZIRREV"):
            mode="ENZIRREV"
        elif(fields[0]=="-METINT"):
            mode="METINT"
        elif(fields[0]=="-METEXT"):
            mode="METEXT"
        elif(fields[0]=="-CAT"):
            mode="CAT"

    # return result
    return result

##################################################
def print_fba_network(metatool_info):

    # Print header
    print "digraph word_graph {"
    print "rankdir=LR;"

    ## Set representation for the different nodes

    # Set representation for reactions
    print "node [shape = box]; ",

    for key in metatool_info.reactin.keys():
        print key,";",
    print ""

    # Set representation for metabolites
    print "node [shape = circle];"

    ## Process stochiometric relations
    for (key,value) in metatool_info.reactin.iteritems():
        for n in range(len(value)):
            print value[n][1],"->",key, "[ label= \"",value[n][0],"\",","color = black ];"

    for (key,value) in metatool_info.reactout.iteritems():
        for n in range(len(value)):
            print key,"->",value[n][1], "[ label= \"",value[n][0],"\",","color = black ];"

    # Print footer
    print "}"

##################################################
def plot_network(metatoolf):

    # load metatoolf
    metatool_info=load_metatool_file(metatoolf)

    # print fba network
    print_fba_network(metatool_info)

##################################################
def main(argv):
    # take parameters
    f_given=False
    metatoolf=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:",["help","metatoolf="])
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
            elif opt in ("-f", "--metatoolf"):
                metatoolf = arg
                f_given=True

    # print parameters
    if(f_given==True):
        print >> sys.stderr, "f is %s" % (metatoolf)
    else:
        print >> sys.stderr, "Error: -f option not given"
        sys.exit(2)

    # create lp file according to selected criterion
    plot_network(metatoolf)
        
if __name__ == "__main__":
    main(sys.argv)
