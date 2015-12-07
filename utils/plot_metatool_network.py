# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, fba

##################################################
def print_help():
    print >> sys.stderr, "plot_metatool_network -f <string> [-o] [-t <int>] [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-f <string> :    metatool model file" 
    print >> sys.stderr, "-o          :    ommit arc label" 
    print >> sys.stderr, "-t <int>    :    plot type"
    print >> sys.stderr, "                 0 -> enzymes + metabolites (default option)"  
    print >> sys.stderr, "                 1 -> enzymes + internal metabolites"  
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
class Metatool_info:
    def __init__(self):
        self.enzrev=set()
        self.enzirrev=set()
        self.metint=set()
        self.metext=set()
        self.react_to_input_metab={}
        self.react_to_output_metab={}
        self.input_metab_to_react={}
        self.output_metab_to_react={}

##################################################
def process_cat_line(result,fields):
    # Initialize variables for reaction
    react=fields[0]
    result.react_to_input_metab[react]={}
    result.react_to_output_metab[react]={}
    cat_mode="in"
    i=2
    
    # Process entry
    while i<len(fields)-1:
        if(fields[i]=="="):
            cat_mode="out"
            i=i+1
        else:
            # Process input metabolites
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
                result.react_to_input_metab[react][metab]=coef
                if(not metab in result.input_metab_to_react.keys()):
                    result.input_metab_to_react[metab]={}
                result.input_metab_to_react[metab][react]=coef
                i=i+1
            # Process output metabolites
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
                result.react_to_output_metab[react][metab]=coef
                if(not metab in result.output_metab_to_react.keys()):
                    result.output_metab_to_react[metab]={}
                result.output_metab_to_react[metab][react]=coef
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
def print_enzymes_metab(metatool_info,int_metab_only,ommit_arc_label):

    # Print header
    print "digraph word_graph {"
    print "rankdir=LR;"
    print "splines=ortho;"
    print "K=2;"

    ## Set representation for the different nodes

    # Set representation for reactions
    print "node [shape = box]; ",
#    print "node [shape = point]; ",

    for key in metatool_info.react_to_input_metab.keys():
        print key,";",
    print ""

    # Set representation for metabolites
    print "node [shape = circle];"

    ## Process stoichiometric relations
    for (key,metabdict) in metatool_info.react_to_input_metab.iteritems():
        for metab in metabdict.keys():
            if(ommit_arc_label==0):
                arc_label=metabdict[metab]
            else:
                arc_label=""
            if(int_metab_only==0):
                print metab,"->",key, "[ label= \"",arc_label,"\",","color = black ];"
            else:
                if(metab in metatool_info.metint):
                    print metab,"->",key, "[ label= \"",arc_label,"\",","color = black ];"

    for (key,metabdict) in metatool_info.react_to_output_metab.iteritems():
        for metab in metabdict.keys():
            if(ommit_arc_label==0):
                arc_label=metabdict[metab]
            else:
                arc_label=""
            if(int_metab_only==0):
                print key,"->",metab, "[ label= \"",arc_label,"\",","color = black ];"
            else:
                if(metab in metatool_info.metint):
                    print key,"->",metab, "[ label= \"",arc_label,"\",","color = black ];"

    # Print footer
    print "}"

##################################################
def plot_network(metatoolf,plottype,ommit_arc_label):

    # load metatoolf
    metatool_info=load_metatool_file(metatoolf)

    # print fba network
    if(plottype==0):
        print_enzymes_metab(metatool_info,0,ommit_arc_label)
    elif(plottype==1):
        print_enzymes_metab(metatool_info,1,ommit_arc_label)

##################################################
def main(argv):
    # take parameters
    f_given=False
    metatoolf=""
    t_given=0
    plottype=0
    ommit_arc_label=0
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:t:o",["help","metatoolf=","plottype=","ommit-arc-label"])
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
            elif opt in ("-t", "--metatoolf"):
                plottype= int(arg)
                t_given=True
            if opt in ("-o", "--ommit-arc-label"):
                ommit_arc_label=1

    # print parameters
    if(f_given==True):
        print >> sys.stderr, "f is %s" % (metatoolf)
    else:
        print >> sys.stderr, "Error: -f option not given"
        sys.exit(2)

    if(t_given==True):
        print >> sys.stderr, "t is %s" % (plottype)

    # create lp file according to selected criterion
    plot_network(metatoolf,plottype,ommit_arc_label)
        
if __name__ == "__main__":
    main(sys.argv)
