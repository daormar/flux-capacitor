"""
Flux Capacitor package
Copyright (C) 2015-2018 Daniel Ortiz-Mart\'inez
 
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
"""
 
# *- python -*

import sys, getopt, fba

##################################################
def print_help():
    print("plot_metatool_network -f <string> [-o] [-t <int>] [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-f <string> :    metatool model file", file=sys.stderr) 
    print("-o          :    ommit arc label", file=sys.stderr) 
    print("-t <int>    :    plot type", file=sys.stderr)
    print("                 0 -> enzymes + metabolites (default option)", file=sys.stderr)  
    print("                 1 -> enzymes + internal metabolites", file=sys.stderr)  
    print("                 2 -> reactions", file=sys.stderr)
    print("                 3 -> reactions (exclude external metabolites)", file=sys.stderr)
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

##################################################
class Metatool_info:
    def __init__(self):
        self.enzrev=set()
        self.enzirrev=set()
        self.metint=set()
        self.metext=set()
        self.reactions=[]
        self.react_to_input_metab={}
        self.react_to_output_metab={}
        self.input_metab_to_react={}
        self.output_metab_to_react={}

##################################################
def process_cat_line(result,fields):
    # Initialize variables for reaction
    react=fields[0]
    result.reactions.append(react)
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
                if(not metab in list(result.input_metab_to_react.keys())):
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
                if(not metab in list(result.output_metab_to_react.keys())):
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
        fields=line.split()            

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
        if(len(fields)>0):
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
    print("digraph enzymes_metab_graph {")
    print("rankdir=LR;")
    print("overlap=false;")
    print("splines=true;")
#    print "splines=ortho;"
    print("K=1;")

    ## Set representation for the different nodes

    # Set representation for reactions
    print("node [shape = box]; ", end=' ')
#    print "node [shape = point]; ",

    for key in list(metatool_info.react_to_input_metab.keys()):
        print(key,";", end=' ')
    print("")

    # Set representation for external metabolites
    if(int_metab_only==0):
        print("node [shape = none]; ", end=' ')
        for (key,metabdict) in metatool_info.react_to_input_metab.items():
            for metab in list(metabdict.keys()):
                if(metab in metatool_info.metext):
                    print(key+"_"+metab,";", end=' ')
        for (key,metabdict) in metatool_info.react_to_output_metab.items():
            for metab in list(metabdict.keys()):
                if(metab in metatool_info.metext):
                    print(key+"_"+metab,";", end=' ')
        print("")

    # Set representation for internal metabolites
    print("node [shape = circle];")

    ## Process stoichiometric relations
    for (key,metabdict) in metatool_info.react_to_input_metab.items():
        if(key in metatool_info.enzrev):
            direction="both"
        else:
            direction="single"
        for metab in list(metabdict.keys()):
            if(ommit_arc_label==0):
                arc_label=metabdict[metab]
            else:
                arc_label=""
            if(int_metab_only==0):
                if(metab in metatool_info.metint):
                    print(metab,"->",key, "[ label= \"",arc_label,"\",","color = black , dir = ","\""+direction+"\""," ];")
                else:
                    print("{"+key+"_"+metab+" [label=\""+metab+"\"]}","->",key, "[ label= \"",arc_label,"\",","color = black , dir = ","\""+direction+"\""," ];")
            else:
                if(metab in metatool_info.metint):
                    print(metab,"->",key, "[ label= \"",arc_label,"\",","color = black, dir = ","\""+direction+"\""," ];")

    for (key,metabdict) in metatool_info.react_to_output_metab.items():
        if(key in metatool_info.enzrev):
            direction="both"
        else:
            direction="single"
        for metab in list(metabdict.keys()):
            if(ommit_arc_label==0):
                arc_label=metabdict[metab]
            else:
                arc_label=""
            if(int_metab_only==0):
                if(metab in metatool_info.metint):
                    print(key,"->",metab, "[ label= \"",arc_label,"\",","color = black , dir = ","\""+direction+"\""," ];")
                else:
                    print(key,"->","{"+key+"_"+metab+" [label=\""+metab+"\"]}", "[ label= \"",arc_label,"\",","color = black , dir = ","\""+direction+"\""," ];")
            else:
                if(metab in metatool_info.metint):
                    print(key,"->",metab, "[ label= \"",arc_label,"\",","color = black , dir = ","\""+direction+"\""," ];")

    # Print footer
    print("}")

##################################################
def print_react(metatool_info,int_metab_only,ommit_arc_label):

    # Print header
    print("graph react_graph {")
    print("rankdir=LR;")
    print("overlap=false;")
    print("splines=true;")
#    print "splines=ortho;"
    print("K=1;")

    ## Set representation for the different nodes

    # Set representation for reactions
    print("node [shape = box];")

    ## Process stoichiometric relations
    for i in range(len(metatool_info.reactions)):
        for j in range(i+1,len(metatool_info.reactions)):
            # Obtain reaction identifiers
            key1=metatool_info.reactions[i]
            key2=metatool_info.reactions[j]
            # Obtain input and output metabolites for key 1
            in_metabdict1=metatool_info.react_to_input_metab[key1]
            out_metabdict1=metatool_info.react_to_output_metab[key1]
            # Look for common metabolites
            found=0
            for m in list(in_metabdict1.keys()):
                if(int_metab_only==0 or m in metatool_info.metint):
                    if m in list(metatool_info.react_to_input_metab[key2].keys()):
                        found=1
                    if m in list(metatool_info.react_to_output_metab[key2].keys()):
                        found=1
            if(found==0):
                if(int_metab_only==0 or m in metatool_info.metint):
                    for m in list(out_metabdict1.keys()):
                        if m in list(metatool_info.react_to_input_metab[key2].keys()):
                            found=1
                        if m in list(metatool_info.react_to_output_metab[key2].keys()):
                            found=1
            if(found==1):
                print(key1,"--",key2, "[ label= \"\" ,","color = black ];")

    # Print footer
    print("}")

##################################################
def plot_network(metatoolf,plottype,ommit_arc_label):

    # load metatoolf
    metatool_info=load_metatool_file(metatoolf)

    # print fba network
    if(plottype==0):
        print_enzymes_metab(metatool_info,0,ommit_arc_label)
    elif(plottype==1):
        print_enzymes_metab(metatool_info,1,ommit_arc_label)
    elif(plottype==2):
        print_react(metatool_info,0,ommit_arc_label)
    elif(plottype==3):
        print_react(metatool_info,1,ommit_arc_label)

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
        print("f is %s" % (metatoolf), file=sys.stderr)
    else:
        print("Error: -f option not given", file=sys.stderr)
        sys.exit(2)

    if(t_given==True):
        print("t is %s" % (plottype), file=sys.stderr)

    # create lp file according to selected criterion
    plot_network(metatoolf,plottype,ommit_arc_label)
        
if __name__ == "__main__":
    main(sys.argv)
