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
    print("plot_metab_network -s <string> -d <string> -f <string> [-e <string>] [-t <int>]", file=sys.stderr)
    print("                   [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-s <string> :    prefix of SBML info files", file=sys.stderr)
    print("-d <string> :    csv file with data for reaction id's (fluxes, p-values, etc.)", file=sys.stderr)
    print("-f <string> :    file containing the reaction id's to be included in the plot", file=sys.stderr)
    print("-e <string> :    file containing the metabolite id's considered to be external", file=sys.stderr)
    print("                 (if not provided, all metabolites are considered internal)", file=sys.stderr)
    print("-t <int>    :    plot type", file=sys.stderr)
    print("                 0 -> reactions + colored reaction senses + metabolites +", file=sys.stderr)
    print("                      stoichiometric coefs. (default option)", file=sys.stderr)
    print("                 1 -> reactions + colored reaction senses + metabolites", file=sys.stderr)
    print("                 2 -> reactions:value + colored reaction senses + metabolites", file=sys.stderr)
    print("                 3 -> reactions:value + colored reaction p-value + metabolites", file=sys.stderr)
    print("                 4 -> all reactions (selected labels) + all metabolites", file=sys.stderr)
    print("                      (no labels)", file=sys.stderr)
    print("                 5 -> all reactions (no labels) + all metabolites (no labels)", file=sys.stderr)
    print("                 6 -> all reactions (all labels) + all metabolites (no labels)", file=sys.stderr)
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

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
        value=float(fields[2])
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
    file = open(extermf, 'r')
    for line in file:
        line=line.strip("\n")
        fields=line.split(" ")
        extern_metab_set.add(int(fields[0]))

    # return result
    return extern_metab_set

##################################################
def assign_color(react_data,rid,crit):

    if(rid in react_data):
        value=react_data[rid]
    else:
        value=None
        
    # Set value of result variable
    if(crit==0):
        # Color assignment criterion for flux values
        if(value==None):
            result="gray"
        elif(value>-1 and value<1):
            result="gray"
        elif(value>=1):
            result="red"
        elif(value<=-1):
            result="blue"
    elif(crit==1):
        # Color assignment criterion for p-values
        if(value==None):
            result="blue"  
        elif(value>0.05):
            result="blue"
        elif(value<=0.05):
            result="red"
            
    # Return result
    return result


##################################################
def print_header():

    print("digraph metab_network {")
#    print "rankdir=LR;"
    print("overlap=false;")
    print("splines=true;")
#    print "splines=ortho;"
    print("K=1;")

##################################################
def print_header_large_graphs():

    print("graph metab_network {")
    print('size="10,10";')
    print('ratio="fill";')

##################################################
def reaction_representation(sbmli,included_rids,shape):

    print("node [shape =",shape,"]; ", end=' ')
    for rid in sbmli.rlowbndmap:
        if(rid in included_rids):
            vname=fba.gen_vname(rid)
            reactname=sbmli.reactmap[rid]
            print("_"+fba.clean_string(reactname),";", end=' ')
    print("")

##################################################
def metab_representation(shape):

    print("node [shape =",shape,"];")

##################################################
def gen_node_id_for_metab(extern_metab_set,mid,metabname,reactname):

    if(mid in extern_metab_set):
        return metabname+"_"+reactname
    else:
        return metabname

##################################################
def print_arc_zero(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids):
    if(rid_in_inclrids==True):
        # Obtain reaction and metabolite names
        metabname=sbmli.metabmap[mid]
        reactname=sbmli.reactmap[rid]
        clreactname=fba.clean_string(reactname)
        clmetabname=fba.clean_string(metabname)

        # Determine arc color
        color=assign_color(reactdata,rid,0)

        # Obtain node string for metabolite and reaction
        nodeid=gen_node_id_for_metab(extern_metab_set,mid,clmetabname,clreactname)
        metabnode_string="{"+"_"+nodeid +" [label=\""+metabname+"\"]}"
        reactnode_string="{"+"_"+clreactname +" [label=\""+reactname+"\"]}"

        # Print arc
        if(stoich_coef>=0):
            print(reactnode_string,"->", metabnode_string,"[ label = \"",stoich_coef,"\", color =",color,"];")
        else:
            print(metabnode_string,"->",reactnode_string,"[ label = \"",stoich_coef,"\", color =",color,"];")

##################################################
def print_arc_one(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids):
    if(rid_in_inclrids==True):
        # Obtain reaction and metabolite names
        metabname=sbmli.metabmap[mid]
        reactname=sbmli.reactmap[rid]
        clreactname=fba.clean_string(reactname)
        clmetabname=fba.clean_string(metabname)

        # Determine arc color
        color=assign_color(reactdata,rid,0)

        # Obtain node string for metabolite and reaction
        nodeid=gen_node_id_for_metab(extern_metab_set,mid,clmetabname,clreactname)
        metabnode_string="{"+"_"+nodeid +" [label=\""+metabname+"\"]}"
        reactnode_string="{"+"_"+clreactname +" [xlabel=< <b>"+reactname+"</b> >]}"

        # Print arc
        if(stoich_coef>=0):
            print(reactnode_string,"->",metabnode_string,"[ color =",color,", penwidth = 3 ];")
        else:
            print(metabnode_string,"->",reactnode_string,"[ color =",color,", penwidth = 3 ];")

##################################################
def print_arc_two(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids):
    if(rid_in_inclrids==True):
        # Obtain reaction and metabolite names
        metabname=sbmli.metabmap[mid]
        reactname=sbmli.reactmap[rid]
        clreactname=fba.clean_string(reactname)
        clmetabname=fba.clean_string(metabname)

        # Determine arc color
        color=assign_color(reactdata,rid,0)

        # Obtain node string for metabolite and reaction
        nodeid=gen_node_id_for_metab(extern_metab_set,mid,clmetabname,clreactname)
        metabnode_string="{"+"_"+nodeid +" [xlabel=\""+metabname+"\", color= darkorange ]}"
        if(rid in reactdata):
            reactnode_string="{"+"_"+clreactname +" [xlabel=< <b>"+reactname+":"+format(reactdata[rid],'.1f')+"</b> >]}"
        else:
            reactnode_string="{"+"_"+clreactname +" [xlabel=< <b>"+reactname+"</b> >]}"

        # Print arc
        if(stoich_coef>=0):
            print(reactnode_string,"->",metabnode_string,"[ color =",color,", penwidth = 3 ];")
        else:
            print(metabnode_string,"->",reactnode_string,"[ color =",color,", penwidth = 3 ];")

##################################################
def print_arc_three(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids):
    if(rid_in_inclrids==True):
        # Obtain reaction and metabolite names
        metabname=sbmli.metabmap[mid]
        reactname=sbmli.reactmap[rid]
        clreactname=fba.clean_string(reactname)
        clmetabname=fba.clean_string(metabname)

        # Determine arc color
        color=assign_color(reactdata,rid,1)

        # Obtain node string for metabolite and reaction
        nodeid=gen_node_id_for_metab(extern_metab_set,mid,clmetabname,clreactname)
        metabnode_string="{"+"_"+nodeid +" [xlabel=\""+metabname+"\", color= darkorange ]}"
        if(rid in reactdata):
            reactnode_string="{"+"_"+clreactname +" [xlabel=< <b>"+reactname+":"+format(reactdata[rid],'.3g')+"</b> >]}"
        else:
            reactnode_string="{"+"_"+clreactname +" [xlabel=< <b>"+reactname+"</b> >]}"

        # Print arc
        if(stoich_coef>=0):
            print(reactnode_string,"->",metabnode_string,"[ color =",color,", penwidth = 3 ];")
        else:
            print(metabnode_string,"->",reactnode_string,"[ color =",color,", penwidth = 3 ];")

##################################################
def print_arc_four(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids):
    # Obtain reaction and metabolite names
    metabname=sbmli.metabmap[mid]
    reactname=sbmli.reactmap[rid]
    clreactname=fba.clean_string(reactname)
    clmetabname=fba.clean_string(metabname)

    # Determine arc color
    if(rid_in_inclrids):
        color=assign_color(reactdata,rid,1)
    else:
        color="gray"

    # Obtain node string for metabolite and reaction
    nodeid=gen_node_id_for_metab(extern_metab_set,mid,clmetabname,clreactname)
    metabnode_string="{"+"_"+nodeid +" [label=\"\", color= darkorange ]}"
    if(rid_in_inclrids and rid in reactdata):
        reactnode_string="{"+"_"+clreactname +" [xlabel=< <b>"+reactname+":"+format(reactdata[rid],'.3g')+"</b> >]}"
    else:
        reactnode_string="{"+"_"+clreactname +" [label=\"\" ]}"

    # Print arc
    if(stoich_coef>=0):
        print(reactnode_string,"--",metabnode_string,"[ color =",color,", penwidth = 3 ];")
    else:
        print(metabnode_string,"--",reactnode_string,"[ color =",color,", penwidth = 3 ];")

##################################################
def print_arc_five(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids):
    # Obtain reaction and metabolite names
    metabname=sbmli.metabmap[mid]
    reactname=sbmli.reactmap[rid]
    clreactname=fba.clean_string(reactname)
    clmetabname=fba.clean_string(metabname)

    # Determine arc color
    if(rid_in_inclrids):
        color=assign_color(reactdata,rid,1)
    else:
        color="gray"

    # Obtain node string for metabolite and reaction
    nodeid=gen_node_id_for_metab(extern_metab_set,mid,clmetabname,clreactname)
    metabnode_string="{"+"_"+nodeid +" [label=\"\", color= darkorange ]}"
    reactnode_string="{"+"_"+clreactname +" [label=\"\" ]}"

    # Print arc
    if(stoich_coef>=0):
        print(reactnode_string,"--",metabnode_string,"[ color =",color,", penwidth = 3 ];")
    else:
        print(metabnode_string,"--",reactnode_string,"[ color =",color,", penwidth = 3 ];")

##################################################
def print_arc_six(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids):
    # Obtain reaction and metabolite names
    metabname=sbmli.metabmap[mid]
    reactname=sbmli.reactmap[rid]
    clreactname=fba.clean_string(reactname)
    clmetabname=fba.clean_string(metabname)

    # Determine arc color
    if(rid_in_inclrids):
        color=assign_color(reactdata,rid,1)
    else:
        color="gray"

    # Obtain node string for metabolite and reaction
    nodeid=gen_node_id_for_metab(extern_metab_set,mid,clmetabname,clreactname)
    metabnode_string="{"+"_"+nodeid +" [label=\"\", color= darkorange ]}"
    if(rid in reactdata):
        reactnode_string="{"+"_"+clreactname +" [xlabel=< <b>"+reactname+":"+format(reactdata[rid],'.3g')+"</b> >]}"
    else:
        reactnode_string="{"+"_"+clreactname +" [label=\"\" ]}"

    # Print arc
    if(stoich_coef>=0):
        print(reactnode_string,"--",metabnode_string,"[ color =",color,", penwidth = 3 ];")
    else:
        print(metabnode_string,"--",reactnode_string,"[ color =",color,", penwidth = 3 ];")

##################################################
def process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation):
    # Iterate over metabolites
    for mid in sbmli.metabmap:
        if(mid in sbmli.stoicheqdict):
            for i in range(len(sbmli.stoicheqdict[mid])):
                rid=sbmli.stoicheqdict[mid][i].v
                if(rid in included_rids):
                    rid_in_inclrids=True;
                else:
                    rid_in_inclrids=False;
                stoich_coef=sbmli.stoicheqdict[mid][i].coef
                if(arc_representation==0):
                    print_arc_zero(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids)
                elif(arc_representation==1):
                    print_arc_one(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids)
                elif(arc_representation==2):
                    print_arc_two(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids)
                elif(arc_representation==3):
                    print_arc_three(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids)
                elif(arc_representation==4):
                    print_arc_four(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids)
                elif(arc_representation==5):
                    print_arc_five(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids)
                elif(arc_representation==6):
                    print_arc_six(sbmli,extern_metab_set,reactdata,stoich_coef,mid,rid,rid_in_inclrids)

##################################################
def print_footer():
    print("}")

##################################################
def print_metab_network_type_zero(sbmli,extern_metab_set,reactdata,included_rids):

    # Print header
    print_header()

    ## Set representation for the different nodes

    # Set representation for reactions
    reaction_representation(sbmli,included_rids,"box")

    # Set representation for metabolites
    metab_representation("circle")

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
    reaction_representation(sbmli,included_rids,"point")

    # Set representation for metabolites
    metab_representation("none")

    ## Process stochiometric relations
    arc_representation=1
    process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation)

    # Print footer
    print_footer()

##################################################
def print_metab_network_type_two(sbmli,extern_metab_set,reactdata,included_rids):

    # Print header
    print_header()

    ## Set representation for the different nodes

    # Set representation for reactions
    reaction_representation(sbmli,included_rids,"point")

    # Set representation for metabolites
    metab_representation("point")

    ## Process stochiometric relations
    arc_representation=2
    process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation)

    # Print footer
    print_footer()

##################################################
def print_metab_network_type_three(sbmli,extern_metab_set,reactdata,included_rids):

    # Print header
    print_header()

    ## Set representation for the different nodes

    # Set representation for reactions
    reaction_representation(sbmli,included_rids,"point")

    # Set representation for metabolites
    metab_representation("point")

    ## Process stochiometric relations
    arc_representation=3
    process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation)

    # Print footer
    print_footer()

##################################################
def print_metab_network_type_four(sbmli,extern_metab_set,reactdata,included_rids):

    # Print header
    print_header_large_graphs()

    ## Set representation for the different nodes

    # Set representation for reactions
    reaction_representation(sbmli,included_rids,"point")

    # Set representation for metabolites
    metab_representation("point")

    ## Process stochiometric relations
    arc_representation=4
    process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation)

    # Print footer
    print_footer()

##################################################
def print_metab_network_type_five(sbmli,extern_metab_set,reactdata,included_rids):

    # Print header
    print_header_large_graphs()

    ## Set representation for the different nodes

    # Set representation for reactions
    reaction_representation(sbmli,included_rids,"point")

    # Set representation for metabolites
    metab_representation("point")

    ## Process stochiometric relations
    arc_representation=5
    process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation)

    # Print footer
    print_footer()

##################################################
def print_metab_network_type_six(sbmli,extern_metab_set,reactdata,included_rids):

    # Print header
    print_header_large_graphs()

    ## Set representation for the different nodes

    # Set representation for reactions
    reaction_representation(sbmli,included_rids,"point")

    # Set representation for metabolites
    metab_representation("point")

    ## Process stochiometric relations
    arc_representation=6
    process_stoich_relations(sbmli,extern_metab_set,reactdata,included_rids,arc_representation)

    # Print footer
    print_footer()

##################################################
def plot_network(sbmlf,dataf,filterf,extermf,pltype):
    # load sbml info
    sbmli=fba.extract_sbml_info(sbmlf)

    # read reaction data from csv file
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
    elif(pltype==2):
        print_metab_network_type_two(sbmli,extern_metab_set,reactdata,included_rids)
    elif(pltype==3):
        print_metab_network_type_three(sbmli,extern_metab_set,reactdata,included_rids)
    elif(pltype==4):
        print_metab_network_type_four(sbmli,extern_metab_set,reactdata,included_rids)
    elif(pltype==5):
        print_metab_network_type_five(sbmli,extern_metab_set,reactdata,included_rids)
    elif(pltype==6):
        print_metab_network_type_six(sbmli,extern_metab_set,reactdata,included_rids)

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
        print("s is %s" % (sbmlf), file=sys.stderr)
    else:
        print("Error: -s option not given", file=sys.stderr)
        sys.exit(2)

    if(d_given==True):
        print("d is %s" % (dataf), file=sys.stderr)
    else:
        print("Error: -d option not given", file=sys.stderr)
        sys.exit(2)

    if(f_given==True):
        print("f is %s" % (filterf), file=sys.stderr)
    else:
        print("Error: -f option not given", file=sys.stderr)
        sys.exit(2)

    if(e_given==True):
        print("e is %s" % (extermf), file=sys.stderr)

    # create lp file according to selected criterion
    plot_network(sbmlf,dataf,filterf,extermf,pltype)
        
if __name__ == "__main__":
    main(sys.argv)
