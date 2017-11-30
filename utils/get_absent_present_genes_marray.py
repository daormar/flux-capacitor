# Author: Daniel Ortiz Mart\'inez
# *- python -*

import sys, getopt

##################################################
def load_gene_dict(filename):
    genedict={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        genedict[fields[0]]=fields[1]
    return genedict

##################################################
def load_array_list(filename):
    arraylist=[]
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        arraylist.append(line)
    return arraylist

##################################################
def proc_line_votes(genedict,fields):

    # Count ones and zeros
    num_zeros=0
    num_ones=0
    for i in fields[1:]:
        if(i=="0"):
            num_zeros=num_zeros+1
        elif(i=="1"):
            num_ones=num_ones+1

    # Print output entry
    if(genedict[fields[0]]!="NA"):
        if(num_zeros>num_ones):
            print genedict[fields[0]]+",0"
        elif(num_ones>num_zeros):
            print genedict[fields[0]]+",1"                
        else:
            print genedict[fields[0]]+",NA"  

##################################################
def print_help():
    print >> sys.stderr, "get_absent_present_genes_marray -d <string> -p <string>"
    print >> sys.stderr, "                                [-c <int> | -l <string>] [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-d <string> :    file with eset gene ids to entrez ids"
    print >> sys.stderr, "-p <string> :    file with data generated by panp R library"
    print >> sys.stderr, "-c <int>    :    get abs/pres data for <int>'th column"
    print >> sys.stderr, "-l <string> :    file with list of arrays to be taken into account"
    print >> sys.stderr, "                 for the generation of abs/pres data"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
def main(argv):
    # take parameters
    d_given=False
    p_given=False
    c_given=False
    l_given=False
    dictf = ""
    panpf= ""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hd:p:c:l:",["help","dictf=","panpf=","col=","listf="])
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
            elif opt in ("-d", "--dictf"):
                dictf = arg
                d_given=True
            elif opt in ("-p", "--panpf"):
                panpf = arg
                p_given=True
            elif opt in ("-c", "--col"):
                col = int(arg)
                c_given=True
            elif opt in ("-l", "--listf"):
                listf = arg
                l_given=True


    # print parameters
    if(d_given==True):
        print >> sys.stderr, "d is %s" % (dictf)
    else:
        print >> sys.stderr, "Error: -d option not given"
        sys.exit(2)

    if(p_given==True):
        print >> sys.stderr, "p is %s" % (panpf)
    else:
        print >> sys.stderr, "Error: -p option not given"
        sys.exit(2)

    if(c_given==True and l_given==True):
        print >> sys.stderr, "Error: -c and -l options cannot be given simultaneously"
        sys.exit(2)
        
    if(c_given==True):
        print >> sys.stderr, "c is %d" % (col)

    if(l_given==True):
        print >> sys.stderr, "l is %s" % (listf)

    # load gene dictionary
    genedict=load_gene_dict(dictf)

    # load list of arrays if given
    if(l_given==True):
        arraylist=load_array_list(listf)

    # Process absent/present genes info
    file = open(panpf, 'r')
    # read file line by line
    lineno=1
    for line in file:
        # Extract line fields
        line=line.strip("\n")
        fields=line.split(",")

        # Process line
        if(lineno==1):
            array_names=[]
            array_names.append("")
            for i in range(len(fields)):
                array_names.append(fields[i])
        else:
            if(c_given==True):
                if(genedict[fields[0]]!="NA"):
                    if(col > 0 and col<len(fields)):
                        print genedict[fields[0]]+","+fields[col]
            else:
                if(l_given==True):
                    filtered_fields=[]
                    filtered_fields.append(fields[0])
                    for i in range(1,len(fields)):
                        if(array_names[i] in arraylist):
                            filtered_fields.append(fields[i])
                    proc_line_votes(genedict,filtered_fields)
                else:
                    proc_line_votes(genedict,fields)

        # Increase lineno
        lineno=lineno+1
    
if __name__ == "__main__":
    main(sys.argv)
