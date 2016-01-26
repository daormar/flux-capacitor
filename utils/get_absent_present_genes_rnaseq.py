# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, numpy

##################################################
class rnaseqc_info:
    def __init__(self):
        self.sampleidmap={}
        self.genelist=[]
        self.samplecounts=[]

##################################################
class key_percs:
    def __init__(self):
        self.low_perc=0
        self.high_perc=0

##################################################
def load_rnaseqc_info(filename):
    rscinfo=rnaseqc_info()
    file = open(filename, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        if(lineno==1):
            # Obtain sample ids
            for i in range(1,len(fields)):
                rscinfo.sampleidmap[fields[i]]=i-1
                rscinfo.samplecounts.append([])
        else:
            # Process gene counts
            rscinfo.genelist.append(fields[0])
            for i in range(1,len(fields)):
                rscinfo.samplecounts[i-1].append(float(fields[i]))

        lineno=lineno+1

    return rscinfo

##################################################
def load_sample_list(filename):
    samplelist=[]
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        samplelist.append(line)
    return samplelist

##################################################
def proc_votes(fields):
    # Count ones and zeros
    num_zeros=0
    num_ones=0
    for i in fields[1:]:
        if(i=="0"):
            num_zeros=num_zeros+1
        elif(i=="1"):
            num_ones=num_ones+1

    # Print output entry
    if(num_zeros>num_ones):
        return "0"
    elif(num_ones>num_zeros):
        return "1"
    else:
        return "NA"  

##################################################
def obtain_percs(numlist):
    percs={}
    for i in range(1,100):
        percs[i]=numpy.percentile(numlist,i)
    return percs

##################################################
def determine_status(kp,percs,val):
    if(val<percs[kp.low_perc]):
        return "0"
    elif(val<percs[kp.high_perc]):
        return "NA"
    elif(val>=percs[kp.high_perc]):
        return "1"

##################################################
def get_abs_pres_genes_given_col(col,kp,rscinfo):
    if(col > 0 and col-1<len(rscinfo.samplecounts)):
        # Obtain percentiles
        percs=obtain_percs(rscinfo.samplecounts[col-1])

        # Process genes
        for i in range(len(rscinfo.samplecounts[col-1])):
            val=rscinfo.samplecounts[col-1][i]
            status=determine_status(kp,percs,val)
            print rscinfo.genelist[i]+","+status

##################################################
def get_abs_pres_genes_given_slist(samplelist,kp,rscinfo):
    # Obtain percentiles
    percs_list={}
    for i in range(len(samplelist)):
        sampleid=rscinfo.sampleidmap[samplelist[i]]
        percs_list[sampleid]=obtain_percs(rscinfo.samplecounts[sampleid])

    # Process genes
    for i in range(len(rscinfo.samplecounts[0])):
        statuslist=[]
        for j in range(len(samplelist)):
            sampleid=rscinfo.sampleidmap[samplelist[j]]
            val=rscinfo.samplecounts[sampleid][i]
            status=determine_status(kp,percs_list[sampleid],val)
            statuslist.append(status)

        # Obtain consensus status
        consensus_status=proc_votes(statuslist)

        # Print gene plus status
        print rscinfo.genelist[i]+","+status

##################################################
def print_help():
    print >> sys.stderr, "get_absent_present_genes_rnaseq -r <string> [-a <int>] [-b <int>]"
    print >> sys.stderr, "                                [-c <int> | -l <string>] [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-r <string> :    file with rna-seq counts"
    print >> sys.stderr, "-a <int>    :    lower percentile used to determine gene presence/absence"
    print >> sys.stderr, "-b <int>    :    higher percentile used to determine gene presence/absence"
    print >> sys.stderr, "-c <int>    :    get abs/pres data for <int>'th column"
    print >> sys.stderr, "-l <string> :    file with list of samples to be taken into account"
    print >> sys.stderr, "                 for the generation of abs/pres data"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
def main(argv):
    # take parameters
    r_given=False
    a_given=False
    a_val=30
    b_given=False
    b_val=40
    c_given=False
    l_given=False
    rnaseqcf = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"ha:b:r:c:l:",["help","lowerp=","higherp=","rnaseqcf=","col=","listf="])
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
            elif opt in ("-r", "--rnaseqcf"):
                rnaseqcf = arg
                r_given=True
            elif opt in ("-a", "--lowerp"):
                a_val = int(arg)
                a_given=True
            elif opt in ("-b", "--higherp"):
                b_val = int(arg)
                b_given=True
            elif opt in ("-c", "--col"):
                col = int(arg)
                c_given=True
            elif opt in ("-l", "--listf"):
                listf = arg
                l_given=True

    # print parameters
    if(r_given==True):
        print >> sys.stderr, "r is %s" % (rnaseqcf)
    else:
        print >> sys.stderr, "Error: -r option not given"
        sys.exit(2)

    if(c_given==True and l_given==True):
        print >> sys.stderr, "Error: -c and -l options cannot be given simultaneously"
        sys.exit(2)
        
    if(c_given==True):
        print >> sys.stderr, "c is %d" % (col)

    if(l_given==True):
        print >> sys.stderr, "l is %s" % (listf)

    # load rna-seq information
    rscinfo=load_rnaseqc_info(rnaseqcf)

    # load list of arrays if given
    if(l_given==True):
        samplelist=load_sample_list(listf)

    # Define key percentiles
    kp=key_percs()
    kp.low_perc=a_val
    kp.high_perc=b_val

    # Process absent/present genes info
    if(c_given==True):
        get_abs_pres_genes_given_col(col,kp,rscinfo)
    else:
        if(l_given==True):
            get_abs_pres_genes_given_slist(samplelist,kp,rscinfo)
        else:
            # Give the whole list of samples
            samplelist=[]
            for k in rscinfo.sampleidmap.keys():
                samplelist.append(k)
            get_abs_pres_genes_given_slist(samplelist,kp,rscinfo)
    
if __name__ == "__main__":
    main(sys.argv)
