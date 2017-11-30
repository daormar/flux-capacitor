# Author: Daniel Ortiz Mart\'inez
# *- python -*

import sys, getopt, numpy
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity

##################################################
class rnaseqc_info:
    def __init__(self):
        self.sampleidmap={}
        self.genelist=[]
        self.samplecounts=[]
        self.rawcounts=[]
        self.log_nonzero_counts=[]

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
                count=float(fields[i])
                rscinfo.rawcounts.append(count)
                if(count!=0):
                    rscinfo.log_nonzero_counts.append(numpy.log2(count))
                rscinfo.samplecounts[i-1].append(count)

        lineno=lineno+1

    return rscinfo

##################################################
def perform_kde(rscinfo):
    logcounts_arr=numpy.asarray(rscinfo.log_nonzero_counts).reshape(-1,1)
    kernel = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(logcounts_arr)
    return kernel

##################################################
def plot_density(kernel):
    x = numpy.arange(-20, 20, .1)
    scores=numpy.exp(kernel.score_samples(x.reshape(-1,1)))
    plt.plot(x,scores)
    plt.ylabel('Density')
    plt.xlabel('log2(RPKM)')
    plt.show()

##################################################
def print_help():
    print >> sys.stderr, "plot_kde_rnaseq -r <string> [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-r <string> :    file with rna-seq counts"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
def main(argv):
    # take parameters
    r_given=False
    rnaseqcf = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hr:",["help","rnaseqcf="])
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

    # print parameters
    if(r_given==True):
        print >> sys.stderr, "r is %s" % (rnaseqcf)
    else:
        print >> sys.stderr, "Error: -r option not given"
        sys.exit(2)

    # load rna-seq information
    rscinfo=load_rnaseqc_info(rnaseqcf)
    
    # Perform kernel density estimation
    kernel=perform_kde(rscinfo)

    # Plot density
    plot_density(kernel)

if __name__ == "__main__":
    main(sys.argv)
