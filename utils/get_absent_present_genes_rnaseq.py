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

import sys, getopt, numpy
from sklearn import mixture

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
class key_percs:
    def __init__(self):
        self.low_perc=0
        self.high_perc=0

##################################################
def fit_gaussian_mix_model(rscinfo):
    print("Fitting gaussian mixture model with 2 components...", file=sys.stderr)
    gmm = mixture.GaussianMixture(n_components=2,covariance_type='full',n_init=1)
    logcounts_arr=numpy.asarray(rscinfo.log_nonzero_counts).reshape(-1,1)
    gmm.fit(logcounts_arr)
    print("Mixture model means:", gmm.means_[0][0],gmm.means_[1][0], file=sys.stderr)
    return gmm

##################################################
def plot_mixture(gmm):
    x = numpy.arange(-20, 20, .1)
    scores = gmm.score(x.reshape(-1,1))
    pdf = numpy.exp(scores)
    plt.plot(x,pdf)
    plt.show()

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
def mixtureid_to_exprlevel(gmm,mixid):
    if(gmm.means_[0][0]<gmm.means_[1][0]):
        return mixid
    else:
        if(mixid==0):
            return 1
        else:
            return 0

##################################################
def determine_status(gmm,val):
    if(val==0):
        return "0"
    else:
        logval=numpy.log2(val)
        mixid=gmm.predict(numpy.asarray(logval).reshape(-1,1))
        expr_level=mixtureid_to_exprlevel(gmm,mixid[0])
        return str(expr_level)

##################################################
def get_abs_pres_genes_given_col(col,gmm,rscinfo):
    if(col > 0 and col-1<len(rscinfo.samplecounts)):
        # Process genes
        for i in range(len(rscinfo.samplecounts[col-1])):
            val=rscinfo.samplecounts[col-1][i]
            status=determine_status(gmm,val)
            print(rscinfo.genelist[i]+","+status)

##################################################
def get_abs_pres_genes_given_slist(samplelist,gmm,rscinfo):
    # Process genes
    for i in range(len(rscinfo.samplecounts[0])):
        statuslist=[]
        for j in range(len(samplelist)):
            sampleid=rscinfo.sampleidmap[samplelist[j]]
            val=rscinfo.samplecounts[sampleid][i]
            status=determine_status(gmm,val)
            statuslist.append(status)

        # Obtain consensus status
        consensus_status=proc_votes(statuslist)

        # Print gene plus status
        print(rscinfo.genelist[i]+","+status)

##################################################
def print_help():
    print("get_absent_present_genes_rnaseq -r <string>", file=sys.stderr)
    print("                                [-c <int> | -l <string>] [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-r <string> :    file with rna-seq counts", file=sys.stderr)
    print("-c <int>    :    get abs/pres data for <int>'th column", file=sys.stderr)
    print("-l <string> :    file with list of samples to be taken into account", file=sys.stderr)
    print("                 for the generation of abs/pres data", file=sys.stderr)
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

##################################################
def main(argv):
    # take parameters
    r_given=False
    c_given=False
    l_given=False
    rnaseqcf = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hr:c:l:",["help","rnaseqcf=","col=","listf="])
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
            elif opt in ("-c", "--col"):
                col = int(arg)
                c_given=True
            elif opt in ("-l", "--listf"):
                listf = arg
                l_given=True

    # print parameters
    if(r_given==True):
        print("r is %s" % (rnaseqcf), file=sys.stderr)
    else:
        print("Error: -r option not given", file=sys.stderr)
        sys.exit(2)

    if(c_given==True and l_given==True):
        print("Error: -c and -l options cannot be given simultaneously", file=sys.stderr)
        sys.exit(2)
        
    if(c_given==True):
        print("c is %d" % (col), file=sys.stderr)

    if(l_given==True):
        print("l is %s" % (listf), file=sys.stderr)

    # load rna-seq information
    rscinfo=load_rnaseqc_info(rnaseqcf)

    # load sample list if given
    if(l_given==True):
        samplelist=load_sample_list(listf)

    # Fit gaussian mixture model
    gmm=fit_gaussian_mix_model(rscinfo)
        
    # Process absent/present genes info
    if(c_given==True):
        get_abs_pres_genes_given_col(col,gmm,rscinfo)
    else:
        if(l_given==True):
            get_abs_pres_genes_given_slist(samplelist,gmm,rscinfo)
        else:
            # Give the whole list of samples
            samplelist=[]
            for k in list(rscinfo.sampleidmap.keys()):
                samplelist.append(k)
            get_abs_pres_genes_given_slist(samplelist,gmm,rscinfo)
    
if __name__ == "__main__":
    main(sys.argv)
