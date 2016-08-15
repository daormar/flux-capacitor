# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, numpy
import statsmodels.stats.weightstats
import scipy.stats
import math

##################################################
class phenotype_info:
    def __init__(self):
        self.sample_to_pheno={}

##################################################
def load_phenotype_info(filename):
    phenoinfo=phenotype_info()
    file = open(filename, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        if(lineno!=1):
            phenoinfo.sample_to_pheno[fields[1]]=fields[2]

        lineno=lineno+1

    return phenoinfo

##################################################
def process_sample_fluxes(phenoinfo,fluxes_file,u_opt):
    file = open(fluxes_file, 'r')
    test_results={}
    # read file line by line
    lineno=1
    first_line_fields=[]
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        if(lineno==1):
            for i in range(len(fields)):
                first_line_fields.append(fields[i])
        else:
            # Obtain cases and controls measurements
            controls=[]
            cases=[]
            react=fields[0]
            for i in range(1,len(fields)):
                sample=first_line_fields[i]
                if(phenoinfo.sample_to_pheno[sample]=="Normal"):
                    controls.append(float(fields[i]))
                else:
                    cases.append(float(fields[i]))

            if(u_opt==False):
                # Perform t-test
                tstat,pvalue,df = statsmodels.stats.weightstats.ttest_ind(numpy.asarray(controls),numpy.asarray(cases),usevar='unequal')
            else:
                # Perform Mann-Whitney U test
                try:
                    tstat,pvalue = scipy.stats.mannwhitneyu(x=numpy.asarray(controls),y=numpy.asarray(cases),alternative='two-sided')
                except ValueError:
                    tstat = 0
                    pvalue = 1
                    print >> sys.stderr, "ValueError exception raised in line",lineno
                        
            # Store result
            test_results[react]=tstat,pvalue
            
        lineno=lineno+1

    return test_results

##################################################
def print_test_results_json(test_results,tstat_opt):
    print "{"
    for react in test_results.keys():
        tstat,pvalue=test_results[react]
        if(not math.isnan(tstat) and not math.isnan(pvalue)):
            if(tstat_opt):
                print '"'+react+'": '+str(tstat)+','
            else:
                print '"'+react+'": '+str(pvalue)+','
    print "}"

##################################################
def print_test_results_csv(test_results,tstat_opt):
    for react in test_results.keys():
        tstat,pvalue=test_results[react]
        if(not math.isnan(tstat) and not math.isnan(pvalue)):
            if(tstat_opt):
                print react+","+str(tstat)
            else:
                print react+","+str(pvalue)
            
##################################################
def print_help():
    print >> sys.stderr, "test_sample_fluxes -p <string> -f <string> --tstat --json [-u] [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-p <string> :    file with phenotype data"
    print >> sys.stderr, "-f <string> :    file with sample fluxes"
    print >> sys.stderr, "--tstat     :    print test statistic instead of p-value"    
    print >> sys.stderr, "--json      :    generate output in json format instead of csv"
    print >> sys.stderr, "-u          :    perform mann-whitney U test instead of t test."    
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

##################################################
def main(argv):
    # take parameters
    p_given=False
    p_val = ""
    f_given=False
    f_val=""
    tstat_opt=False
    json_opt=False
    u_opt=False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hp:f:tju",["help","pheno=","fluxes=","tstat","json","u"])
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
            elif opt in ("-p", "--pheno"):
                p_val = arg
                p_given=True
            elif opt in ("-f", "--fluxes"):
                f_val = arg
                f_given=True
            if opt in ("-t", "--tstat"):
                tstat_opt=True
            if opt in ("-j", "--json"):
                json_opt=True
            if opt in ("-u", "--u"):
                u_opt=True

    # print parameters
    if(p_given==True):
        print >> sys.stderr, "p is %s" % (p_val)
    else:
        print >> sys.stderr, "Error: -p option not given"
        sys.exit(2)

    if(f_given==True):
        print >> sys.stderr, "f is %s" % (f_val)
    else:
        print >> sys.stderr, "Error: -f option not given"
        sys.exit(2)

    # load phenotype file
    phenoinfo=load_phenotype_info(p_val)

    # Obtain test results
    test_results=process_sample_fluxes(phenoinfo,f_val,u_opt)

    # Print results
    if(json_opt):
        print_test_results_json(test_results,tstat_opt)
    else:
        print_test_results_csv(test_results,tstat_opt)
    
if __name__ == "__main__":
    main(sys.argv)
