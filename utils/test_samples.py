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

import sys, getopt, numpy, math
import statsmodels.stats.weightstats
import scipy.stats
from json import JSONEncoder

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
def process_samples(phenoinfo,samples_file,u_opt):
    file = open(samples_file, 'r')
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

                # Store result
                test_results[react]=tstat,pvalue
            else:
                # Perform Mann-Whitney U test
                try:
                    tstat,pvalue = scipy.stats.mannwhitneyu(x=numpy.asarray(controls),y=numpy.asarray(cases),alternative='two-sided')
                    # Store result
                    test_results[react]=tstat,pvalue
                    
                except ValueError:
                    tstat = 0
                    pvalue = 1
                    print("ValueError exception raised in line",lineno,"(reaction:",react,")", file=sys.stderr)

            print(react,"; Controls:",controls,"; Cases:",cases,"; test stat:",tstat,"; p-value:",pvalue, file=sys.stderr)

        lineno=lineno+1

    return test_results

##################################################
def print_test_results_json(test_results,tstat_opt):
    processed_test_results={}
    for react in list(test_results.keys()):
        tstat,pvalue=test_results[react]
        if(not math.isnan(tstat) and not math.isnan(pvalue)):
            if(tstat_opt):
                processed_test_results[react]=str(tstat)
            else:
                processed_test_results[react]=str(pvalue)
                
    print(JSONEncoder(indent=0).encode(processed_test_results))
    
##################################################
def print_test_results_csv(test_results,tstat_opt):
    for react in list(test_results.keys()):
        tstat,pvalue=test_results[react]
        if(not math.isnan(tstat) and not math.isnan(pvalue)):
            if(tstat_opt):
                print(react+","+str(tstat))
            else:
                print(react+","+str(pvalue))
            
##################################################
def print_help():
    print("test_samples -p <string> -s <string> --tstat --json [-u] [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-p <string> :    csv file with phenotype data", file=sys.stderr)
    print("-s <string> :    csv file with sample data", file=sys.stderr)
    print("--tstat     :    print test statistic instead of p-value", file=sys.stderr)    
    print("--json      :    generate output in json format instead of csv", file=sys.stderr)
    print("-u          :    perform Mann-Whitney's U-test instead of t-test.", file=sys.stderr)    
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

##################################################
def main(argv):
    # take parameters
    p_given=False
    p_val = ""
    s_given=False
    s_val=""
    tstat_opt=False
    json_opt=False
    u_opt=False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hp:s:tju",["help","pheno=","samples=","tstat","json","u"])
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
            elif opt in ("-s", "--samples"):
                s_val = arg
                s_given=True
            if opt in ("-t", "--tstat"):
                tstat_opt=True
            if opt in ("-j", "--json"):
                json_opt=True
            if opt in ("-u", "--u"):
                u_opt=True

    # print parameters
    if(p_given==True):
        print("p is %s" % (p_val), file=sys.stderr)
    else:
        print("Error: -p option not given", file=sys.stderr)
        sys.exit(2)

    if(s_given==True):
        print("s is %s" % (s_val), file=sys.stderr)
    else:
        print("Error: -s option not given", file=sys.stderr)
        sys.exit(2)

    # load phenotype file
    phenoinfo=load_phenotype_info(p_val)

    # Obtain test results
    test_results=process_samples(phenoinfo,s_val,u_opt)

    # Print results
    if(json_opt):
        print_test_results_json(test_results,tstat_opt)
    else:
        print_test_results_csv(test_results,tstat_opt)
    
if __name__ == "__main__":
    main(sys.argv)
