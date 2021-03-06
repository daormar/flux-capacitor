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
from statsmodels.sandbox.stats.multicomp import multipletests

##################################################
def load_pvalues_file(filename):
    labels=[]
    pvalues=[]
    file = open(filename, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        labels.append(fields[0])
        pvalues.append(float(fields[1]))
        
        lineno=lineno+1

    return labels,pvalues
    
##################################################
def print_corrected_pvalues(labels,corrected_pvalues):
    for i in range(len(labels)):
        print(labels[i]+","+str(corrected_pvalues[i]))
            
##################################################
def print_help():
    print("correct_pvalues -p <string> [-a <float>] [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-p <string> :    file with p-values", file=sys.stderr)
    print("-a <float>  :    alpha value (0.05 by default)", file=sys.stderr)
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

##################################################
def main(argv):
    # take parameters
    p_given=False
    p_val = ""
    a_given=False
    a_val=0.05
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hp:a:",["help","pval=","alpha"])
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
            elif opt in ("-p", "--pval"):
                p_val = arg
                p_given=True
            elif opt in ("-p", "--alpha"):
                a_val = float(arg)
                a_given=True

    # Print parameters
    if(p_given==True):
        print("p is %s" % (p_val), file=sys.stderr)
    else:
        print("Error: -p option not given", file=sys.stderr)
        sys.exit(2)

    # Load file with p-values
    labels,pvalues=load_pvalues_file(p_val)

    # Obtain correction
    reject,corrected_pvalues,asidak,abonf=multipletests(pvalues,alpha=a_val,method='fdr_bh')
    
    # Print results
    print_corrected_pvalues(labels,corrected_pvalues)
    
if __name__ == "__main__":
    main(sys.argv)
