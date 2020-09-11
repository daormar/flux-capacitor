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

import sys, getopt

##################################################
def show_data(csvfname,col):
    csvf = open(csvfname, 'r')
    # read file line by line
    lineno=1
    for line in csvf:
        line=line.strip("\n")
        fields=line.split(",")
        if(lineno==1):
            # Obtain sample ids
            for i in range(1,len(fields)):
                if(col==fields[i]):
                    colnum=i
        else:
            # Print data
            print(fields[0],fields[colnum])

        lineno=lineno+1

##################################################
class key_percs:
    def __init__(self):
        self.low_perc=0
        self.high_perc=0

##################################################
def print_help():
    print("show_data_for_column -f <string> -c <string> [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-f <string> :    csv file", file=sys.stderr)
    print("-c <int>    :    column name", file=sys.stderr)
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

##################################################
def main(argv):
    # take parameters
    f_given=False
    c_given=False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:c:",["help","csvfname=","col="])
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
            elif opt in ("-f", "--csvfname"):
                csvfname = arg
                f_given=True
            elif opt in ("-c", "--col"):
                col = arg
                c_given=True

    # print parameters
    if(f_given==True):
        print("f is %s" % (csvfname), file=sys.stderr)
    else:
        print("Error: -f option not given", file=sys.stderr)
        sys.exit(2)

    if(c_given==True):
        print("c is %s" % (col), file=sys.stderr)
    else:
        print("Error: -c option not given", file=sys.stderr)
        sys.exit(2)

    # process parameters
    show_data(csvfname,col)
    
if __name__ == "__main__":
    main(sys.argv)
