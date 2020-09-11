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
def load_mapping_info(filename):
    mapinfo={}
    file = open(filename, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        mapinfo[fields[1]]=fields[0]
        
        lineno=lineno+1

    return mapinfo
    
##################################################
def process_file(mapinfo,filename):
    file = open(filename, 'r')
    # read file line by line
    lineno=1
    for line in file:
        line=line.strip("\n")
        fields=line.split(",")
        if(fields[0] in mapinfo):
            print(mapinfo[fields[0]])
        else:
            print("Warning: symbol",fields[0],"not present in dictionary", file=sys.stderr)
            
        lineno=lineno+1
            
##################################################
def print_help():
    print("symbol2id -m <string> -f <string> [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-m <string> :    csv file with mapping", file=sys.stderr)
    print("-f <string> :    file to be processed", file=sys.stderr)
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

##################################################
def main(argv):
    # take parameters
    m_given=False
    m_val = ""
    f_given=False
    f_val=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hm:f:",["help","map=","file="])
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
            elif opt in ("-m", "--map"):
                m_val = arg
                m_given=True
            elif opt in ("-f", "--file"):
                f_val = arg
                f_given=True

    # print parameters
    if(m_given==True):
        print("m is %s" % (m_val), file=sys.stderr)
    else:
        print("Error: -m option not given", file=sys.stderr)
        sys.exit(2)

    if(f_given==True):
        print("s is %s" % (f_val), file=sys.stderr)
    else:
        print("Error: -f option not given", file=sys.stderr)
        sys.exit(2)

    # load mapping file
    mapinfo=load_mapping_info(m_val)

    # process file
    process_file(mapinfo,f_val)
    
if __name__ == "__main__":
    main(sys.argv)
