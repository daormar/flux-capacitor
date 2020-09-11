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
import json

##################################################
def print_react_info(decoded_file,react_file_name):
    file=open(react_file_name, 'w')
    for reactid,react in decoded_file[1]['reactions'].items():
        entry=reactid+" "+react['bigg_id']+" ; rev: "+str(react['reversibility'])+" ; coeffs:"
        for i in range(len(react['metabolites'])):
            entry=entry+" "+react['metabolites'][i]['bigg_id']+" "+str(react['metabolites'][i]['coefficient'])+" ,"
        print(entry, file=file)

##################################################
def print_arc_info(decoded_file,arc_file_name):
    file=open(arc_file_name, 'w')
    for reactid,react in decoded_file[1]['reactions'].items():
        for segmentid,segment in react['segments'].items():
            print(segment['from_node_id'],segment['to_node_id'],react['bigg_id'], file=file)

##################################################
def print_node_info(decoded_file,node_file_name):
    file=open(node_file_name, 'w')
    for nodeid,node in decoded_file[1]['nodes'].items():
        if(node['node_type']=='metabolite'):
            print(nodeid,node['bigg_id'], file=file)

##################################################
def print_help():
    print("extract_info_from_escher_json -f <string> -o <string> [--help]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-f <string> :    Escher file in json format", file=sys.stderr)
    print("-o <string> :    Prefix of output files", file=sys.stderr)
    print("--help      :    print this help message", file=sys.stderr) 
    print("", file=sys.stderr)

##################################################
def main(argv):
    # take parameters
    f_given=False
    f_val=""
    o_given=False
    o_val=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:o:",["help","file=","outpref="])
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
            elif opt in ("-f", "--file"):
                f_val = arg
                f_given=True
            elif opt in ("-o", "--outpref"):
                o_val = arg
                o_given=True

    # print parameters
    if(f_given==True):
        print("f is %s" % (f_val), file=sys.stderr)
    else:
        print("Error: -f option not given", file=sys.stderr)
        sys.exit(2)

    if(o_given==True):
        print("o is %s" % (o_val), file=sys.stderr)
    else:
        print("Error: -o option not given", file=sys.stderr)
        sys.exit(2)

    # Open file
    file = open(f_val, 'r')

    # Decode file
    decoded_file=json.load(file)

    # Print react information
    react_file_name=o_val+".reacts"
    print_react_info(decoded_file,react_file_name)

    # Print arc information
    arc_file_name=o_val+".arcs"
    print_arc_info(decoded_file,arc_file_name)

    # Print node information
    node_file_name=o_val+".nodes"
    print_node_info(decoded_file,node_file_name)    

if __name__ == "__main__":
    main(sys.argv)
