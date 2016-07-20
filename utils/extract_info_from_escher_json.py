# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
import sys, getopt, numpy
import json

##################################################
def print_react_info(decoded_file,react_file_name):
    file=open(react_file_name, 'w')
    for reactid,react in decoded_file[1]['reactions'].iteritems():
        entry=reactid+" "+react['bigg_id']+" ; rev: "+str(react['reversibility'])+" ; coeffs:"
        for i in range(len(react['metabolites'])):
            entry=entry+" "+react['metabolites'][i]['bigg_id']+" "+str(react['metabolites'][i]['coefficient'])+" ,"
        print >> file, entry

##################################################
def print_arc_info(decoded_file,arc_file_name):
    file=open(arc_file_name, 'w')
    for reactid,react in decoded_file[1]['reactions'].iteritems():
        for segmentid,segment in react['segments'].iteritems():
            print >> file, segment['from_node_id'],segment['to_node_id'],react['bigg_id']

##################################################
def print_node_info(decoded_file,node_file_name):
    file=open(node_file_name, 'w')
    for nodeid,node in decoded_file[1]['nodes'].iteritems():
        if(node['node_type']=='metabolite'):
            print >> file, nodeid,node['bigg_id']

##################################################
def print_help():
    print >> sys.stderr, "extract_info_from_escher_json -f <string> -o <string> [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-f <string> :    Escher file in json format"
    print >> sys.stderr, "-o <string> :    Prefix of output files"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

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
        print >> sys.stderr, "f is %s" % (f_val)
    else:
        print >> sys.stderr, "Error: -f option not given"
        sys.exit(2)

    if(o_given==True):
        print >> sys.stderr, "o is %s" % (o_val)
    else:
        print >> sys.stderr, "Error: -o option not given"
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
