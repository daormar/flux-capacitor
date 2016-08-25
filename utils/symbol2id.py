# Author: Daniel Ortiz Mart\'inez
# *- python -*

# import modules
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
        print mapinfo[fields[0]]
        
        lineno=lineno+1
            
##################################################
def print_help():
    print >> sys.stderr, "symbol2id -m <string> -f <string> [--help]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-m <string> :    csv file with mapping"
    print >> sys.stderr, "-f <string> :    file to be processed"
    print >> sys.stderr, "--help      :    print this help message" 
    print >> sys.stderr, ""

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
        print >> sys.stderr, "m is %s" % (m_val)
    else:
        print >> sys.stderr, "Error: -m option not given"
        sys.exit(2)

    if(f_given==True):
        print >> sys.stderr, "s is %s" % (f_val)
    else:
        print >> sys.stderr, "Error: -f option not given"
        sys.exit(2)

    # load mapping file
    mapinfo=load_mapping_info(m_val)

    # process file
    process_file(mapinfo,f_val)
    
if __name__ == "__main__":
    main(sys.argv)
