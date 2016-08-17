# Author: Daniel Ortiz Mart\'inez
# *- bash -*

########
function extract_fluxnum_val()
{
    $AWK -F "\"" \
     '{
       if(substr($2,1,1)=="v") 
       {
        for(i=3;i<=NF;++i)
         if($i==" value=")
         {
           printf"%d %s\n",substr($2,2),$(i+1)
           break
         }
       }
      }'
}

########
function print_fluxnum_rid_val()
{
    $AWK -v rids=${_mfilepref}_reaction_ids.csv \
     'BEGIN{
             while( (getline < rids) > 0)
             {
               commapos=index($0,",")
               rnum=substr($0,1,commapos-1)
               rid=substr($0,commapos+1)
               ridmap[rnum]=rid
             }
      }
      {
        printf"%d,%s,%s\n",$1,ridmap[$1],$2
      }'
}

########
function get_cplex_fluxes_csv()
{
    # Take parameters
    _cplexfile=$1
    _mfilepref=$2

    # Get fluxes
    $GREP "variable name" ${_cplexfile} | \
        extract_fluxnum_val | print_fluxnum_rid_val
}

########
function get_cplex_fluxes_json()
{
    # Take parameters
    _cplexfile=$1
    _mfilepref=$2

    # Obtain number of fluxes
    nfluxes=`$GREP "variable name" ${_cplexfile} | $AWK -F "\"" '{if(substr($2,1,1)=="v") printf"$s\n"}' | wc -l`

    # Get fluxes
    echo "{"
    get_cplex_fluxes_csv ${_cplexfile} ${_mfilepref} | \
        $AWK -F "," -v nl=$nfluxes '{printf"\"%s\": %s",$2,$3; if(NR<nl) printf","; printf"\n"}'
    echo "}"
}

########
if [ $# -lt 1 ]; then
    echo "Use: get_cplex_fluxes -f <string> [-of <int>]"
    echo ""
    echo "-f <string>   : path to file with cplex solution"
    echo "-of <int>     : output format, can be one of following,"    
    echo "                0 -> csv"    
    echo "                1 -> json"    
    echo ""
else
    
    # Read parameters
    f_given=0
    m_given=0
    of_given=0
    oformat=0
    while [ $# -ne 0 ]; do
        case $1 in
        "-f") shift
            if [ $# -ne 0 ]; then
                cplexfile=$1
                f_given=1
            fi
            ;;
        "-m") shift
            if [ $# -ne 0 ]; then
                mfilepref=$1
                m_given=1
            fi
            ;;
        "-of") shift
            if [ $# -ne 0 ]; then
                oformat=$1
                of_given=1
            fi
            ;;
        esac
        shift
    done

    # Check parameters
    if [ ${f_given} -eq 0 ]; then
        echo "Error! -f parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${cplexfile} ]; then
        echo "Error! ${file} file does not exist" >&2
        exit 1
    fi

    if [ ${m_given} -eq 0 ]; then
        echo "Error! -m parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${mfilepref}_reaction_ids.csv ]; then
        echo "Error! ${mfilepref}_reaction_ids.csv file does not exist" >&2
        exit 1
    fi

    ### Process parameters
    
    # Generate basic statistics
    case $oformat in
        0)
            get_cplex_fluxes_csv $cplexfile $mfilepref
            ;;
        1)
            get_cplex_fluxes_json $cplexfile $mfilepref
            ;;
    esac

fi
