# Author: Daniel Ortiz Mart\'inez
# *- bash -*

########
function filter_file()
{
    ff_par=$1
    cat ${ff_par} | $AWK -F "," -v gfile=$gfile 'BEGIN{
                                               while( (getline <gfile) > 0)
                                               {
                                                filt_eids[$1]=$3
                                               }
                                              }
                                              {
                                               if(FNR==1)
                                                 printf"%s\n",$0
                                               else
                                               {
                                                if($1 in filt_eids)
                                                  printf"%s\n",$0                                               
                                               }
                                              }'
}

########
if [ $# -lt 1 ]; then
    echo "Use: filter_panp_results -p <string> -g <string>"
    echo ""
    echo "-p <string>   :  file with panp results (generated using exec_panp_eset)"
    echo "-g <string>   :  file with greatest j-scores for probesets"
    echo ""
else
    
    # Read parameters
    p_given=0
    g_given=0
    while [ $# -ne 0 ]; do
        case $1 in
        "-p") shift
            if [ $# -ne 0 ]; then
                pfile=$1
                p_given=1
            fi
            ;;
        "-g") shift
            if [ $# -ne 0 ]; then
                gfile=$1
                g_given=1
            fi
            ;;
        esac
        shift
    done

    # Check parameters
    if [ ${p_given} -eq 0 ]; then
        echo "Error! -p parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${pfile} ]; then
        echo "Error! ${pfile} file does not exist" >&2
        exit 1
    fi

    if [ ${g_given} -eq 0 ]; then
        echo "Error! -g parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${gfile} ]; then
        echo "Error! ${gfile} file does not exist" >&2
        exit 1
    fi

    # Process parameters
    filter_file $pfile

fi
