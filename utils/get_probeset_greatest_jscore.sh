# Author: Daniel Ortiz Mart\'inez
# *- bash -*

########
if [ $# -lt 1 ]; then
    echo "Use: get_probeset_greatest_jscore -f <string>"
    echo ""
    echo "-f <string>   :  file with jetset scores"
    echo ""
else
    
    # Read parameters
    f_given=0
    while [ $# -ne 0 ]; do
        case $1 in
        "-f") shift
            if [ $# -ne 0 ]; then
                jsfile=$1
                f_given=1
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

    if [ ! -f ${jsfile} ]; then
        echo "Error! ${jsfile} file does not exist" >&2
        exit 1
    fi

    # Process parameters
    cat ${jsfile} | $AWK -F ',' '{if($9!="NA") printf"%s\n",$0}' | LC_ALL=C $SORT -t ',' -k9,9 -k8,8gr |\
        $SED 's/,/ /g' | $UNIQ -f 8 | $SED 's/ /,/g'
fi
