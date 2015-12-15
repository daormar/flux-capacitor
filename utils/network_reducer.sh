# Author: Daniel Ortiz Mart\'inez
# *- bash -*

########
function create_out_dir()
{
    _dir=$1
    if [ ! -d ${_dir} ]; then
        mkdir ${_dir} || { echo "Error while creating output directories: ${_dir}" >&2; exit 1; }
    fi
}

########
function netred()
{
    echo "TBD"
}

########
if [ $# -lt 1 ]; then
    echo "Use: network_reducer -m <string> -pm <string> -pr <string>"
    echo "                     -md <int> -mr <int> -o <string>"
    echo "                     [-qs <string>] [-sdir <string>]"
    echo ""
    echo "-m <string>    : path to file with metabolic model in SBML format"
    echo "-pm <string>   : file with list of protected metabolites"
    echo "-pr <string>   : file with list of protected reactions"
    echo "-md <int>      : minimum degrees of freedom"
    echo "-mr <int>      : minimum number of reactions"
    echo "-o <string>    : output directory"
    echo "-qs <string>   : specific options to be given to the qsub command"
    echo "                 (example: -qs \"-l pmem=1gb\")."
    echo "-sdir <string> : absolute path of a directory common to all"
    echo "                 processors. If not given, \$HOME will be used."
    echo "                 NOTES:"
    echo "                  a) give absolute paths when using pbs clusters"
    echo "                  b) ensure there is enough disk space in the partition"
    echo ""
else
    
    # Read parameters
    m_given=0
    pm_given=0
    pr_given=0
    md_given=0
    mr_given=0
    o_given=0
    sdir=$HOME
    while [ $# -ne 0 ]; do
        case $1 in
        "-m") shift
            if [ $# -ne 0 ]; then
                mfile=$1
                m_given=1
            fi
            ;;
        "-pm") shift
            if [ $# -ne 0 ]; then
                pmfile=$1
                pm_given=1
            fi
            ;;
        "-r") shift
            if [ $# -ne 0 ]; then
                prfile=$1
                pr_given=1
            fi
            ;;
        "-md") shift
            if [ $# -ne 0 ]; then
                mdval=$1
                md_given=1
            fi
            ;;
        "-mr") shift
            if [ $# -ne 0 ]; then
                mrval=$1
                mr_given=1
            fi
            ;;
        "-o") shift
            if [ $# -ne 0 ]; then
                outd=$1
                o_given=1
            fi
            ;;
        "-qs") shift
            if [ $# -ne 0 ]; then
                qs_opt="-qs"
                qs_par=$1
                qs_given=1
            else
                qs_given=0
            fi
            ;;
        "-sdir") shift
            if [ $# -ne 0 ]; then
                sdir=$1                
            fi
            ;;
        esac
        shift
    done

    # Check parameters
    if [ ${m_given} -eq 0 ]; then
        echo "Error! -m parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${mfile} ]; then
        echo "Error! ${mfile} file does not exist" >&2
        exit 1
    fi

    if [ ${pm_given} -eq 0 ]; then
        echo "Error! -pm parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${pmfile} ]; then
        echo "Error! ${pmfile} file does not exist" >&2
        exit 1
    fi

    if [ ${pr_given} -eq 0 ]; then
        echo "Error! -pr parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${prfile} ]; then
        echo "Error! ${prfile} file does not exist" >&2
        exit 1
    fi

    if [ ${md_given} -eq 0 ]; then
        echo "Error! -md parameter not given" >&2
        exit 1
    fi

    if [ ${mr_given} -eq 0 ]; then
        echo "Error! -mr parameter not given" >&2
        exit 1
    fi

    if [ ${o_given} -eq 0 ]; then
        echo "Error! -o parameter not given" >&2
        exit 1
    fi

    if [ -d ${outd} ]; then
        echo "Warning! ${outd} directory already exists" >&2
    else
        mkdir ${outd} || { echo "Error! cannot create output directory" >&2; exit 1; }
    fi

    ### Print parameters
    if [ ${m_given} -eq 1 ]; then
        echo "-m parameter is ${mfile}" > ${outd}/params.txt
    fi

    if [ ${pm_given} -eq 1 ]; then
        echo "-pm parameter is ${pmfile}" >> ${outd}/params.txt
    fi

    if [ ${pr_given} -eq 1 ]; then
        echo "-pr parameter is ${prfile}" >> ${outd}/params.txt
    fi

    if [ ${md_given} -eq 1 ]; then
        echo "-md parameter is ${mdval}" >> ${outd}/params.txt
    fi

    if [ ${mr_given} -eq 1 ]; then
        echo "-mr parameter is ${mrval}" >> ${outd}/params.txt
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o parameter is ${outd}" >> ${outd}/params.txt
    fi

    ### Process parameters

    # create shared directory
    SDIR="${sdir}/network_reducer_$$"
    mkdir $SDIR || { echo "Error: shared directory cannot be created"  >&2 ; exit 1; }

    # remove temp directories on exit
    if [ $debug -eq 0 ]; then
        trap "rm -rf $SDIR 2>/dev/null" EXIT
    fi

    # Create output directory
    create_out_dir ${outd}/fba

    # Execute network reducer algorithm
    netred || exit 1

fi
