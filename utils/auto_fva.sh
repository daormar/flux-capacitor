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
function extract_fvars_from_lpf()
{
    # Initialize variables
    _fba_file=$1
    
    # Extract flux variables
    cat ${_fba_file} | $AWK '{for(i=1;i<=NF;++i) if(match($i,"v")==1) printf"%s\n",$i}' | $SORT | $UNIQ
}

########
if [ $# -lt 1 ]; then
    echo "Use: auto_fva [-pr <int>] -l <string> -o <string> [-v <string>] [-g <float>]"
    echo "              [-rt <float>] [-qs <string>] [-sdir <string>]"
    echo ""
    echo "-pr <int>      : number of processors (1 by default)"
    echo "-l <string>    : prefix of lp files"
    echo "-o <string>    : output directory"
    echo "-v <string>    : file with variables to be analyzed (if not given, the whole"
    echo "                 set of variables contained in the lp file are analyzed)"
    echo "-g <float>     : value of the gamma parameter (between 0 and 1, 0.9 by default)"
    echo "-rt <float>    : relative tolerance gap (0.01 by default)"
    echo "-qs <string>   : specific options to be given to the qsub command"
    echo "                 (example: -qs \"-l pmem=1gb\")"
    echo "-sdir <string> : absolute path of a directory common to all"
    echo "                 processors. If not given, \$HOME will be used."
    echo "                 NOTES:"
    echo "                  a) give absolute paths when using pbs clusters"
    echo "                  b) ensure there is enough disk space in the partition"
    echo ""
else
    
    # Read parameters
    pr_given=0
    nprocs=1
    l_given=0
    o_given=0
    v_given=0
    g_given=0
    g_val=0.9
    rt_given=0
    rt_val=0.01
    sdir=$HOME
    while [ $# -ne 0 ]; do
        case $1 in
        "-pr") shift
            if [ $# -ne 0 ]; then
                nprocs=$1
                pr_given=1
            fi
            ;;
        "-l") shift
            if [ $# -ne 0 ]; then
                pref=$1
                l_given=1
            fi
            ;;
        "-o") shift
            if [ $# -ne 0 ]; then
                outd=$1
                o_given=1
            fi
            ;;
        "-v") shift
            if [ $# -ne 0 ]; then
                v_val=$1
                v_given=1
            fi
            ;;
        "-g") shift
            if [ $# -ne 0 ]; then
                g_val=$1
                g_given=1
            fi
            ;;
        "-rt") shift
            if [ $# -ne 0 ]; then
                rt_val=$1
                rt_given=1
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
    if [ ${l_given} -eq 0 ]; then
        echo "Error! -l parameter not given" >&2
        exit 1
    fi

    fba_file=${pref}.lp
    fva_templ=${pref}_fva_template.lp

    if [ ! -f ${fba_file} ]; then
        echo "Error! ${fba_file} file does not exist" >&2
        exit 1
    fi

    if [ ! -f ${fva_templ} ]; then
        echo "Error! ${fva_templ} file does not exist" >&2
        exit 1
    fi

    if [ ${o_given} -eq 0 ]; then
        echo "Error! -o parameter not given" >&2
        exit 1
    fi

    if [ ${v_given} -eq 1 ]; then
        if [ ! -f ${v_val} ]; then
            echo "Error! ${v_val} file does not exist" >&2
            exit 1
        fi
    fi

    if [ -d ${outd} ]; then
        echo "Warning! ${outd} directory already exists" >&2
    else
        mkdir ${outd} || { echo "Error! cannot create output directory" >&2; exit 1; }
    fi

    if [ ! -d ${sdir} ]; then
        echo "Error! ${sdir} directory does not exist" >&2
        exit 1
    fi

    ### Print parameters
    if [ ${l_given} -eq 1 ]; then
        echo "-l parameter is ${pref}" > ${outd}/params.txt
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o parameter is ${outd}" >> ${outd}/params.txt
    fi

    if [ ${g_given} -eq 1 ]; then
        echo "-g parameter is ${g_val}" >> ${outd}/params.txt
    fi

    if [ ${rt_given} -eq 1 ]; then
        echo "-rt parameter is ${rt_val}" >> ${outd}/params.txt
    fi

    # check presence of cplex
    if [ ! -f ${CPLEX_BINARY_DIR}/cplex ]; then
        echo "Error, CPLEX binary not found (shell variable CPLEX_BINARY_DIR should be defined)">&2
        exit 1
    fi

    ### Process parameters

    # Copying initial lp files
    echo "* Copying initial lp files..." >&2
    echo "" >&2
    create_out_dir ${outd}/initial_lp
    cp ${fba_file} ${outd}/initial_lp
    cp ${fva_templ} ${outd}/initial_lp

    # Generate solution and MIP start file
    echo "* Generating initial solution and MIP start file..." >&2
    echo "" >&2
    create_out_dir ${outd}/fba
    ${CPLEX_BINARY_DIR}/cplex -c "read ${fba_file}" "set mip tolerances mipgap ${rt_val}" \
        "optimize" "write ${outd}/fba/fba.sol" \
        "write ${outd}/fba/fba.mst all" > ${outd}/fba/cplex.log || exit 1
    $bindir/gen_fba_stats -f ${outd}/fba/fba.sol -c 0 > ${outd}/fba/fba.sol.stats
    fba_sol=`$GREP "Objective value" ${outd}/fba/fba.sol.stats | $AWK '{printf"%d\n",int($4)}'`

    # Define -m option for solve_fva_for_vlist program
    if [ -f ${outd}/fba/fba.mst ]; then
        m_opt="-m ${outd}/fba/fba.mst"
    else
        m_opt=""
    fi

    # Generate list of flux variables to be studied
    create_out_dir ${outd}/fvars

    if [ ${v_given} ]; then
        cp ${v_val} ${outd}/fvars/fvars.txt
    else
        echo "* Generating list of flux variables to be studied..." >&2
        echo "" >&2
        extract_fvars_from_lpf ${fba_file} > ${outd}/fvars/fvars.txt
    fi

    # Solve lp problems for flux variables
    echo "* Solving lp problems for flux variables..." >&2
    echo "" >&2
    create_out_dir ${outd}/fvar_lp
    $bindir/solve_fva_for_vlist -pr ${nprocs} -f ${outd}/fvars/fvars.txt \
        -t ${fva_templ} -s ${fba_sol} -g ${g_val} ${m_opt} \
        -rt ${rt_val} -o ${outd}/fvar_lp ${qs_opt} "${qs_par}" -sdir $sdir || exit 1

fi
