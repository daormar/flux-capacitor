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
function solve_fvar_lp_prob()
{
    # Initialize variables
    _fvars=$1
    _fva_templ=$2
    _fba_sol=$3
    _g_val=$4
    _mst=$5
    _rt_val=0.01

    # Process flux variables
    cat ${_fvars} | while read fvar; do
        # Instantiate fva templates
        $bindir/instantiate_fva_templ -f ${_fva_templ} -d 0 -v ${fvar} \
            -s ${_fba_sol} -g ${_g_val} > ${outd}/fvar_lp/${fvar}_min.lp 
        $bindir/instantiate_fva_templ -f ${_fva_templ} -d 1 -v ${fvar} \
            -s ${_fba_sol} -g ${_g_val} > ${outd}/fvar_lp/${fvar}_max.lp 

        # Solve lp problems
        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/fvar_lp/${fvar}_min.lp" "set mip tolerances mipgap ${_rt_val}" \
            "read ${_mst}" "optimize" "write ${outd}/fvar_lp/${fvar}_min.sol" \
             "write ${outd}/fvar_lp/${fvar}_min.mst all" > ${outd}/fvar_lp/${fvar}_min.log 2>&1 || exit 1

        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/fvar_lp/${fvar}_max.lp" "set mip tolerances mipgap ${_rt_val}" \
            "read ${_mst}" "optimize" "write ${outd}/fvar_lp/${fvar}_max.sol" \
            "write ${outd}/fvar_lp/${fvar}_max.mst all" > ${outd}/fvar_lp/${fvar}_max.log 2>&1 || exit 1
    done
}

########
if [ $# -lt 1 ]; then
    echo "Use: auto_fva -p <string> -o <string> [-g <float>] [-rt <float>]"
    echo ""
    echo "-p <string>   : prefix of lp files"
    echo "-o <string>   : output directory"
    echo "-g <float>    : value of the gamma parameter (between 0 and 1, 1 by default)"
    echo "-rt <float>   : relative tolerance gap for initial fba (0.01 by default)"
    echo ""
else
    
    # Read parameters
    p_given=0
    o_given=0
    g_given=0
    g_val=1
    rt_given=0
    rt_val=0.01
    while [ $# -ne 0 ]; do
        case $1 in
        "-p") shift
            if [ $# -ne 0 ]; then
                pref=$1
                p_given=1
            fi
            ;;
        "-o") shift
            if [ $# -ne 0 ]; then
                outd=$1
                o_given=1
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
        esac
        shift
    done

    # Check parameters
    if [ ${p_given} -eq 0 ]; then
        echo "Error! -p parameter not given" >&2
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

    if [ -d ${outd} ]; then
        echo "Warning! ${outd} directory already exists" >&2
    else
        mkdir ${outd} || { echo "Error! cannot create output directory" >&2; exit 1; }
    fi

    ### Print parameters
    if [ ${p_given} -eq 1 ]; then
        echo "-p parameter is ${pfile}" >> ${outd}/params.txt
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
        "write ${outd}/fba/fba.mst all" > ${outd}/fba/cplex.log 2>&1 || exit 1
    $bindir/gen_fba_stats -f ${outd}/fba/fba.sol -c 0 > ${outd}/fba/fba.sol.stats
    fba_sol=`$GREP "Objective value" ${outd}/fba/fba.sol.stats | $AWK '{printf"%s\n",$4}'`
        
    # Generate list of flux variables to be studied
    echo "* Generating list of flux variables to be studied..." >&2
    echo "" >&2
    create_out_dir ${outd}/fvars
    extract_fvars_from_lpf ${fba_file} > ${outd}/fvars/fvars.txt

    # Solve lp problems for flux variables
    echo "* Solving lp problems for flux variables..." >&2
    echo "" >&2
    create_out_dir ${outd}/fvar_lp
    solve_fvar_lp_prob ${outd}/fvars/fvars.txt ${fva_templ} ${fba_sol} \
        ${g_val} ${outd}/fba/fba.mst

fi
