########
if [ $# -lt 1 ]; then
    echo "Use: solve_fva_for_vlist [-pr <string>] -f <string> -t <string> -s <float>"
    echo "                         -o <string> [-g <float>] [-m <string>] [-rt <float>]"
    echo ""
    echo "-pr <int>     : number of processors"
    echo "-f <string>   : file with flux variables to be analyzed"
    echo "-t <string>   : file with fva template"
    echo "-s <float>    : solution of fba problem"
    echo "-o <string>   : output directory"
    echo "-g <float>    : value of the gamma parameter (between 0 and 1, 1 by default)"
    echo "-m <string>   : file with MIP start"
    echo "-rt <float>   : relative tolerance gap provided to lp solver (0.01 by default)"
    echo ""
else
    
    # Read parameters
    pr_given=0
    nprocs=1
    f_given=0
    t_given=0
    s_given=0
    o_given=0
    g_given=0
    g_val=1
    m_given=0
    rt_given=0
    rt_val=0.01
    while [ $# -ne 0 ]; do
        case $1 in
        "-pr") shift
            if [ $# -ne 0 ]; then
                nprocs=$1
                pr_given=1
            fi
            ;;
        "-f") shift
            if [ $# -ne 0 ]; then
                fvars=$1
                f_given=1
            fi
            ;;
        "-t") shift
            if [ $# -ne 0 ]; then
                fva_templ=$1
                t_given=1
            fi
            ;;
        "-s") shift
            if [ $# -ne 0 ]; then
                fba_sol=$1
                s_given=1
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
        "-m") shift
            if [ $# -ne 0 ]; then
                mst=$1
                m_given=1
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
    if [ ${f_given} -eq 0 ]; then
        echo "Error! -f parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${fvars} ]; then
        echo "Error! ${fvars} file does not exist" >&2
        exit 1
    fi

    if [ ${t_given} -eq 0 ]; then
        echo "Error! -t parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${fva_templ} ]; then
        echo "Error! ${fva_templ} file does not exist" >&2
        exit 1
    fi

    if [ ${s_given} -eq 0 ]; then
        echo "Error! -s parameter not given" >&2
        exit 1
    fi

    if [ ${o_given} -eq 0 ]; then
        echo "Error! -o parameter not given" >&2
        exit 1
    fi

    if [ ! -d ${outd} ]; then
        echo "Error: ${outd} directory does not exist" >&2
    fi

    # Print parameters
    if [ ${pr_given} -eq 1 ]; then
        echo "-pr parameter is ${nproc}" >> ${outd}/params.txt
    fi

    if [ ${f_given} -eq 1 ]; then
        echo "-f parameter is ${fvars}" >> ${outd}/params.txt
    fi

    if [ ${t_given} -eq 1 ]; then
        echo "-t parameter is ${fva_templ}" >> ${outd}/params.txt
    fi

    if [ ${s_given} -eq 1 ]; then
        echo "-s parameter is ${fba_sol}" >> ${outd}/params.txt
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o parameter is ${outd}" >> ${outd}/params.txt
    fi

    if [ ${g_given} -eq 1 ]; then
        echo "-g parameter is ${g_val}" >> ${outd}/params.txt
    fi

    if [ ${m_given} -eq 1 ]; then
        echo "-m parameter is ${mst}" >> ${outd}/params.txt
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

    # Process flux variables
    cat ${fvars} | while read fvar; do
        # Instantiate fva templates
        $bindir/instantiate_fva_templ -f ${fva_templ} -d 0 -v ${fvar} \
            -s ${fba_sol} -g ${g_val} > ${outd}/${fvar}_min.lp 
        $bindir/instantiate_fva_templ -f ${fva_templ} -d 1 -v ${fvar} \
            -s ${fba_sol} -g ${g_val} > ${outd}/${fvar}_max.lp 

        # Solve lp problems
        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/${fvar}_min.lp" "set mip tolerances mipgap ${rt_val}" \
            "read ${mst}" "optimize" "write ${outd}/${fvar}_min.sol" \
             "write ${outd}/${fvar}_min.mst all" > ${outd}/${fvar}_min.log 2>&1 || exit 1

        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/${fvar}_max.lp" "set mip tolerances mipgap ${rt_val}" \
            "read ${mst}" "optimize" "write ${outd}/${fvar}_max.sol" \
            "write ${outd}/${fvar}_max.mst all" > ${outd}/${fvar}_max.log 2>&1 || exit 1
    done

fi
