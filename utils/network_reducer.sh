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
function obtain_removable_reac()
{
    cat $SDIR/reacs $prfile | $SORT -n | $UNIQ -u > $SDIR/removable_reacs_aux
    cat $SDIR/reacs $SDIR/removable_reacs_aux | $SORT -n | $UNIQ -d
}

########
function obtain_nremreac()
{
    _nremreac=`wc -l $SDIR/removable_reacs | $AWK '{printf"%s\n",$1}'`
    echo ${_nremreac}
}

########
function obtain_reacs()
{
    tail -n +3 $SDIR/curr_minfo/model_sparse_st_matrix.csv | $AWK '{printf"%s\n",$2}' | $SORT -n | $UNIQ
}

########
function obtain_nreac()
{
    _nreac=`wc -l $SDIR/reacs | $AWK '{printf"%s\n",$1}'`
    echo ${_nreac}
}

########
function obtain_matrix_rank()
{
    _rank=`$bindir/calc_matrix_rank -m $SDIR/curr_minfo/model_sparse_st_matrix.csv 2>/dev/null | $AWK '{printf"%s\n",$2}'`
    echo ${_rank}
}

########
function netred()
{
    # Initialize variables
    niter=1
    maxiters=100
    end=0

    # Create temporary directory with metabolic model information
    cp -r ${auto_fba_outdir}/minfo $SDIR/curr_minfo
    
    # Execute network reduction loop
    echo "Executing network reduction loop...">&2

    while [ $end -eq 0 ]; do
        ## Check ending conditions
        
        # Check maximum number of iterations
        if [ $niter -gt $maxiters ]; then
            echo "Warning: maximum number of iterations exceeded!" >&2
            end=1
        fi

        # Obtain number of reactions
        obtain_reacs > $SDIR/reacs
        nreac=`obtain_nreac`
        if [ $nreac -lt $mrval ]; then
            echo "Number of reactions ($nreac) are lower than minimum ($mrval)" >&2
            end=1
        fi

        # Check degrees of freedom
        rank=`obtain_matrix_rank`
        dof=`expr $nreac - $rank`
        if [ $dof -lt $mdval ]; then
            echo "Degrees of freedom ($dof) are lower than minimum ($mdval)" >&2
            end=1
        fi

        # Obtain set of removable reactions
        obtain_removable_reac > $SDIR/removable_reacs
        nremreac=`obtain_nremreac`
        if [ $nremreac -eq 0 ]; then
            echo "There are not any removable reactions" >&2
            end=1
        fi

        # Print iteration information
        echo "* Iteration $niter (nreac: $nreac ; nremreac: $nremreac ; dof: $dof)" >&2

        # Prune network if ending conditions were not met
        if [ $end -eq 0 ]; then

            # Create lp file
            $bindir/create_lp_file -s $SDIR/curr_minfo/model \
                -a ${auto_fba_outdir}/abs_pres_info/abs_pres_genes_filt_${cel_file}.csv \
                -m ${auto_fba_outdir}/esetdir/esetgenes_to_entrezids.csv \
                -c 0 > $SDIR/curr_fba_problem.lp 2> $SDIR/create_lp_file.log

            # Execute fva
            # TBD

            # Internal while loop
            # TBD

            # Increase number of iterations
            niter=`expr $niter + 1`
        fi

    done
}

########
if [ $# -lt 1 ]; then
    echo "Use: network_reducer -a <string> -c <string> -pm <string> -pr <string>"
    echo "                     -md <int> -mr <int> -o <string> [-np <int>]"
    echo "                     [-qs <string>] [-sdir <string>] [-debug]"
    echo ""
    echo "-a <string>    : directory storing the output of auto_fba tool"
    echo "-c <string>    : CEL file name chosen from those analyzed with auto_fba"
    echo "-pm <string>   : file with list of protected metabolites"
    echo "-pr <string>   : file with list of protected reactions"
    echo "-md <int>      : minimum degrees of freedom"
    echo "-mr <int>      : minimum number of reactions"
    echo "-o <string>    : output directory"
    echo "-np <int>      : number of processors"
    echo "-qs <string>   : specific options to be given to the qsub command"
    echo "                 (example: -qs \"-l pmem=1gb\")."
    echo "-sdir <string> : absolute path of a directory common to all"
    echo "                 processors. If not given, \$HOME will be used."
    echo "                 NOTES:"
    echo "                  a) give absolute paths when using pbs clusters"
    echo "                  b) ensure there is enough disk space in the partition"
    echo "-debug         : After ending, do not delete temporary files"
    echo "                 (for debugging purposes)"
    echo ""
else
    
    # Read parameters
    a_given=0
    c_given=0
    pm_given=0
    pr_given=0
    md_given=0
    mr_given=0
    np_given=0
    nprocs=1
    o_given=0
    sdir=$HOME
    debug=0
    while [ $# -ne 0 ]; do
        case $1 in
        "-a") shift
            if [ $# -ne 0 ]; then
                auto_fba_outdir=$1
                a_given=1
            fi
            ;;
        "-c") shift
            if [ $# -ne 0 ]; then
                cel_file=$1
                c_given=1
            fi
            ;;
        "-pm") shift
            if [ $# -ne 0 ]; then
                pmfile=$1
                pm_given=1
            fi
            ;;
        "-pr") shift
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
        "-np") shift
            if [ $# -ne 0 ]; then
                nprocs=$1
                np_given=1
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
        "-debug") debug=1
            ;;
        esac
        shift
    done

    # Check parameters
    if [ ${a_given} -eq 0 ]; then
        echo "Error! -a parameter not given" >&2
        exit 1
    fi

    if [ ! -d ${auto_fba_outdir} ]; then
        echo "Error! ${auto_fba_outdir} directory does not exist" >&2
        exit 1
    fi

    if [ ${c_given} -eq 0 ]; then
        echo "Error! -c parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${auto_fba_outdir}/abs_pres_info/abs_pres_genes_filt_${cel_file}.csv ]; then
        echo "Error! file with absent/present genes for CEL file ${cel_file} does not exist" >&2
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
    if [ ${a_given} -eq 1 ]; then
        echo "-a parameter is ${auto_fba_outdir}" > ${outd}/params.txt
    fi

    if [ ${c_given} -eq 1 ]; then
        echo "-c parameter is ${cel_file}" > ${outd}/params.txt
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

    if [ ${np_given} -eq 1 ]; then
        echo "-np parameter is ${nprocs}" >> ${outd}/params.txt
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
    create_out_dir ${outd}

    # Execute network reducer algorithm
    netred || exit 1

fi
