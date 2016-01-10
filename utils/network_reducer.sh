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
    cat $SDIR/reacs $lprfile | $SORT -n | $UNIQ -u > $SDIR/removable_reacs_aux
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
function obtain_fba_criterion()
{
    _auto_fba_outdir=$1
    grep "\-c parameter is" ${_auto_fba_outdir}/params.txt | $AWK '{printf"%s\n",$4}'
}

########
function shlomi_fva()
{
    $bindir/create_lp_file -s $SDIR/curr_minfo/model \
        -a ${auto_fba_outdir}/abs_pres_info/abs_pres_genes_filt_${cel_file}.csv \
        -m ${auto_fba_outdir}/esetdir/esetgenes_to_entrezids.csv \
        -c 1 > $SDIR/lp/${cel_file}.lp 2> $SDIR/lp/${cel_file}.log

    # Generate template for fva analysis in lp format
    $bindir/create_lp_file -s ${SDIR}/curr_minfo/model \
        -a ${auto_fba_outdir}/abs_pres_info/abs_pres_genes_filt_${cel_file}.csv \
        -m ${auto_fba_outdir}/esetdir/esetgenes_to_entrezids.csv \
        -c 1 --fva > $SDIR/lp/${cel_file}_fva_template.lp \
        2> $SDIR/lp/${cel_file}_fva_template.log || exit 1
    
    echo "- Executing fva..." >&2
                    
    # Execute fva
    $bindir/auto_fva -l $SDIR/lp/${cel_file} -o $SDIR/fva \
        -g ${g_val} -rt ${rt_val} ${qs_opt} "${qs_par}" -sdir ${sdir} 2> $SDIR/fva.log
}

########
function obtain_cand_reac_for_removal()
{
    # Obtain flux differences in ascending order
    $AWK '{printf"%s %g\n",$1,$15}' $SDIR/fva/fvar_lp/results | $SORT -gk2  > $SDIR/sorted_flux_diff

    # Obtain first removable reaction
    cat $SDIR/sorted_flux_diff | while read fluxdiff; do
        _reac=`echo $fluxdiff | $AWK '{printf"%d\n",substr($1,2)}'`
        _not_remov=`$GREP ${_reac} $SDIR/removable_reacs | wc -l`
        if [ ${_not_remov} -eq 0 ]; then
            break
        fi
    done

    # Return reaction
    echo ${_reac}
}

########
function remove_reac_from_remov()
{
    # Take parameters
    _reac=$1
    
    # Remove reaction from $SDIR/removable_reacs
    $AWK -v reac=${_reac} '{if($1!=reac) printf"%s\n",$1}' $SDIR/removable_reacs > $SDIR/removable_reacs_aux
    cp $SDIR/removable_reacs_aux $SDIR/removable_reacs
}

########
remove_reac_from_nw()
{
    # TBD

    # Take parameters
    _reac=$1
    currmi_dir=$2
    currmi_aux_dir=$3
    
    # Copy current model to auxiliary model
    cp ${currmi_dir}/* ${currmi_dir_aux}

    # 
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
            echo "Number of reactions ($nreac) is lower than minimum ($mrval)" >&2
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

            echo "- Creating lp files..." >&2

            # Create lp file
            case $crit in
                1)
                    shlomi_fva
                    ;;
            esac

            # Internal while loop
            success=1

            while [ $success -eq 1 -a $nremreac -ne 0 ]; do

                # Obtain candidate reaction for removal
                reac=`obtain_cand_reac_for_removal`

                # Remove reaction from set of removables
                remove_reac_from_remov $reac

                # Obtain number of removable reactions
                nremreac=`obtain_nremreac`

                # Remove reaction from network and store it in an
                # auxiliary directory
                remove_reac_from_nw $reac $SDIR/curr_minfo $SDIR/curr_minfo_aux

                # Check protected functions
#                success=`check_protected_functions`

                # Check success
                if [ $success -eq 1 ]; then
                    # Replace current network with auxiliary network
                    rm -rf $SDIR/curr_minfo
                    cp -r $SDIR/curr_minfo_aux $SDIR/curr_minfo
                fi
                
                # Clear curr_minfo_aux
                rm $SDIR/curr_minfo_aux/*

            done

            # Increase number of iterations
            niter=`expr $niter + 1`
        fi

    done
}

########
if [ $# -lt 1 ]; then
    echo "Use: network_reducer [-pr <int>] -a <string> -lpm <string> -lpr <string>"
    echo "                     -md <int> -mr <int> -o <string>"
    echo "                     [-cf <string>] [-g <float>] [-rt <float>]"
    echo "                     [-qs <string>] [-sdir <string>] [-debug]"
    echo ""
    echo "-pr <int>      : number of processors"
    echo "-a <string>    : directory storing the output of auto_fba tool"
    echo "-lpm <string>  : file with list of protected metabolites"
    echo "-lpr <string>  : file with list of protected reactions"
    echo "-md <int>      : minimum degrees of freedom"
    echo "-mr <int>      : minimum number of reactions"
    echo "-o <string>    : output directory"
    echo "-cf <string>   : CEL file name chosen from those analyzed with auto_fba"
    echo "                 (required if auto_fba was executed with -c 1 option)"
    echo "-g <float>     : value of the gamma parameter (between 0 and 1, 1 by default)"
    echo "-rt <float>    : relative tolerance gap (0.01 by default)"
    echo "-qs <string>   : specific options to be given to the qsub command"
    echo "                 (example: -qs \"-l pmem=1gb\")."
    echo "-sdir <string> : absolute path of a directory common to all"
    echo "                 processors. If not given, \$HOME will be used."
    echo "                 NOTES:"
    echo "                  a) give absolute paths when using pbs clusters"
    echo "                  b) ensure there is enough disk space in the partition"
    echo "-debug         : after ending, do not delete temporary files"
    echo "                 (for debugging purposes)"
    echo ""
else
    
    # Read parameters
    a_given=0
    lpm_given=0
    lpr_given=0
    md_given=0
    mr_given=0
    pr_given=0
    nprocs=1
    o_given=0
    c_given=0
    crit=0
    cf_given=0
    g_given=0
    g_val=1
    rt_given=0
    rt_val=0.01
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
        "-lpm") shift
            if [ $# -ne 0 ]; then
                lpmfile=$1
                lpm_given=1
            fi
            ;;
        "-lpr") shift
            if [ $# -ne 0 ]; then
                lprfile=$1
                lpr_given=1
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
        "-cf") shift
            if [ $# -ne 0 ]; then
                cel_file=$1
                cf_given=1
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
        "-pr") shift
            if [ $# -ne 0 ]; then
                nprocs=$1
                pr_given=1
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

    if [ ! -f ${auto_fba_outdir}/abs_pres_info/abs_pres_genes_filt_${cel_file}.csv ]; then
        echo "Error! file with absent/present genes for CEL file ${cel_file} does not exist" >&2
        exit 1
    fi

    if [ ${lpm_given} -eq 0 ]; then
        echo "Error! -lpm parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${lpmfile} ]; then
        echo "Error! ${lpmfile} file does not exist" >&2
        exit 1
    fi

    if [ ${lpr_given} -eq 0 ]; then
        echo "Error! -lpr parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${lprfile} ]; then
        echo "Error! ${lprfile} file does not exist" >&2
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

    if [ ${lpm_given} -eq 1 ]; then
        echo "-lpm parameter is ${lpmfile}" >> ${outd}/params.txt
    fi

    if [ ${lpr_given} -eq 1 ]; then
        echo "-lpr parameter is ${lprfile}" >> ${outd}/params.txt
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

    if [ ${cf_given} -eq 1 ]; then
        echo "-cf parameter is ${cel_file}" >> ${outd}/params.txt
    fi

    if [ ${g_given} -eq 1 ]; then
        echo "-g parameter is ${g_val}" >> ${outd}/params.txt
    fi

    if [ ${rt_given} -eq 1 ]; then
        echo "-rt parameter is ${rt_val}" >> ${outd}/params.txt
    fi

    if [ ${pr_given} -eq 1 ]; then
        echo "-pr parameter is ${nprocs}" >> ${outd}/params.txt
    fi

    # check presence of cplex
    if [ ! -f ${CPLEX_BINARY_DIR}/cplex ]; then
        echo "Error, CPLEX binary not found (shell variable CPLEX_BINARY_DIR should be defined)">&2
        exit 1
    fi

    ### Process parameters

    # create shared directory
    SDIR="${sdir}/network_reducer_$$"
    mkdir $SDIR || { echo "Error: shared directory cannot be created"  >&2 ; exit 1; }
 
    # create shared subdirectories
    mkdir $SDIR/lp
    mkdir $SDIR/curr_minfo_aux

    # remove temp directories on exit
    if [ $debug -eq 0 ]; then
        trap "rm -rf $SDIR 2>/dev/null" EXIT
    fi

    # Create output directory
    create_out_dir ${outd}

    # Obtain fba criterion
    crit=`obtain_fba_criterion ${auto_fba_outdir}`

    # Execute network reducer algorithm
    netred || exit 1

fi
