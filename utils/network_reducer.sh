# Author: Daniel Ortiz Mart\'inez
# *- bash -*

########
create_out_dir()
{
    _dir=$1
    if [ ! -d ${_dir} ]; then
        mkdir ${_dir} || { echo "Error while creating output directories: ${_dir}" >&2; exit 1; }
    fi
}

########
obtain_removable_reac()
{
    obtain_reacs > $SDIR/reacs
    cat $SDIR/reacs $lprfile | LC_ALL=C $SORT -n | $UNIQ -u > $SDIR/removable_reacs_aux
    cat $SDIR/reacs $SDIR/removable_reacs_aux | LC_ALL=C $SORT -n | $UNIQ -d
}

########
obtain_nremreac()
{
    _nremreac=`wc -l $SDIR/removable_reacs | $AWK '{printf"%s\n",$1}'`
    echo ${_nremreac}
}

########
obtain_reacs()
{
    tail -n +3 $SDIR/curr_minfo/model_sparse_st_matrix.csv | $AWK '{printf"%s\n",$2}' | LC_ALL=C $SORT -n | $UNIQ
}

########
obtain_nreac()
{
    _nreac=`wc -l $SDIR/reacs | $AWK '{printf"%s\n",$1}'`
    echo ${_nreac}
}

########
obtain_matrix_rank()
{
    _rank=`$bindir/calc_matrix_rank -m $SDIR/curr_minfo/model_sparse_st_matrix.csv 2>/dev/null | $AWK '{printf"%s\n",$2}'`
    echo ${_rank}
}

########
obtain_fba_criterion()
{
    _auto_fba_outdir=$1
    grep "\-c parameter is" ${_auto_fba_outdir}/params.txt | $AWK '{printf"%s\n",$4}'
}

########
extract_fvars_from_lpf()
{
    # Initialize variables
    _fba_file=$1
    
    # Extract flux variables
    cat ${_fba_file} | $AWK '{for(i=1;i<=NF;++i) if(match($i,"v")==1) printf"%s\n",$i}' | LC_ALL=C $SORT | $UNIQ
}

########
obtain_flux_ranges_file()
{
    # Initialize variables
    _fva_result_file=$1

    # Obtain flux ranges file
    $AWK '{printf"%d %s\n",substr($1,2),$NF}' ${_fva_result_file} > $SDIR/flux_ranges
}

########
obtain_fva_vars()
{
    # Initialize variables
    _reacs_file=$1

    # Obtain fva variables
    if [ ${li_given} -eq 0 ]; then
        $AWK '{printf "v%06d\n",$1}' ${_reacs_file}
    else
        _rnd_num=$RANDOM
        $AWK '{printf "v%06d\n",$1}' ${_reacs_file} | \
            $bindir/shuffle ${_rnd_num} | $HEAD -${li_val}
    fi
}

########
biomass_fva()
{
    echo "- Creating lp files..." >&2

    $bindir/create_lp_file -s $SDIR/curr_minfo/model \
        -c 0 > $SDIR/lp/biomass.lp 2> $SDIR/lp/biomass.log || exit 1

    # Generate template for fva analysis in lp format
    $bindir/create_lp_file -s ${SDIR}/curr_minfo/model \
        -c 0 --fva > $SDIR/lp/biomass_fva_template.lp \
        2> $SDIR/lp/biomass_fva_template.log || exit 1

    # Obtain file with variables to be analyzed
    obtain_fva_vars $SDIR/removable_reacs > $SDIR/fva_vars

    # Remove $SDIR/fva if exists
    if [ -d $SDIR/fva ]; then
        rm -rf $SDIR/fva
    fi
    
    # Execute fva
    echo "- Executing fva..." >&2
    $bindir/auto_fva -l $SDIR/lp/biomass -o $SDIR/fva -v $SDIR/fva_vars \
                     -g ${g_val} -rt ${rt_val} ${qs_opt} "${qs_par}" -sdir ${sdir} 2> $SDIR/fva.log || exit 1

    # Create file with flux ranges for variable numbers
    obtain_flux_ranges_file $SDIR/fva/fvar_lp/results > $SDIR/flux_ranges
}

########
shlomi_fva()
{
    echo "- Creating lp files..." >&2

    $bindir/create_lp_file -s $SDIR/curr_minfo/model \
        -a ${auto_fba_outdir}/abs_pres_info/abs_pres_genes_${sample_file}.csv \
        -c 1 > $SDIR/lp/${sample_file}.lp 2> $SDIR/lp/${sample_file}.log || exit 1

    # Generate template for fva analysis in lp format
    $bindir/create_lp_file -s ${SDIR}/curr_minfo/model \
        -a ${auto_fba_outdir}/abs_pres_info/abs_pres_genes_${sample_file}.csv \
        -c 1 --fva > $SDIR/lp/${sample_file}_fva_template.lp \
        2> $SDIR/lp/${sample_file}_fva_template.log || exit 1

    # Obtain file with variables to be analyzed
    obtain_fva_vars $SDIR/removable_reacs > $SDIR/fva_vars

    # Remove $SDIR/fva if exists
    if [ -d $SDIR/fva ]; then
        rm -rf $SDIR/fva
    fi

    # Execute fva
    echo "- Executing fva..." >&2
    $bindir/auto_fva -l $SDIR/lp/${sample_file} -o $SDIR/fva -v $SDIR/fva_vars \
                     -g ${g_val} -rt ${rt_val} ${qs_opt} "${qs_par}" -sdir ${sdir} 2> $SDIR/fva.log || exit 1

    # Create file with flux ranges for variable numbers
    obtain_flux_ranges_file $SDIR/fva/fvar_lp/results > $SDIR/flux_ranges
}

########
obtain_sorted_flux_ranges()
{
    echo "- Sorting flux ranges..." >&2

    case ${sort_crit} in
        0)
            # Obtain flux differences in ascending order
            LC_ALL=C $SORT -gk2  $SDIR/flux_ranges > $SDIR/sorted_flux_ranges
            ;;
        1)
            # Obtain flux differences in descending order
            LC_ALL=C $SORT -rgk2  $SDIR/flux_ranges > $SDIR/sorted_flux_ranges
            ;;
        *)
            # Obtain flux differences in ascending order
            LC_ALL=C $SORT -gk2  $SDIR/flux_ranges > $SDIR/sorted_flux_ranges
            ;;
    esac
}

########
obtain_cand_reac_for_removal()
{
    # Check if sorted flux ranges file is empty
    lenfile=`wc -l $SDIR/sorted_flux_ranges | $AWK '{printf"%s",$1}'`
    if [ ${lenfile} -eq 0 ]; then
        echo "NONE"
    else
        # Obtain first reaction contained in sorted flux ranges file
        _reac=`$AWK '{if(NR==1) printf"%s\n",$1}' $SDIR/sorted_flux_ranges`
        
        # Remove first reaction from sorted flux ranges file
        $TAIL -n +2 $SDIR/sorted_flux_ranges > $SDIR/tmp
        mv $SDIR/tmp $SDIR/sorted_flux_ranges
        
        echo ${_reac}
    fi
}

########
remove_reac_from_remov()
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
    # Take parameters
    _reac=$1
    _currmi_dir=$2
    _currmi_aux_dir=$3
    
    # Copy current model to auxiliary model
    cp ${_currmi_dir}/* ${_currmi_aux_dir}

    # Filter reaction from model files
    for f in model_gpr_rules.csv model_reaction_ids.csv model_reaction_lowbnds.csv model_reaction_lowbnds.csv; do
        $AWK -F "," -v reac=${_reac} '{if($1!=reac) printf"%s\n",$0}' ${_currmi_aux_dir}/$f > $SDIR/tmp
        mv $SDIR/tmp ${_currmi_aux_dir}/$f
    done

    # Filter reaction from stoichiometric matrix file
    $AWK -v reac=${_reac} '{if($2!=reac) printf"%s\n",$0}' ${_currmi_aux_dir}/model_sparse_st_matrix.csv > $SDIR/tmp
    mv $SDIR/tmp ${_currmi_aux_dir}/model_sparse_st_matrix.csv
}

########
check_feasibility()
{
    # Take parameters
    _modelinfo_dir=$1

    # Check feasibility of protected reactions
    # TBD

    # Check feasibility of protected metabolites
    feasibility_met=1
    while read metab; do
        met_found=`$AWK -v m=${metab} '{if($1==m) printf"%s\n",$0}' ${_modelinfo_dir}/model_sparse_st_matrix.csv | wc -l`
        if [ ${met_found} -eq 0 ]; then
            feasibility_met=0
            break
        fi
    done < ${lpmfile}

    if [ ${feasibility_met} -eq 0 ]; then
        echo 0
    fi
    
    # Feasibility check successful
    echo 1
}

########
netred()
{
    # Initialize variables
    niter=1
    niter_store=${st_val}
    maxiters=100
    end=0

    # Create temporary directory with metabolic model information
    cp -r ${auto_fba_outdir}/minfo/*.csv $SDIR/curr_minfo
    
    # Execute network reduction loop
    echo "Executing network reduction loop...">&2

    # Obtain set of removable reactions
    obtain_removable_reac > $SDIR/removable_reacs

    while [ $end -eq 0 ]; do
        ## Check ending conditions
        
        # Check maximum number of iterations
        if [ $niter -gt $maxiters ]; then
            echo "Process finished: maximum number of iterations exceeded! (nreac: $nreac ; nremreac: $nremreac ; dof: $dof)" >&2
            echo "" >&2
            end=1
        fi

        # Obtain number of reactions
        obtain_reacs > $SDIR/reacs
        nreac=`obtain_nreac`
        if [ $nreac -eq $mrval ]; then
            echo "Process finished: number of reactions is equal to minimum (nreac: $nreac ; nremreac: $nremreac ; dof: $dof)" >&2
            echo "" >&2
            end=1
        fi

        # Check degrees of freedom
        rank=`obtain_matrix_rank`
        dof=`expr $nreac - $rank`
        if [ $dof -eq $mdval ]; then
            echo "Process finished: degrees of freedom are equal to minimum (nreac: $nreac ; nremreac: $nremreac ; dof: $dof)" >&2
            echo "" >&2
            end=1
        fi

        # Check number of removable reactions
        nremreac=`obtain_nremreac`
        if [ $nremreac -eq 0 ]; then
            echo "Process finished: there are not any removable reactions (nreac: $nreac ; nremreac: $nremreac ; dof: $dof)" >&2
            echo "" >&2
            end=1
        fi

        # Print iteration information
        if [ $end -eq 0 ]; then
            echo "* Iteration $niter (nreac: $nreac ; nremreac: $nremreac ; dof: $dof)" >&2
        fi

        # Prune network if ending conditions were not met
        if [ $end -eq 0 ]; then

            # Create lp file
            case $crit in
                0)
                    biomass_fva
                    ;;
                1)
                    shlomi_fva
                    ;;
            esac

            # Obtain file with sorted flux ranges
            obtain_sorted_flux_ranges

            # Internal while loop
            echo "- Executing internal loop..." >&2
            success=0

            while [ $success -eq 0 -a $nremreac -ne 0 ]; do

                # Obtain candidate reaction for removal
                reac=`obtain_cand_reac_for_removal`
                if [ $reac = "NONE" ]; then
                    echo "... No additional candidates are available" >&2
                    break
                else
                    echo "... Trying to remove reaction $reac..." >&2
                fi
                
                # Remove reaction from set of removables
                remove_reac_from_remov $reac

                # Obtain number of removable reactions
                nremreac=`obtain_nremreac`

                # Remove reaction from network and store it in an
                # auxiliary directory
                remove_reac_from_nw $reac $SDIR/curr_minfo $SDIR/curr_minfo_aux

                # Check protected functions
                success=`check_feasibility $SDIR/curr_minfo_aux`

                # Check success
                if [ $success -eq 1 ]; then
                    # Replace current network with auxiliary network
                    rm -rf $SDIR/curr_minfo
                    cp -r $SDIR/curr_minfo_aux $SDIR/curr_minfo
                fi
                                
                # Clear curr_minfo_aux directory
                rm $SDIR/curr_minfo_aux/*

            done

            echo "" >&2

            # Check whether to store partial result
            if [ `expr $niter % ${niter_store}` -eq 0 ]; then
                cp -r $SDIR/curr_minfo $outd/minfo_iter${niter}
            fi

            # Increase number of iterations
            niter=`expr $niter + 1`
        fi
        
    done

    # Copy result to output directory
    cp -r $SDIR/curr_minfo/* $outd/
}

########
if [ $# -lt 1 ]; then
    echo "Use: network_reducer [-pr <int>] -a <string> -lpm <string> -lpr <string>"
    echo "                     -md <int> -mr <int> -o <string> [-li <int>]"
    echo "                     [-sf <string>] [-g <float>] [-rt <float>] [-st <int>]"
    echo "                     [-c <int>] [-qs <string>] [-sdir <string>] [-debug]"
    echo ""
    echo "-pr <int>      : number of processors"
    echo "-a <string>    : directory storing the output of auto_fba tool"
    echo "-lpm <string>  : file with list of protected metabolites (id's should be given)"
    echo "-lpr <string>  : file with list of protected reactions (id's should be given)"
    echo "-md <int>      : minimum degrees of freedom"
    echo "-mr <int>      : minimum number of reactions"
    echo "-o <string>    : output directory"
    echo "-li <int>      : execute fva at each iteration over a list of <int> randomly"
    echo "                 selected reactions"
    echo "-sf <string>   : Sample file name chosen from those analyzed with auto_fba"
    echo "                 (required if auto_fba was executed with -c 1 option)"
    echo "-g <float>     : value of the gamma parameter (between 0 and 1, 1 by default)"
    echo "-rt <float>    : relative tolerance gap (0.01 by default)"
    echo "-st <int>      : store model after <int> iterations (1000 by default)"
    echo "-c <int>       : sorting criterion for flux ranges (ascending order by default)"
    echo "                 0 -> ascending order"
    echo "                 1 -> descending order"
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
    li_given=0
    c_given=0
    crit=0
    sf_given=0
    g_given=0
    g_val=1
    rt_given=0
    rt_val=0.01
    st_given=0
    st_val=1000
    c_given=0
    sort_crit=0
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
        "-li") shift
            if [ $# -ne 0 ]; then
                li_val=$1
                li_given=1
            fi
            ;;
        "-sf") shift
            if [ $# -ne 0 ]; then
                sample_file=$1
                sf_given=1
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
        "-st") shift
            if [ $# -ne 0 ]; then
                st_val=$1
                st_given=1
            fi
            ;;
        "-c") shift
            if [ $# -ne 0 ]; then
                sort_crit=$1
                c_given=1
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

    if [ ${li_given} -eq 1 ]; then
        echo "-li parameter is ${li_val}" >> ${outd}/params.txt
    fi

    if [ ${sf_given} -eq 1 ]; then
        echo "-sf parameter is ${sample_file}" >> ${outd}/params.txt
    fi

    if [ ${g_given} -eq 1 ]; then
        echo "-g parameter is ${g_val}" >> ${outd}/params.txt
    fi

    if [ ${rt_given} -eq 1 ]; then
        echo "-rt parameter is ${rt_val}" >> ${outd}/params.txt
    fi

    if [ ${st_given} -eq 1 ]; then
        echo "-st parameter is ${st_val}" >> ${outd}/params.txt
    fi

    if [ ${c_given} -eq 1 ]; then
        echo "-c parameter is ${sort_crit}" >> ${outd}/params.txt
    fi

    if [ ${pr_given} -eq 1 ]; then
        echo "-pr parameter is ${nprocs}" >> ${outd}/params.txt
    fi

    # check presence of cplex
    if [ ! -f ${CPLEX_BINARY_DIR}/cplex ]; then
        echo "Error, CPLEX binary not found (shell variable CPLEX_BINARY_DIR should be defined)">&2
        exit 1
    else
        echo "CPLEX_BINARY_DIR=${CPLEX_BINARY_DIR}" >> ${outd}/params.txt
    fi

    ### Process parameters

    # create shared directory
    SDIR="${sdir}/network_reducer_$$"
    mkdir $SDIR || { echo "Error: shared directory cannot be created"  >&2 ; exit 1; }
 
    # create shared subdirectories
    mkdir $SDIR/lp
    mkdir $SDIR/curr_minfo
    mkdir $SDIR/curr_minfo_aux

    # remove temp directories on exit
    if [ $debug -eq 0 ]; then
        trap "rm -rf $SDIR 2>/dev/null" EXIT
    fi

    # Create output directory
    create_out_dir ${outd}

    # Obtain fba criterion
    crit=`obtain_fba_criterion ${auto_fba_outdir}`

    # Initialize random numbers seed
    RANDOM=31415
    
    # Execute network reducer algorithm
    netred || exit 1

fi
