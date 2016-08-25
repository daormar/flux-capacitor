# Author: Daniel Ortiz Mart\'inez
# *- bash -*

########
is_absolute_path()
{
    case $1 in
        /*) echo 1 ;;
        *) echo 0 ;;
    esac
}

########
get_absolute_path()
{
    file=$1
    # Check if an absolute path was given
    absolute=`is_absolute_path $file`
    if [ $absolute -eq 1 ]; then
        echo $file
    else
        oldpwd=$PWD
        basetmp=`$BASENAME $PWD/$file`
        dirtmp=`$DIRNAME $PWD/$file`
        cd $dirtmp
        result=${PWD}/${basetmp}
        cd $oldpwd
        echo $result
    fi
}

########
function obtain_array_names()
{
    # Initialize parameters
    _file=$1

    # Obtain array names
    cat ${_file} | $AWK -F "\"" '{
                              if(NR>1)
                              {
                                printf"%s\n",$2
                              }
                             }'
}

########
function obtain_array_types()
{
    # Initialize parameters
    _file=$1

    # Obtain column number for "type" variable
    _coln=`head -1 ${file} | awk '{for(i=1;i<=NF;++i){if(index($i,"type")!=0) printf"%d",i}}'`
    
    # Obtain array types
    cat ${_file} | $AWK -F "\"" -v cn=${_coln} '{
                              if(NR>1)
                              {
                                printf"%s\n",$(cn*3)
                              }
                             }' | LC_ALL=C $SORT | $UNIQ
}

########
function obtain_arrays_of_type()
{
    # Initialize parameters
    _file=$1
    _type=$2

    # Obtain arrays of type
    $GREP ${_type} ${_file} | $AWK -F "\"" '{
                                printf"%s\n",$2
                             }'    
}

########
function obtain_rnaseq_sample_names()
{
    # Initialize parameters
    _file=$1

    # Obtain sample names
    cat ${_file} | $AWK -F "," '{if(FNR>1) printf"%s\n",$2}'
}

########
function create_out_dir()
{
    _dir=$1
    if [ ! -d ${_dir} ]; then
        mkdir ${_dir} || { echo "Error while creating output directories: ${_dir}" >&2; exit 1; }
    fi
}

########
function biomass_crit()
{
    # Obtain model information
    echo "* Generating metabolic model information..." >&2
    echo "" >&2
    if [ ! -d ${outd}/minfo ]; then
        create_out_dir ${outd}/minfo
        $bindir/extract_sbml_model_info -m $mfile -o ${outd}/minfo/model > ${outd}/minfo/extract_sbml_model_info.log 2>&1 || exit 1
    else
        echo "Warning: metabolic model information extraction skipped, $outd/minfo directory already exists" >&2
        echo "" >&2
    fi

    # Create directories required for the rest of the process
    create_out_dir ${outd}/lp
    create_out_dir ${outd}/sol
    create_out_dir ${outd}/stats

    # Create temporary file
    TMPFILE=`mktemp`
    trap "rm -f $TMPFILE 2>/dev/null" EXIT

    # take parameters
    _sample_file=$TMPFILE
    _exp_name="biomass"

    echo "* Processing experiment with label \"${_exp_name}\"..." >&2
    echo "" >&2

    # generate linear programming problem in lp format
    echo "** Generating linear programming problem in lp format..." >&2
    echo "" >&2
    $bindir/create_lp_file -s ${outd}/minfo/model -c 0 > ${outd}/lp/${_exp_name}.lp \
        2> ${outd}/lp/create_lp_file_${_exp_name}.log || exit 1

    # generate template for fva analysis in lp format
    $bindir/create_lp_file -s ${outd}/minfo/model -c 0 --fva > ${outd}/lp/${_exp_name}_fva_template.lp \
        2> ${outd}/lp/create_lp_file_${_exp_name}_fva_template.log || exit 1

    # check presence of cplex
    if [ -f ${CPLEX_BINARY_DIR}/cplex ]; then
        # solve linear programming problem
        echo "** Solving linear programming problem..." >&2
        echo "" >&2

        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/lp/${_exp_name}.lp" \
            "optimize" "write ${outd}/sol/${_exp_name}.sol" \
            > ${outd}/sol/cplex_${_exp_name}.log || exit 1

        # Command line for clp tool:
        # ${CBC_BINARY_DIR}/clp -import ${outd}/lp/${_exp_name}.lp

        # obtain statistics about solution
        echo "** Obtaining solution statistics..." >&2
        echo "" >&2

        if [ -f ${outd}/sol/${_exp_name}.sol ]; then
            $bindir/gen_fba_stats -f ${outd}/sol/${_exp_name}.sol -c 0 > ${outd}/stats/${_exp_name}.md
        else
            echo "** Error, solution file ${outd}/sol/${_exp_name}.sol does not exist" >&2
            echo "" >&2          
        fi

    else
        echo "Warning: cplex binary not found, linear programming problem will not be solved. Please define CPLEX_BINARY_DIR shell variable" >&2
    fi

    echo "" >&2
}

########
function fba_exp_shlomi_marray()
{
    # take parameters
    _sample_file=$1
    _exp_name=$2

    echo "* Processing experiment with label \"${_exp_name}\"..." >&2
    echo "" >&2

    # obtain absent/present genes
    echo "** Obtaining absent/present genes..." >&2
    echo "" >&2
    $bindir/get_absent_present_genes_marray -d  ${outd}/esetdir/probesets_to_entrezids.csv \
        -p ${outd}/esetdir/panp_results_filt.csv -l ${_sample_file} > ${outd}/abs_pres_info/abs_pres_genes_${_exp_name}.csv \
        2>${outd}/abs_pres_info/get_absent_present_genes_marray_${_exp_name}.log || exit 1

    # generate linear programming problem in lp format
    echo "** Generating linear programming problem in lp format..." >&2
    echo "" >&2
    $bindir/create_lp_file -s ${outd}/minfo/model -a ${outd}/abs_pres_info/abs_pres_genes_${_exp_name}.csv \
        -c 1 > ${outd}/lp/${_exp_name}.lp \
        2> ${outd}/lp/create_lp_file_${_exp_name}.log || exit 1

    # generate template for fva analysis in lp format
    $bindir/create_lp_file -s ${outd}/minfo/model -a ${outd}/abs_pres_info/abs_pres_genes_${_exp_name}.csv \
        -c 1 --fva > ${outd}/lp/${_exp_name}_fva_template.lp \
        2> ${outd}/lp/create_lp_file_${_exp_name}_fva_template.log || exit 1

    # check presence of cplex
    if [ -f ${CPLEX_BINARY_DIR}/cplex ]; then

        # solve linear programming problem
        echo "** Solving linear programming problem..." >&2
        echo "" >&2

        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/lp/${_exp_name}.lp" "set mip tolerances mipgap ${rt_val}" \
            "optimize" "write ${outd}/sol/${_exp_name}.sol" \
            "write ${outd}/sol/${_exp_name}.mst all" > ${outd}/sol/cplex_${_exp_name}.log || exit 1

        # Command line for cbc tool:
        # ${CBC_BINARY_DIR}/cbc -import ${outd}/lp/${_exp_name}.lp -ratio ${rt_val} -branchAnd

        # obtain csv and json files with fluxes for solution
        echo "** Obtaining csv and json files with fluxes for solution..." >&2
        echo "" >&2
        $bindir/get_cplex_fluxes -f ${outd}/sol/${_exp_name}.sol -m ${outd}/minfo/model \
            -of 0 > ${outd}/sol/fluxes_${_exp_name}.csv || exit 1
        $bindir/get_cplex_fluxes -f ${outd}/sol/${_exp_name}.sol -m ${outd}/minfo/model \
            -of 1 > ${outd}/sol/fluxes_${_exp_name}.json || exit 1

        # obtain statistics about solution
        echo "** Obtaining solution statistics..." >&2
        echo "" >&2
        if [ -f ${outd}/sol/${_exp_name}.sol ]; then
            $bindir/gen_fba_stats -f ${outd}/sol/${_exp_name}.sol -c 1 > ${outd}/stats/${_exp_name}.md
        else
            echo "** Error, solution file ${outd}/sol/${_exp_name}.sol does not exist" >&2
            echo "" >&2          
        fi

    else
        echo "Warning: cplex binary not found, linear programming problem will not be solved. Please define CPLEX_BINARY_DIR shell variable" >&2
    fi

    echo "" >&2
}

########
function shlomi_crit_marray()
{
    # obtain model information
    echo "* Generating metabolic model information..." >&2
    echo "" >&2
    if [ ! -d ${outd}/minfo ]; then
        create_out_dir ${outd}/minfo
        $bindir/extract_sbml_model_info -m $mfile -o ${outd}/minfo/model > ${outd}/minfo/extract_sbml_model_info.log 2>&1 || exit 1
    else
        echo "Warning: metabolic model information extraction skipped, $outd/minfo directory already exists" >&2
        echo "" >&2
    fi

    # generate expression set
    echo "* Generating expression set..." >&2
    echo "" >&2
    create_out_dir ${outd}/esetdir
    $bindir/affy_to_eset -d $cdir -p $pfile -o ${outd}/esetdir/eset.rda > ${outd}/esetdir/affy_to_eset.log 2>&1 || exit 1

    # obtain entrezid's for probe set ids
    echo "* Obtaining entrez ids for probe set ids..." >&2
    echo "" >&2
    $bindir/get_entrezid_for_probesets -f ${outd}/esetdir/eset.rda \
        -o ${outd}/esetdir/probesets_to_entrezids.csv > ${outd}/esetdir/eids.log 2>&1 || exit 1

    # obtain expression information using panp
    echo "* Obtaining expression information using panp library..." >&2
    echo "" >&2
    $bindir/exec_panp_eset -f ${outd}/esetdir/eset.rda \
        -o ${outd}/esetdir/panp_results.csv > ${outd}/esetdir/exec_panp_eset.log 2>&1 || exit 1

    # obtain file with jetset scores
    echo "* Obtaining file with jetset scores for probesets..." >&2
    echo "" >&2
    annot=`$GREP "Annotation:" ${outd}/esetdir/affy_to_eset.log | $AWK '{printf"%s",$3}'`
    $bindir/get_jetset_scores -c ${annot} \
        -o ${outd}/esetdir/${annot}_jscores.csv > ${outd}/esetdir/get_jetset_scores.log 2>&1

    # filter panp results using jscores
    echo "* Filtering panp expression information using jetset scores..." >&2
    echo "" >&2
    $bindir/filter_panp_results -p ${outd}/esetdir/panp_results.csv \
        -g ${outd}/esetdir/${annot}_jscores.csv > ${outd}/esetdir/panp_results_filt.csv || exit 1

    ##########

    # Create directories required for the rest of the process

    create_out_dir ${outd}/abs_pres_info
    create_out_dir ${outd}/lp
    create_out_dir ${outd}/sol
    create_out_dir ${outd}/stats

    ## Carry out an FBA experiment for all samples:

    array_names=`obtain_array_names $pfile`

    # Create temporary file
    TMPFILE=`mktemp`
    trap "rm -f $TMPFILE 2>/dev/null" EXIT

    if [ ! -f ${CPLEX_BINARY_DIR}/cplex ]; then
        echo "Warning, CPLEX binary not found (shell variable CPLEX_BINARY_DIR should be defined)">&2
    fi

    for arrn in ${array_names}; do

        # write sample info to file
        echo $arrn > $TMPFILE

        # carry out fba experiment
        fba_exp_shlomi_marray $TMPFILE $arrn

    done
}

########
function fba_exp_shlomi_rnaseq()
{
    # take parameters
    _sample_file=$1
    _exp_name=$2

    echo "* Processing experiment with label \"${_exp_name}\"..." >&2
    echo "" >&2

    # obtain absent/present genes
    echo "** Obtaining absent/present genes..." >&2
    echo "" >&2
    $bindir/get_absent_present_genes_rnaseq -r ${rnaseq_cfile} \
        -l ${_sample_file} > ${outd}/abs_pres_info/abs_pres_genes_${_exp_name}.csv \
        2>${outd}/abs_pres_info/get_absent_present_genes_rnaseq_${_exp_name}.log || exit 1

    # generate linear programming problem in lp format
    echo "** Generating linear programming problem in lp format..." >&2
    echo "" >&2
    $bindir/create_lp_file -s ${outd}/minfo/model -a ${outd}/abs_pres_info/abs_pres_genes_${_exp_name}.csv \
        -c 1 > ${outd}/lp/${_exp_name}.lp \
        2> ${outd}/lp/create_lp_file_${_exp_name}.log || exit 1

    # generate template for fva analysis in lp format
    $bindir/create_lp_file -s ${outd}/minfo/model -a ${outd}/abs_pres_info/abs_pres_genes_${_exp_name}.csv \
        -c 1 --fva > ${outd}/lp/${_exp_name}_fva_template.lp \
        2> ${outd}/lp/create_lp_file_${_exp_name}_fva_template.log || exit 1

    # check presence of cplex
    if [ -f ${CPLEX_BINARY_DIR}/cplex ]; then

        # solve linear programming problem
        echo "** Solving linear programming problem..." >&2
        echo "" >&2

        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/lp/${_exp_name}.lp" "set mip tolerances mipgap ${rt_val}" \
            "optimize" "write ${outd}/sol/${_exp_name}.sol" \
            "write ${outd}/sol/${_exp_name}.mst all" > ${outd}/sol/cplex_${_exp_name}.log || exit 1

        # Command line for cbc tool:
        # ${CBC_BINARY_DIR}/cbc -import ${outd}/lp/${_exp_name}.lp -ratio ${rt_val} -branchAnd

        # obtain csv and json files with fluxes for solution
        echo "** Obtaining csv and json files with fluxes for solution..." >&2
        echo "" >&2
        $bindir/get_cplex_fluxes -f ${outd}/sol/${_exp_name}.sol -m ${outd}/minfo/model \
            -of 0 > ${outd}/sol/fluxes_${_exp_name}.csv || exit 1
        $bindir/get_cplex_fluxes -f ${outd}/sol/${_exp_name}.sol -m ${outd}/minfo/model \
            -of 1 > ${outd}/sol/fluxes_${_exp_name}.json || exit 1

        # obtain statistics about solution
        echo "** Obtaining solution statistics..." >&2
        echo "" >&2
        if [ -f ${outd}/sol/${_exp_name}.sol ]; then
            $bindir/gen_fba_stats -f ${outd}/sol/${_exp_name}.sol -c 1 > ${outd}/stats/${_exp_name}.md
        else
            echo "** Error, solution file ${outd}/sol/${_exp_name}.sol does not exist" >&2
            echo "" >&2          
        fi

    else
        echo "Warning: cplex binary not found, linear programming problem will not be solved. Please define CPLEX_BINARY_DIR shell variable" >&2
    fi

    echo "" >&2
}

########
function shlomi_crit_rnaseq()
{
    # obtain model information
    echo "* Generating metabolic model information..." >&2
    echo "" >&2
    if [ ! -d ${outd}/minfo ]; then
        create_out_dir ${outd}/minfo
        $bindir/extract_sbml_model_info -m $mfile -o ${outd}/minfo/model > ${outd}/minfo/extract_sbml_model_info.log 2>&1 || exit 1
    else
        echo "Warning: metabolic model information extraction skipped, $outd/minfo directory already exists" >&2
        echo "" >&2
    fi


    ##########

    # Create directories required for the rest of the process

    create_out_dir ${outd}/abs_pres_info
    create_out_dir ${outd}/lp
    create_out_dir ${outd}/sol
    create_out_dir ${outd}/stats

    ## Carry out an FBA experiment for all samples:

    sample_names=`obtain_rnaseq_sample_names $pfile`

    # Create temporary file
    TMPFILE=`mktemp`
    trap "rm -f $TMPFILE 2>/dev/null" EXIT

    if [ ! -f ${CPLEX_BINARY_DIR}/cplex ]; then
        echo "Warning, CPLEX binary not found (shell variable CPLEX_BINARY_DIR should be defined)">&2
    fi

    for samplen in ${sample_names}; do

        # write sample info to file
        echo $samplen > $TMPFILE

        # carry out fba experiment
        fba_exp_shlomi_rnaseq $TMPFILE $samplen

    done
}

########
if [ $# -lt 1 ]; then
    echo "Use: auto_fba -m <string> [-c <int>] -o <string>"
    echo "              [-d <string>|-r <string> -p <string>] [-rt <float>]"
    echo ""
    echo "-m <string>   : path to file with metabolic model in SBML format"
    echo "-c <int>      : fba criterion used to generate the lp file. The criteria"
    echo "                can be selected from the following list,"    
    echo "                0 -> Maximize biomass"    
    echo "                1 -> Shlomi et al. 2008"    
    echo "-o <string>   : output directory"
    echo "-d <string>   : path to directory with CEL files (required by criterion 1)"
    echo "-r <string>   : file with rna-seq counts (required by criterion 1)"
    echo "-p <string>   : file with phenotype data (required by criterion 1)"
    echo "-rt <float>   : relative tolerance gap (0.01 by default)"
    echo ""
else
    
    # Read parameters
    m_given=0
    d_given=0
    r_given=0
    p_given=0
    o_given=0
    c_given=0
    crit=0
    rt_given=0
    rt_val=0.01
    while [ $# -ne 0 ]; do
        case $1 in
        "-m") shift
            if [ $# -ne 0 ]; then
                mfile=$1
                m_given=1
            fi
            ;;
        "-d") shift
            if [ $# -ne 0 ]; then
                cdir=$1
                d_given=1
            fi
            ;;
        "-r") shift
            if [ $# -ne 0 ]; then
                rnaseq_cfile=$1
                r_given=1
            fi
            ;;
        "-p") shift
            if [ $# -ne 0 ]; then
                pfile=$1
                p_given=1
            fi
            ;;
        "-o") shift
            if [ $# -ne 0 ]; then
                outd=$1
                o_given=1
            fi
            ;;
        "-c") shift
            if [ $# -ne 0 ]; then
                crit=$1
                c_given=1
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
    if [ ${m_given} -eq 0 ]; then
        echo "Error! -m parameter not given" >&2
        exit 1
    fi

    if [ ! -f ${mfile} ]; then
        echo "Error! ${mfile} file does not exist" >&2
        exit 1
    fi

    if [ ${crit} -eq 1 ]; then
        if [ ${d_given} -eq 0 -a ${r_given} -eq 0 ]; then
            echo "Error! -d or -r parameters should be given" >&2
            exit 1
        fi

        if [ ${d_given} -eq 1 -a ${r_given} -eq 1 ]; then
            echo "Error! -d and -r cannot be given simultaneously" >&2
            exit 1
        fi

        if [ ${d_given} -eq 1 ]; then
            if [ ! -d ${cdir} ]; then
                echo "Error! ${cdir} directory does not exist" >&2
                exit 1
            fi
        fi

        if [ ${r_given} -eq 1 ]; then
            if [ ! -f ${rnaseq_cfile} ]; then
                echo "Error! ${rnaseq_cfile} file does not exist" >&2
                exit 1
            else
                rnaseq_cfile=`get_absolute_path ${rnaseq_cfile}`
            fi
        fi
    fi

    if [ ${crit} -eq 1 ]; then
        if [ ${p_given} -eq 0 ]; then
            echo "Error! -p parameter not given" >&2
            exit 1
        else
            pfile=`get_absolute_path ${pfile}`
        fi

        if [ ! -f ${pfile} ]; then
            echo "Error! ${pfile} file does not exist" >&2
            exit 1
        fi
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

    if [ ${d_given} -eq 1 ]; then
        echo "-d parameter is ${cdir}" >> ${outd}/params.txt
    fi

    if [ ${r_given} -eq 1 ]; then
        echo "-r parameter is ${rnaseq_cfile}" >> ${outd}/params.txt
    fi

    if [ ${p_given} -eq 1 ]; then
        echo "-p parameter is ${pfile}" >> ${outd}/params.txt
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o parameter is ${outd}" >> ${outd}/params.txt
    fi

    if [ ${rt_given} -eq 1 ]; then
        echo "-rt parameter is ${rt_val}" >> ${outd}/params.txt
    fi

    echo "-c parameter is ${crit}" >> ${outd}/params.txt


    ### Process parameters
    case $crit in
        0)
            biomass_crit
            ;;
        1)
            if [ ${d_given} -eq 1 ]; then 
                shlomi_crit_marray
            fi

            if [ ${r_given} -eq 1 ]; then 
                shlomi_crit_rnaseq
            fi
            ;;
    esac
    
fi
