# Author: Daniel Ortiz Mart\'inez
# *- bash -*

########
function gen_stats_biomass()
{
    # Take parameters
    local_file=$1

    # Print header
    echo "# Stats for file ${local_file}"
    echo ""

    # Obtain stats
    num_metab=`$GREP "constraint name" ${local_file} | $GREP \"_ | wc -l | $AWK '{printf"%s",$1}'`
    num_reactions=`$GREP "variable name" ${local_file} | $GREP \"v | wc -l | $AWK '{printf"%s",$1}'`
    total_num_vars=`$GREP "variable name" ${local_file} | wc -l | $AWK '{printf"%s",$1}'`
    sol_status=`$GREP "solutionStatusString=" ${local_file} | $AWK -F "\"" '{printf"%s",$2}'`
    obj_val=`$GREP "objectiveValue=" ${local_file} | $AWK -F "\"" '{printf"%s",$2}'`

    # Print stats
    echo "## Problem summary"
    echo ""
    echo "- **Number of metabolites**: ${num_metab}"
    echo "- **Number of reactions**: ${num_reactions}"
    echo "- **Total number of variables in the lp problem**: ${total_num_vars}"
    echo ""

    echo "## Solution summary"
    echo ""
    echo "- **Solution status**: ${sol_status}"
    echo "- **Objective value**: ${obj_val}"
}

########
function gen_stats_shlomi()
{
    # Take parameters
    local_file=$1

    # Print header
    echo "# Stats for file ${local_file}"
    echo ""

    # Obtain stats
    num_metab=`$GREP "constraint name" ${local_file} | $GREP \"_ | wc -l | $AWK '{printf"%s",$1}'`
    num_reactions=`$GREP "variable name" ${local_file} | $GREP \"v | wc -l | $AWK '{printf"%s",$1}'`
    hexp_reactions=`$GREP "variable name" ${local_file} | $GREP \"ymh | wc -l | $AWK '{printf"%s",$1}'`
    lexp_reactions=`$GREP "variable name" ${local_file} | $GREP \"ypl | wc -l | $AWK '{printf"%s",$1}'`
    total_num_vars=`$GREP "variable name" ${local_file} | wc -l | $AWK '{printf"%s",$1}'`
    sol_status=`$GREP "solutionStatusString=" ${local_file} | $AWK -F "\"" '{printf"%s",$2}'`
    obj_val=`$GREP "objectiveValue=" ${local_file} | $AWK -F "\"" '{printf"%s",$2}'`
    perc_act_hl_react=`echo "" | $AWK -v ov=${obj_val} -v hr=${hexp_reactions} -v lr=${lexp_reactions} '{printf"%.1f",(ov/(hr+lr))*100}'`

    # Print stats
    echo "## Problem summary"
    echo ""
    echo "- **Number of metabolites**: ${num_metab}"
    echo "- **Number of reactions**: ${num_reactions}"
    echo "- **Highly expressed reactions**: ${hexp_reactions}"
    echo "- **Lowly expressed reactions**: ${lexp_reactions}"
    echo "- **Total number of variables in the lp problem**: ${total_num_vars}"
    echo ""

    echo "## Solution summary"
    echo ""
    echo "- **Solution status**: ${sol_status}"
    echo "- **Objective value**: ${obj_val}"
    echo "- **Percentage of actual hi/lo expressed reactions**: ${perc_act_hl_react}"
}

########
if [ $# -lt 1 ]; then
    echo "Use: gen_fba_stats -f <string> [-c <int>]"
    echo ""
    echo "-f <string>   : path to file with results"
    echo "-c <int>      : fba criterion used to obtain the results. The criteria"
    echo "                used can be one of following,"    
    echo "                0 -> Maximize biomass"    
    echo "                1 -> Shlomi et al. 2008"    
    echo ""
else
    
    # Read parameters
    f_given=0
    c_given=0
    crit=0
    while [ $# -ne 0 ]; do
        case $1 in
        "-f") shift
            if [ $# -ne 0 ]; then
                file=$1
                f_given=1
            fi
            ;;
        "-c") shift
            if [ $# -ne 0 ]; then
                crit=$1
                c_given=1
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

    if [ ! -f ${file} ]; then
        echo "Error! ${file} file does not exist" >&2
        exit 1
    fi

    ### Process parameters
    
    # Generate basic statistics
    case $crit in
        0)
            gen_stats_biomass $file
            ;;
        1)
            gen_stats_shlomi $file
            ;;
    esac


fi
