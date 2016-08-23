# Author: Daniel Ortiz Mart\'inez
# *- bash -*

########
exclude_readonly_vars()
{
    ${AWK} -F "=" 'BEGIN{
                         readonlyvars["BASHOPTS"]=1
                         readonlyvars["BASH_VERSINFO"]=1
                         readonlyvars["EUID"]=1
                         readonlyvars["PPID"]=1
                         readonlyvars["SHELLOPTS"]=1
                         readonlyvars["UID"]=1
                        }
                        {
                         if(!($1 in readonlyvars)) printf"%s\n",$0
                        }'
}

########
exclude_bashisms()
{
    $AWK '{if(index($1,"=(")==0) printf"%s\n",$0}'
}

########
write_functions()
{
    for f in `${AWK} '{if(index($1,"()")!=0) printf"%s\n",$1}' $0`; do
        $SED -n /^$f/,/^}/p $0
    done
}

########
create_script()
{
    # Init variables
    local name=$1
    local command=$2

    # Write environment variables
    set | exclude_readonly_vars | exclude_bashisms > ${name}

    # Write functions if necessary
    $GREP "()" ${name} -A1 | $GREP "{" > /dev/null || write_functions >> ${name}

    # Write PBS directives
    echo "#PBS -o ${name}.o\${PBS_JOBID}" >> ${name}
    echo "#PBS -e ${name}.e\${PBS_JOBID}" >> ${name}
    echo "#$ -cwd" >> ${name}

    # Write command to be executed
    echo "${command}" >> ${name}

    # Give execution permission
    chmod u+x ${name}
}

########
launch()
{
    local file=$1
    local outvar=$2

    ### qsub invocation
    if [ "${QSUB_WORKS}" = "no" ]; then
        $file &
        eval "${outvar}=$!"
    else
        local jid=$($QSUB ${QSUB_TERSE_OPT} ${qs_opts} $file | ${TAIL} -1)
        eval "${outvar}='${jid}'"
    fi
    ###################
}

########
all_procs_ok()
{
    # Init variables
    local files="$1"
    local sync_num_files=`echo "${files}" | $AWK '{printf"%d",NF}'`    
    local sync_curr_num_files=0

    # Obtain number of processes that terminated correctly
    for f in ${files}; do
        if [ -f ${f}_end ]; then
            sync_curr_num_files=`expr ${sync_curr_num_files} + 1`
        fi
    done

    # Return result
    if [ ${sync_num_files} -eq ${sync_curr_num_files} ]; then
        echo "1"
    else
        echo "0"
    fi
}

########
sync()
{
    # Init vars
    local files="$1"
    local job_ids="$2"

    if [ "${QSUB_WORKS}" = "no" ]; then
        wait
        sync_ok=`all_procs_ok "$files"`
        if [ $sync_ok -eq 1 ]; then
            return 0
        else
            return 1
        fi
    else
        pbs_sync "$files" "${job_ids}"
    fi
}

########
job_is_unknown()
{
    nl=`$QSTAT ${QSTAT_J_OPT} ${jid} 2>&1 | $GREP -e "Unknown" -e "do not exist" | wc -l`
    if [ $nl -ne 0 ]; then
        echo 1
    else
        echo 0
    fi
}

########
pbs_sync()
{
    local files="$1"
    local job_ids="$2"
    local num_pending_procs=0
    end=0
    while [ $end -ne 1 ]; do
        sleep 3
        end=1
        # Check if all processes have finished
        for f in ${files}; do
            if [ ! -f ${f}_end ]; then
                num_pending_procs=`expr ${num_pending_procs} + 1`
                end=0
            fi
        done

        # Sanity check
        num_running_procs=0
        for jid in ${job_ids}; do
            job_unknown=`job_is_unknown ${jid}`
            if [ ${job_unknown} -eq 0 ]; then
                num_running_procs=`expr ${num_running_procs} + 1`
            fi
        done
        if [ ${num_running_procs} -eq 0 -a ${num_pending_procs} -ge 1 ]; then
            return 1
        fi
    done
}

########
fva_for_vlist_frag()
{
    echo "** Processing fragment ${fragm} (started at "`date`")..." >> $SDIR/log || \
        { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

    echo "** Processing fragment ${fragm} (started at "`date`")..." >> $SDIR/${fragm}_proc.log || \
        { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

    # Process flux variables
    while read fvar; do
        # Instantiate fva templates
        $bindir/instantiate_fva_templ -f ${fva_templ} -d 0 -v ${fvar} \
            -s ${fba_sol} -g ${g_val} 2>> $SDIR/${fragm}_proc.log > ${outd}/${fvar}_min.lp || \
            { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

        $bindir/instantiate_fva_templ -f ${fva_templ} -d 1 -v ${fvar} \
            -s ${fba_sol} -g ${g_val} 2>> $SDIR/${fragm}_proc.log > ${outd}/${fvar}_max.lp || \
            { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

        # Configure cplex read mst file command
        read_mst_comm=""
        polishing_comm=""
        if [ ${m_given} -eq 1 ]; then
            read_mst_comm="read ${mst}"

            if [ ${po_given} -ne 0 ]; then
                polishing_comm="set mip polishafter solutions ${po_val}"
            fi
        fi

        # Solve lp problems
        workmem_param=2048
        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/${fvar}_min.lp" "set mip tolerances mipgap ${rt_val}" \
            "set workmem ${workmem_param}" "${read_mst_comm}" "${polishing_comm}" "set timelimit ${tl_val}" \
            "optimize" "write ${outd}/${fvar}_min.sol" \
            2>> $SDIR/${fragm}_proc.log > ${outd}/${fvar}_min.log || \
            { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

        ${CPLEX_BINARY_DIR}/cplex -c "read ${outd}/${fvar}_max.lp" "set mip tolerances mipgap ${rt_val}" \
            "set workmem ${workmem_param}" "${read_mst_comm}" "${polishing_comm}" "set timelimit ${tl_val}" \
            "optimize" "write ${outd}/${fvar}_max.sol" \
            2>> $SDIR/${fragm}_proc.log > ${outd}/${fvar}_max.log || \
            { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

        # Check that sol files where generated (cplex seems to return
        # zero even when it fails)
        if [ ! -f ${outd}/${fvar}_max.sol -or ! -f ${outd}/${fvar}_max.sol ]; then
            echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm} (sol files could not be generated)" >> $SDIR/log; 
            return 1
        fi

        # Add ranges and other info to result file
        min_objv=`$GREP "Objective = " ${outd}/${fvar}_min.log | $AWK '{printf"%s\n",$NF}'`
        time_min=`$GREP "Solution time = " ${outd}/${fvar}_min.log | $AWK '{printf"%s\n",$4}'`
        max_objv=`$GREP "Objective = " ${outd}/${fvar}_max.log | $AWK '{printf"%s\n",$NF}'`
        time_max=`$GREP "Solution time = " ${outd}/${fvar}_max.log | $AWK '{printf"%s\n",$4}'`
        diff=`echo "${min_objv} ${max_objv}" | $AWK '{printf"%f",$2-$1}'`
        echo "$fvar min: ${min_objv} (time: ${time_min} s) ; max: ${max_objv} (time: ${time_max} s) ; diff: $diff" \
            2>> $SDIR/${fragm}_proc.log >> $SDIR/${fragm}.results || \
            { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

        # Compress files
        $GZIP -f ${outd}/${fvar}_min.lp ${outd}/${fvar}_max.lp \
              ${outd}/${fvar}_min.sol ${outd}/${fvar}_max.sol || \
            { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

    done < $SDIR/${fragm} || return 1

    # Write date to log file
    echo "Processing of fragment ${fragm} finished ("`date`")" >> $SDIR/log || \
        { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

    echo "Processing of fragment ${fragm} finished ("`date`")" >> $SDIR/${fragm}_proc.log || \
        { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

    echo "" > $SDIR/qs_fva_${fragm}_end || \
        { echo "Error while executing fva_for_vlist_frag for $SDIR/${fragm}" >> $SDIR/log; return 1 ; }

}

########
gen_log_err_files()
{
    # Copy log file to its final location
    if [ -f $SDIR/log ]; then
        cp $SDIR/log ${outd}/log
    fi

    # Generate file for error diagnosing
    if [ -f ${outd}/solve_fva_for_vlist.err ]; then
        rm ${outd}/solve_fva_for_vlist.err
    fi
    for f in $SDIR/*_proc.log; do
        cat $f >> ${outd}/solve_fva_for_vlist.err
    done
}

########
report_errors()
{
    num_err=`$GREP "Error while executing" ${outd}/log | wc -l`
    if [ ${num_err} -gt 0 ]; then
        prog=`$GREP "Error while executing" ${outd}/log | head -1 | $AWK '{printf"%s",$4}'`
        echo "Error during the execution of solve_fva_for_vlist (${prog})" >&2
        echo "File ${outd}/solve_fva_for_vlist.err contains information for error diagnosing" >&2
    else
        echo "Synchronization error" >&2
        echo "File ${outd}/log contains information for error diagnosing" >&2
    fi
}

########
if [ $# -lt 1 ]; then
    echo "Use: solve_fva_for_vlist [-pr <string>] -f <string> -t <string> -s <float>"
    echo "                         -o <string> [-g <float>] [-m <string>] [-rt <float>]"
    echo "                         [-tl <int>] [-po <int>]"
    echo "                         [-qs <string>] [-sdir <string>] [-debug]"
    echo ""
    echo "-pr <int>      : number of processors"
    echo "-f <string>    : file with flux variables to be analyzed"
    echo "-t <string>    : file with fva template"
    echo "-s <float>     : solution of fba problem"
    echo "-o <string>    : output directory"
    echo "-g <float>     : value of the gamma parameter (between 0 and 1, 1 by default)"
    echo "-m <string>    : file with MIP start"
    echo "-rt <float>    : relative tolerance gap provided to lp solver (0.01 by default)"
    echo "-tl <float>    : time limit in seconds (1e6 by default)"
    echo "-po <int>      : use polishing heuristic after obtaining <int> solutions"
    echo "-qs <string>   : specific options to be given to the qsub command"
    echo "                 (example: -qs \"-l pmem=1gb\")"
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
    tl_given=0
    tl_val=1000000
    po_given=0
    sdir=$HOME
    debug=0
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
        "-tl") shift
            if [ $# -ne 0 ]; then
                tl_val=$1
                tl_given=1
            fi
            ;;
        "-po") shift
            if [ $# -ne 0 ]; then
                po_val=$1
                po_given=1
            fi
            ;;
        "-qs") shift
            if [ $# -ne 0 ]; then
                qs_opts=$1
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

    if [ ${m_given} -eq 1 ]; then
        if [ ! -f ${mst} ]; then
            echo "Error! ${mst} file with MIP start does not exist" >&2
            exit 1
        fi
    fi 

    if [ ! -d ${outd} ]; then
        echo "Error: ${outd} directory does not exist" >&2
        exit 1
    fi

    if [ ! -d ${sdir} ]; then
        echo "Error! ${sdir} directory does not exist" >&2
        exit 1
    fi

    # Print parameters
    if [ ${pr_given} -eq 1 ]; then
        echo "-pr parameter is ${nprocs}" > ${outd}/params.txt
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

    # create shared directory
    SDIR="${sdir}/solve_fva_for_vlist_sdir_$$"
    mkdir $SDIR || { echo "Error: shared directory cannot be created"  >&2 ; exit 1; }

    # remove temp directories on exit
    if [ $debug -eq 0 ]; then
        trap "rm -rf $SDIR 2>/dev/null" EXIT
    fi

    # Output info about tracking script progress
    echo "NOTE: see file $SDIR/log to track fva progress" >&2

    # create log file
    echo "*** Parallel process started at: " `date` > $SDIR/log
    echo "">> $SDIR/log

    # process the input

    # check input size
    input_size=`wc ${fvars} 2>/dev/null | ${AWK} '{printf"%d",$1}'`
    if [ ${input_size} -eq 0 ]; then
        echo "Error: input file ${fvars} is empty" >&2
        exit 1
    fi

    if [ ${input_size} -lt ${nprocs} ]; then
        echo "Error: problem too small"  >&2
        exit 1
    fi

    # shuffle input
    echo "Shuffling input: ${fvars}..." >> $SDIR/log
    $bindir/shuffle 31415 ${fvars} > $SDIR/fvars_shuff
    fvars=$SDIR/fvars_shuff

    # split input
    echo "Spliting input: ${fvars}..." >> $SDIR/log
    frag_size=`expr ${input_size} / ${nprocs}`
    ${SPLIT} -l ${frag_size} ${fvars} $SDIR/frag\_ || exit 1

    # Process each fragment
    i=1
    qs_fva=""
    jids=""
    for f in $SDIR/frag\_*; do
        fragm=`${BASENAME} $f`
        
        create_script $SDIR/qs_fva_${fragm} fva_for_vlist_frag || exit 1
        launch $SDIR/qs_fva_${fragm} job_id || exit 1
        qs_fva="${qs_fva} $SDIR/qs_fva_${fragm}"
        jids="${jids} ${job_id}"
        
        i=`expr $i + 1`
    done

    ### Check that all queued jobs are finished
    sync "${qs_fva}" "${jids}" || { gen_log_err_files ; report_errors ; exit 1; }

    # finish log file
    echo "">> $SDIR/log
    echo "*** Parallel process finished at: " `date` >> $SDIR/log

    # Generate file with results
    cat $SDIR/*.results > ${outd}/results

    # Generate log and err files
    gen_log_err_files

fi
