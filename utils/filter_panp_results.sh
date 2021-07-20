# Flux Capacitor package
# Copyright (C) 2015-2018 Daniel Ortiz-Mart\'inez
#  
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.
  
# *- bash -*

########
filter_file()
{
    # Take parameters
    _pfile=$1
    _gfile=$2

    # Filter file
    cat "${_pfile}" | "$AWK" -F "," -v gfile="$_gfile" 'BEGIN{
                                               while( (getline <gfile) > 0)
                                               {
                                                probe_to_eid[$1]=$3
                                                probe_overall_score[$1]=$8
                                                if($8 > best_overall_score[$3])
                                                {
                                                  best_overall_score[$3]=$8
                                                  probe_with_best_overall_scr[$3]=$1
                                                }
                                               }
                                              }
                                              {
                                               if(FNR==1)
                                                 printf"%s\n",$0
                                               else
                                               {
                                                # Print entry if probe mapped in probe_to_eid
                                                if($1 in probe_to_eid)
                                                {
                                                  # Print entry if probe has the best score
                                                  eid=probe_to_eid[$1]
                                                  if(probe_with_best_overall_scr[eid]==$1)
                                                    printf"%s\n",$0
# printf"%s %s %f\n",$1,eid,best_overall_score[eid]
                                                }
                                               }
                                              }'
}

########
if [ $# -lt 1 ]; then
    echo "Use: filter_panp_results -p <string> -g <string>"
    echo ""
    echo "-p <string>   :  file with panp results (generated using exec_panp_eset)"
    echo "-g <string>   :  file with greatest j-scores for probesets"
    echo ""
else
    
    # Read parameters
    p_given=0
    g_given=0
    while [ $# -ne 0 ]; do
        case $1 in
        "-p") shift
            if [ $# -ne 0 ]; then
                pfile=$1
                p_given=1
            fi
            ;;
        "-g") shift
            if [ $# -ne 0 ]; then
                gfile=$1
                g_given=1
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

    if [ ! -f "${pfile}" ]; then
        echo "Error! ${pfile} file does not exist" >&2
        exit 1
    fi

    if [ ${g_given} -eq 0 ]; then
        echo "Error! -g parameter not given" >&2
        exit 1
    fi

    if [ ! -f "${gfile}" ]; then
        echo "Error! ${gfile} file does not exist" >&2
        exit 1
    fi

    # Process parameters
    filter_file "$pfile" "$gfile"

fi
