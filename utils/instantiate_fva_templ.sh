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
if [ $# -lt 1 ]; then
    echo "Use: instantiate_fva_templ -f <string> -d <int> -v <string> -s <float>"
    echo "                           [-g <float>]"
    echo ""
    echo "-f <string>   : fva template file in lp format"
    echo "-d <int>      : optimization direction (0: Min, 1: Max)"
    echo "-v <string>   : variable to be analyzed"
    echo "-s <float>    : solution of the fba problem"
    echo "-g <float>    : value of the gamma parameter (between 0 and 1, 1 by default)"
    echo ""
else
    
    # Read parameters
    f_given=0
    d_given=0
    v_given=0
    s_given=0
    g_given=0
    g_val=1
    while [ $# -ne 0 ]; do
        case $1 in
        "-f") shift
            if [ $# -ne 0 ]; then
                file=$1
                f_given=1
            fi
            ;;
        "-d") shift
            if [ $# -ne 0 ]; then
                d_val=$1
                d_given=1
            fi
            ;;
        "-v") shift
            if [ $# -ne 0 ]; then
                v_val=$1
                v_given=1
            fi
            ;;
        "-s") shift
            if [ $# -ne 0 ]; then
                s_val=$1
                s_given=1
            fi
            ;;
        "-g") shift
            if [ $# -ne 0 ]; then
                g_val=$1
                g_given=1
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

    if [ ${d_given} -eq 0 ]; then
        echo "Error! -d parameter not given" >&2
        exit 1
    fi

    if [ ${v_given} -eq 0 ]; then
        echo "Error! -v parameter not given" >&2
        exit 1
    fi

    if [ ${s_given} -eq 0 ]; then
        echo "Error! -s parameter not given" >&2
        exit 1
    fi

    ## Process parameters

    # Obtain right hand side
    rhv=`echo "" | "$AWK" -v g_val=${g_val} -v s_val=${s_val} '{printf"%f",s_val*g_val}'`

    # Instantiate template
    if [ ${d_val} -eq 0 ]; then
        "$SED" 's/<GOAL>/Minimize/g' "$file" | "$SED" "s/<VAR>/${v_val}/g" | "$SED" "s/<RH>/${rhv}/g"
    fi

    if [ ${d_val} -eq 1 ]; then
        "$SED" 's/<GOAL>/Maximize/g' "$file" | "$SED" "s/<VAR>/${v_val}/g" | "$SED" "s/<RH>/${rhv}/g"
    fi
fi
