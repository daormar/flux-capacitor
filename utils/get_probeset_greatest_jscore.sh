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
    echo "Use: get_probeset_greatest_jscore -f <string>"
    echo ""
    echo "-f <string>   :  file with jetset scores"
    echo ""
else
    
    # Read parameters
    f_given=0
    while [ $# -ne 0 ]; do
        case $1 in
        "-f") shift
            if [ $# -ne 0 ]; then
                jsfile=$1
                f_given=1
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

    if [ ! -f ${jsfile} ]; then
        echo "Error! ${jsfile} file does not exist" >&2
        exit 1
    fi

    # Process parameters
    cat ${jsfile} | $AWK -F ',' '{if($9!="NA") printf"%s\n",$0}' | LC_ALL=C $SORT -t ',' -k9,9 -k8,8gr |\
        $SED 's/,/ /g' | $UNIQ -f 8 | $SED 's/ /,/g'
fi
