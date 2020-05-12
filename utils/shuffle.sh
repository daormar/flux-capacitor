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

# Randomly reorders the lines of a file given a seed

shuffle()
{
    # Initialize variables
    local seed=$1
    local file=$2

    # Shuffle file
    $AWK -v seed=$seed 'BEGIN{srand(seed)}{printf"%f %d %s\n",rand(),NR,$0}' $file \
        | LC_ALL=C $SORT -k1n -k2n | ${CUT} -d' ' -f3-
}

if [ $# -ne 1 -a $# -ne 2 ]; then
    echo "Usage: shuffle <seed> [<file>]"
else

    # Take parameters
    seed=$1

    if [ $# -eq 2 ]; then
        file=$2
    fi

    # Invoke shuffling function
    shuffle $seed $file

fi
