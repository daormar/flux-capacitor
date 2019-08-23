# Author: Daniel Ortiz Mart\'inez
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
