# Author: Daniel Ortiz Mart\'inez
# *- r -*

# Basic libraries
library(R.utils)

## Collect arguments
args <- commandArgs(asValue=TRUE, excludeReserved=TRUE)[-1]
 
# Turn arguments into R variables
keys <- attachLocally(args)

# Check arguments

## Check --help option
if("help" %in% keys || length(keys)==0)
{    
    cat("calc_matrix_rank [arguments]
 
Arguments:
-m     <string>        path to file with R sparse matrix
--help                 print this text\n")
  
    q(save="no",status=0)
}

## Check -m option
if(!"m" %in% keys)
{
    print("Error: -m option not provided")
    q(save="no",status=1)
}

# Specific libraries
library(Matrix)
library(methods)

# Process parameters

# Read matrix
mat=readMM(m)

# Calculate rank
rank=rankMatrix(mat,tol=NULL,method="qr.R")

# Print result
print(rank[1])
