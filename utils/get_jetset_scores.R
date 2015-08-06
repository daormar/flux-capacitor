# Author: Daniel Ortiz Mart\'inez
# *- r -*

# Basic libraries
library(R.utils)

# Set echo to TRUE
# options(echo=TRUE) 

## Collect arguments
args <- commandArgs(asValue=TRUE, excludeReserved=TRUE)[-1]
 
# Turn arguments into R variables
keys <- attachLocally(args)
#cat("Command-line arguments attached to global environment:\n");
#print(keys);
#str(mget(keys, envir=globalenv()))

# Check arguments

## Check --help option
if("help" %in% keys || length(keys)==0)
{    
    cat("get_jetset_scores [arguments]
 
Arguments:
-c     <string>        chip name (hgu133a, hgu133plus2, hgu95av2, u133x3p)
-o     <string>        output file
--help                 print this text\n")
  
    q(save="no",status=0)
}

## Check -c option
if(!"c" %in% keys)
{
    print("Error: -c option not provided")
    q(save="no",status=1)
}

## Check -o option
if(!"o" %in% keys)
{
    print("Error: -o option not provided")
    q(save="no",status=1)
}

# Specific libraries
## library(Biobase)
library(jetset)
## library(annotate)
## library(hgu133plus2.db)

# Process parameters

## Obtain jetset scores
jscr=jscores(c)

## Write result
write.table(jscr, o, sep=",",col.names=TRUE,quote=FALSE,row.names=TRUE)
#rm(list=ls())
