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
    cat("exec_panp_eset [arguments]
 
Arguments:
-f     <string>        path to rda file
-o     <string>        output file
--help                 print this text\n")
  
    q(save="no",status=0)
}

## Check -f option
if(!"f" %in% keys)
{
    print("Error: -f option not provided")
    q(save="no",status=1)
}

## Check -o option
if(!"o" %in% keys)
{
    print("Error: -o option not provided")
    q(save="no",status=1)
}

# Specific libraries
library(Biobase)
library(affy)
library(gcrma)
library(panp)

# Process parameters

## Load ExpressionSet file
esetname=load(f)

## Obtain highly and lowly expressed genes
esetPaCalls <- pa.calls(get(esetname))
paCalls <- esetPaCalls$Pcalls
paVals <- esetPaCalls$Pvals
paCalls <- ifelse(paCalls=="P","1",paCalls)
paCalls <- ifelse(paCalls=="A","0",paCalls)
paCalls <- ifelse(paCalls=="M","NA",paCalls)

## Write result
write.table(paCalls, o, sep=",",,col.names=TRUE,quote=FALSE,row.names=TRUE)
