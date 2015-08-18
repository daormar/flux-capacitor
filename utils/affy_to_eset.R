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
    cat("affy_to_eset [arguments]
 
Arguments:
-d     <string>        path to directory with CEL files
-p     <string>        file with phenotype data
-o     <string>        output file
--help                 print this text\n")
  
    q(save="no",status=0)
}

## Check -d option
if(!"d" %in% keys)
{
    print("Error: -d option not provided")
    q(save="no",status=1)
}

## Check -p option
if(!"p" %in% keys)
{
    print("Warning: -p option not provided")
    p_option=FALSE
} else
{
    p_option=TRUE
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
library(annotate)

# Process parameters

## Get current wd
currwd=getwd()

## Load CEL files
setwd(d)
samples=ReadAffy()

## Generate ExpressionSet with normalized values
eset=rma(samples)

## Add phenotype data if provided
if(p_option)
{
    pdata=read.table(p)
    pData(eset) = pdata
}

## Obtain annotation
annotation = annotation(eset)

## Print annotation
cat("\n")
cat("* Annotation: ", annotation, "\n" )

## Save result
setwd(currwd)
save(eset,file=o)
