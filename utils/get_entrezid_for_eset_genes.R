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
    cat("get_entrezid_for_eset_genes [arguments]
 
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
library(annotate)
library(hgu133plus2.db)

# Process parameters

## Load ExpressionSet file
esetname=load(f)

## Obtain entrezid names (currently works for hgu133plus2 annotation)
#annotation(get(esetname))
ID = featureNames(get(esetname))
entrezid =as.character(lookUp(ID, "hgu133plus2.db", "ENTREZID"))

## Obtain correspondence matrix
idmap=cbind(ID,entrezid)

## Write result
write.table(idmap, o, sep=",",col.names=FALSE,quote=FALSE,row.names=FALSE)
#rm(list=ls())
