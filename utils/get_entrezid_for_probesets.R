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
    cat("get_entrezid_for_probesets [arguments]
 
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

# Process parameters

## Load ExpressionSet file
esetname=load(f)

## Obtain annotation
annotation = annotation(get(esetname))

## Load library for annotation
annot_db=sprintf("%s.db",annotation)
library(annot_db,character.only=TRUE)

## Obtain entrezid names
ID = featureNames(get(esetname))
entrezid =as.character(lookUp(ID, annot_db, "ENTREZID"))

## Obtain correspondence matrix
idmap=cbind(ID,entrezid)

## Write result
write.table(idmap, o, sep=",",col.names=FALSE,quote=FALSE,row.names=FALSE)
