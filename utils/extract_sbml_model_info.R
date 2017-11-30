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
    cat("extract_sbml_model_info [arguments]
 
Arguments:
-m     <string>        path to file with metabolic model in SBML format
-o     <string>        prefix for output files
--help                 print this text\n")
  
    q(save="no",status=0)
}

## Check -m option
if(!"m" %in% keys)
{
    print("Error: -m option not provided")
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
library(sybilSBML)

# Process parameters

## Read model in SBML format
sbml_model=readSBMLmod(m)

# Write results

## Write compartments
write.table(sbml_model@mod_compart, file=paste(o,"_comp.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## Write gene id's
write.table(allGenes(sbml_model), file=paste(o,"_gene_ids.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## ## Write reaction ids
write.table(sbml_model@react_id, file=paste(o,"_reaction_ids.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## Write gpr boolean vector
write.table(gpr(sbml_model)!="", file=paste(o,"_gpr_bool_vec.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## Write gpr rules
write.table(gprRules(sbml_model), file=paste(o,"_gpr_rules.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## Write metabolite ids
write.table(sbml_model@met_id, file=paste(o,"_metabolite_ids.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## Write metabolite compartments
write.table(sbml_model@met_comp, file=paste(o,"_metabolite_comp.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## Write stoichiometric matrix
writeMM(sbml_model@S, file=paste(o,"_sparse_st_matrix.csv",sep=""))

## Write lower bounds for reactions
write.table(sbml_model@lowbnd, file=paste(o,"_reaction_lowbnds.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## Write upper bounds for reactions
write.table(sbml_model@uppbnd, file=paste(o,"_reaction_uppbnds.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=TRUE)

## Write reaction id of objective function
objf=which(sbml_model@obj_coef==1)
write.table(objf, file=paste(o,"_obj_fun.csv",sep=""), sep=",",col.names=FALSE,quote=FALSE,row.names=FALSE)
