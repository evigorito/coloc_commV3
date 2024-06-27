library(data.table)

## Get inputs and output:
fm <- snakemake@input[['fm']]
out <- snakemake@output[['out']]

(load(fm))
## get block name from fm
block <- gsub("fm_(.*)\\.RData", "\\1", basename(fm))

########################
## loads logpp and ind
########################
## FM was done by Chris. She said:
## ind shows how each dataset was finemapped - might be susie, cojo or single causal variant. (I prefer susie if it gives results, but sometimes the results look dodgy, then I prefer cojo, unless dodgy, otherwise single. If neither susie nor cojo find multiple causal variants, I use single).

## logpp is a list, one element per fine mapped disease in the block
## each element in the list is a matrix.

##' QC fine mapped blocks, no more than 10 traits allowed
##'
##' @param logpp list with fine mapped SNPs for IMDs within a given block
##' @param block name of fm block
##' @param out file name to save output
##' @export
##' @return saves  a data table with column trait, the name of the trait, and second element a list of data tables. Each data table has the overlaps with the chromatine states for the snps in the cred set(s) per chrom state each sample
##' chrom_fm_cs()

chrom_fm_cs <- function(chromstates, logpp, credsize=0.95, block, out){

    ## ## QC ind: Chris said: If you have 10 or more traits then something is dodgy and I missed it.
    ## traits <- names(ind)

    ## dt <- data.table(trait = sub(".tsv.gz", "", traits))
    
    ## if(length(traits) >= 10 || length(traits) == 0){
        
    ##     saveRDS(list(traits=dt), out)
        
    ## }
