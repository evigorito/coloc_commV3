library(data.table)
library(GenomicRanges)

source("/home/ev250/communities/Functions/aux.R")

## Get inputs and output:
chromstates <- snakemake@input[['chromstates']]
fm <- snakemake@input[['fm']]
credsize <- snakemake@params[['credsize']]
out <- snakemake@output[['out']]

## fm <- "~/rds/rds-cew54-wallace-share/People/Chris/Elena/fine-mapping/fm_chr10_block8.RData"

## credsize <-  0.95
## chromstates <- list.files("/home/ev250/rds/rds-mrc-bsu/ev250/communities/chrom18states_b38", "_18_CALLS_segments.bed.gz", full.names=T)

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

## In each logpp matrix there is a column names 'null', which is the hypothesis that there is no causal variant. Chris advised to remove it: I would just do setdiff(., "null") to remove.



##' Map SNPs in credible sets to chromstates from immune cells
##'
##' @param chromstates files with chrom states per cell type per individual
##' @param logpp list with fine mapped SNPs for IMDs within a given block
##' @param credsize size for credible set, defaults to 0.95
##' @param block name of fm block
##' @param out file name to save output
##' @export
##' @return saves list with first element a data table with column trait, the name of the trait, and second element a list of data tables. Each data table has the overlaps with the chromatine states for the snps in the cred set(s) per chrom state each sample
##' chrom_fm_cs()

chrom_fm_cs <- function(chromstates, logpp, credsize=0.95, block, out){

    ## QC ind: Chris said: If you have 10 or more signals for a trait then something is dodgy and I missed it.
    traits <- names(ind)

    dt <- data.table(trait = sub(".tsv.gz", "", traits))

    dt <- rbindlist(lapply(names(logpp), function(i) data.table(trait=sub(".tsv.gz", "",i), signals=nrow(logpp[[i]]))))
    dt[, block := block]

    if(any(dt$signals >= 10)){
        
        saveRDS(list(traits=dt), out)
        
    } else {

        ## Create Granges for the union of snsps across the credible sets in all traits. This is to overlap each snp with chromstates only once.
        ## First define credsets for each imd in the block

        cs <- lapply(logpp, function(i) {
            ## each row in logpp corresponds to an independent signal, so the need to create credsets by row
            ## exclude 'null' column from i
            i <- i[, colnames(i) != "null", drop=FALSE]
            ## transform to PP
            i <- exp(i)
            ## select credset by row, order and cumsum
            lapply(1:nrow(i), function(j) {
                cum <- cumsum(sort(i[j,], decreasing=T))
                ## select credset based on credsize
                w=which(cum >=credsize)[1]
                snps_credset <- names(cum)[ 1:w]
                i[j, snps_credset]
            })
        })
                
        snps <- unique(Reduce(union, lapply(cs, function(i) unlist(lapply(i, names)))))

        ## For each trait, create a data table with cred set snp, block and pp cols

        trait <- setNames(lapply(names(cs), function(n) data.table(trait=n,
                                                                   block=block,
                                                                   signal=rep(1:length(cs[[n]]), sapply(cs[[n]], length)),
                                                                   snp=names(unlist(cs[[n]])),

                                                                   ##get pp 
                                                                   pp= unlist(cs[[n]]))
                                 ),
                          names(cs))
                          
                                 
        ## create Granges object for snps

        g <- GRanges(seqnames=paste0("chr",sub(":.*", "", snps)),
                     ranges=IRanges(as.numeric(sub(".*:([0-9]+)_.*", "\\1", snps)), width=1))

        ## name files for chromstates with sample name
        names(chromstates) <- sub("_.*", "", basename(chromstates))

        message("start overlap")

        
        ## Overlap snps with sample files with chrom states
        o <- rbindlist(lapply(names(chromstates), function(i) state4cs2(snps, g, chromstates[[i]], i)))

        ## For each trait add chrom states info

        trait2 <- lapply(trait, function(i) {
            
            ## having issues merging across many signals for some files, split by signal and rbind
            s <- unique(i$signal)
            rbindlist(lapply(s, function(m) merge(i[signal==m,], o, by="snp", sort=F)[,id := paste(block, signal, sep="_")]))


        })
        
        ## save
        saveRDS(list(trait=dt, chromstates=trait2), out)
        
       
    }
}

chrom_fm_cs(chromstates, logpp, credsize, block, out)
