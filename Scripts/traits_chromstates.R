library(data.table)


chromsets <- snakemake@input[['chrom_credsets']]
out <- snakemake@output[['out']]

## chromsets <- list.files("/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_chrom_stats_b38",full.names=T)

##' Combine chromatin annotations by trait
##'
##' @param chromsets files with chromatin annotation by block
##' @param out file to save output
##' @return savefiles with chrom annotations per trait
##' chrom_trait()

chrom_trait <- function(chromsets, out){

    n <- c()
    c_trait <- list()
    for (f in chromsets){
        l <- readRDS(f)
        if(length(l) == 2){            
            if(is.null(intersect(n, l$trait$trait)) | !length(intersect(n, l$trait$trait))){
                ## add new trait info to n and c_trait
                new <- l$trait$trait                
                full_new <- paste0(new, ".tsv.gz")
                c_trait <- append(c_trait, l$chromstates[full_new])
                n <- c(n, new)
            } else {
                ## look for intersect and append to c_trait  
                i <- intersect(n, l$trait$trait)
                full_i <- paste0(i, ".tsv.gz")
                for(j in full_i) {
                    c_trait[[j]] <- rbind(c_trait[[j]], l$chromstates[[j]])
                }
                ## add no intersects to c_trait, if applicable
                new <- setdiff(l$trait$trait, n)
                if(length(new)){
                    full_new <- paste0(new, ".tsv.gz")
                    c_trait <- append(c_trait, l$chromstates[full_new])
                    n <- c(n, new)
                }
            }
        }
    }
    
    saveRDS(c_trait, out)
            
                
            
}

chrom_trait(chromsets, out)
