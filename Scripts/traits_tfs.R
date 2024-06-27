library(data.table)
library(qusage)
library(parallel)

source("/home/ev250/communities/Functions/aux_leiden.R")

## inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

tf.ids <- args[1]
paths <- args[2:4]
n <- args[5:7]
out <- args[8]
tfs_cs <- args[9:length(args)]


##' Make adj matrix connecting genes through pathways
##'
##' @param tf.ids file with tfs to use (tf_symbol output)
##' @param tfs_cs files with TF binding site disrupted by SNPs within a community credible set (tf_annot output)
##' @param paths files with annotation pathways for "h", "kegg" and "reactome"
##' @param n names for paths files (defaults to h, kegg, reactome)
##' @param out file name to save a data table with traits, signal and block linking TFs binding sites predicted to be disrupted by SNPs using pathways
##' @return saves output 
##' tf2path()

tf2path <- function(tf.ids, tfs_cs, paths, n=c("h", "kegg", "reactome"), out){
    
    genes <- fread(tf.ids)

    ## Get pathways
    pth <- lapply(paths, qusage::read.gmt)

    tf.paths <- lapply(pth, function(i) {
        l= mclapply(genes$SYMBOL, function(g) path(g, i))
        names(l) <- genes$SYMBOL
        return(l)
    }
    )
    names(tf.paths) <- n

    ## Look at the pathways within a credset and make adj matrix. This way avoids loading all the tfs_cs files at the same time in R
    
    l <- lapply(tf.paths, function(p){
        cs.path <- list()
        for(f in tfs_cs) {
            tmp <- readRDS(f)
            block <- gsub(".*(chr.*)\\.rds", "\\1", basename(f))
            tmp2 <- lapply(tmp, function(i){
                ## make compatible to run tfpath
                setnames(i, "snp", "SNP_id")
                ## split credible set by signal
                res <- lapply(unique(i$signal), function(x) {
                    r <- tfpath(i[signal == x,], p)
                    if(is.data.table(r)) return(r[, c("block", "signal", "trait") := list(block, x, unique(i$trait))])}
                    )
                if(any(!is.null(unlist(lapply(res, class))))) {

                    return(rbindlist(res))
                }
            }
            )                      
                    
            if(any(!is.null(unlist(lapply(tmp2, class))))) cs.path[[block]] <- rbindlist(tmp2)
        }
        return(rbindlist(cs.path))
    })
    

    saveRDS(l, out)

}

tf2path(tf.ids, tfs_cs, paths, n, out)

## tf.ids <- "/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_tf/tf_entrez/entrez_traits_tfs.txt"
## tfs_cs <- list.files("/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_tf/tf_sites", full.names=T)
## paths <- c("/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/h.all.v2023.1.Hs.symbols.gmt", "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt", "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
## n=c("h", "kegg", "reactome")
