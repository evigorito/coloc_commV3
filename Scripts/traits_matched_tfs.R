library(data.table)
library(qusage)
library(parallel)

source("/home/ev250/communities/Functions/aux_leiden.R")

## inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

paths <- args[1:3]
n <- args[4:6]
out <- args[7]
tfs_match <- args[8:length(args)]


##' Make adj matrix connecting genes through pathways
##'
##' @param tfs_match files with TF binding site disrupted by SNPs within a block (mat_tf2match output)
##' @param paths files with annotation pathways for "h", "kegg" and "reactome"
##' @param n names for paths files (defaults to h, kegg, reactome)
##' @param out file name to save a data table with traits, signal and block linking TFs binding sites predicted to be disrupted by SNPs using pathways
##' @return saves output 
##' match_tf2path()

match_tf2path <- function(tfs_match, paths, n=c("h", "kegg", "reactome"), out){
    
    genes <- rbindlist(lapply(tfs_match, function(i) {
        dt <- readRDS(i)
        if("effect" %in% names(dt)){
            return(dt[effect == 'strong',])
        } else {
            return(dt)
        }
        }), fill=T)

    ## Get pathways
    pth <- lapply(paths, qusage::read.gmt)
    tfs <- unique(genes$geneSymbol[!is.na(genes$geneSymbol)])
    
    tf.paths <- lapply(pth, function(i) {        
        l= lapply(tfs, function(g) path(g, i))
        names(l) <- tfs
        return(l[sapply(l, function(x) !is.null(x))])
    }
    )
    names(tf.paths) <- n

    ## make dt from tf.paths for easy merging
    tf.dts <- lapply(tf.paths, function(i) rbindlist(lapply(names(i), function(j) data.table(geneSymbol=j, pathway=i[[j]]))))

    ## Look at the pathways associated with a tf predicted to be disrupted by a matched snp within a credset for a given trait and combine
    genes[, sig := paste(block, signal, sep="_")]
    
    l <- lapply(tf.dts, function(p) {
        tmp <- merge(genes, p, by="geneSymbol", allow.cartesian=T)
        tmp[, .(score=sum(pp)), by=.(trait, sig, pathway)]
        })

    saveRDS(l, out)

}

match_tf2path(tfs_match, paths, n, out)

## tfs_match <- list.files("/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_match_snps/tf_sites", full.names=T)
## paths <- c("/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/h.all.v2023.1.Hs.symbols.gmt", "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt", "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
## n=c("h", "kegg", "reactome")
