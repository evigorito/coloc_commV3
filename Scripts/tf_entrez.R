library(data.table)
library(mygene)
library("clusterProfiler")
library(org.Hs.eg.db)
library("DBI")

## inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

tf_dir <- args[1]
out.genes <- args[2]

tf_f <- list.files(tf_dir, full.names=T)

## tf_f <- list.files("/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_tf/tf_sites", full.names=T)

##' Get extrez ids for the intersection of TFs mapping predicted binding sites for community SNPs
##'
##' @param tf_f file names with tf info for each fine mapped region
##' @param out.genes file name to output results
##' enz4tf()
##'
enz4tf <- function(tf_f, out.genes){

    ## Get predicted TFs with strong effect
    l <- lapply(tf_f, function(i) {
        r <- rbindlist(readRDS(i))
        if("effect" %in% names(r)) {                       
            unique(r[effect == "strong" , geneSymbol])
        }
    }
    )
    
    tf <- Reduce(function(x,y) {
        union(x, y)
    }, l)

## Get entrezid
    entrez <- as.data.table(bitr(gene=tf, fromType='SYMBOL', toType='ENTREZID', OrgDb="org.Hs.eg.db"))

    ## Look for missing genes, some genes missing if given a alias symbol

    if(length(tf[! tf %in% entrez$SYMBOL])){
        
        queryGeneNames <- tf[! tf %in% entrez$SYMBOL]

                                        # use sql to get alias table and gene_info table (contains the symbols)
                                        # first open the database connection
        dbCon <- org.Hs.eg_dbconn()
                                        # write your SQL query
        sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
                                        # execute the query on the database
        aliasSymbol <- as.data.table(dbGetQuery(dbCon, sqlQuery))
                                        # subset to get your results
        result <- aliasSymbol[alias_symbol %in% queryGeneNames,]

        if(nrow(result)){

            tf2 <- result$symbol
            entrez2 <- as.data.table(bitr(gene=tf2, fromType='SYMBOL', toType='ENTREZID', OrgDb="org.Hs.eg.db"))

             if(nrow(entrez2)){
                 entrez2 <- merge(entrez2, result[, .(alias_symbol, symbol)], by.x="SYMBOL", by.y="symbol")

                 entrez <- rbindlist(list(entrez, entrez2), fill=T)

             }
        }
    }
        
    write.table(entrez, out.genes, row.names=F)
}
    
enz4tf(tf_f, out.genes)
