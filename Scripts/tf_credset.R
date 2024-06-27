library(data.table)
library(motifbreakR)
library(BSgenome)
library(XtraSNPlocs.Hsapiens.dbSNP144.GRCh38) # look up for rsids
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)     # hg38 genome
library(MotifDb)

## inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

fm <- args[1]
credsize <- as.numeric(args[2])
out <- args[3]


## fm <- "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/fine-mapping/fm_chr3_block20.RData"

## credsize <-  0.95

(load(fm))
## get block name from fm
block <- gsub("fm_(.*)\\.RData", "\\1", basename(fm))

##' Map SNPs in credsets to TF binding sites
##'
##' @param snps dbSNP dataset to llok up for rsid
##' @param logpp list with fine mapped SNPs for IMDs within a given block
##' @param ind character vector with names traits and value the method used for fine mapping
##' @param credsize size for credible set, defaults to 0.95
##' @param block name of fm block
##' @param out file name to save output
##' @export
##' @return saves list with each element the TF for the named credible set
##' tf_cset()

tf_cset <- function(snps=SNPlocs.Hsapiens.dbSNP155.GRCh38, logpp, ind, credsize=0.95, block, out){

    message("Starting")
    traits <- names(ind)

    dt <- data.table(trait = sub(".tsv.gz", "", traits))
    
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
    message("part2")
    snps_cs <- unique(Reduce(union, lapply(cs, function(i) unlist(lapply(i, names)))))
    ## Format for genomic ranges
    snps4gr <- paste0("chr", gsub("_", ":", snps_cs))

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
    message("part3")
    g <- GRanges(seqnames=sub(":.*", "", snps_cs),
                 ranges=IRanges(as.numeric(sub(".*:([0-9]+)_.*", "\\1", snps_cs)), width=1),
                 strand="*",
                 SNP_id=snps4gr,
                 REF=gsub(".*([A-Z])_.*", "\\1", snps_cs),
                 ALT=gsub(".*_([A-Z]$)", "\\1", snps_cs)
                 )

    o <- snpsByOverlaps(snps,g)

    snps.mb <- snps.from.rsid(rsid = as.data.table(o)[,RefSNP_id],
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)
    message("part4")
    snps.dt <- as.data.table(snps.mb)

    ## subset snps.mb to contain only the snps_cs

    snps.dt[, id := paste0(sub("chr", "", seqnames), ":", start, "_", REF, "_", ALT)]
    ## get rows with relevant SNPS
    rows <- snps.dt[id %in% snps_cs, which=T]
    
    message("Finding motifs")
    
    results <- motifbreakR(snpList = snps.mb[rows], filterp = TRUE,
                       pwmList = subset(MotifDb, 
                                        dataSource %in% c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C")),
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam())
    
    r <- as.data.table(results)

    if(nrow(r)){
        r <- r[effect == "strong",]
        message("results")
        if(nrow(r)){

        
            r[, chrom := sub("chr", "", seqnames)][, snp := paste0(chrom, ":", start, "_", REF, "_", ALT)]

            ## Add trait, block and pp to r
            tf.trait <- lapply(trait, function(i)  merge(i, r, by="snp", all.x=T))
            
        } else {

            tf.trait <- trait

        }
    } else {

        tf.trait <- trait
    }
    
        
    
    saveRDS(tf.trait, file=out)
                    
}



tf_cset(snps=SNPlocs.Hsapiens.dbSNP155.GRCh38, logpp, ind, credsize, block, out)




