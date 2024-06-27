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

match_fm <- args[1]
out <- args[2]


## match_fm <- "/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_match_snps/matched_snps/match4credset_fm_chr11_block63.txt"


## get block name from fm
## block <- gsub("fm_(.*)\\.RData", "\\1", basename(fm))

##' Map SNPs in credsets to TF binding sites
##'
##' @param snps dbSNP dataset to look up for rsid
##' @param match_fm file name with output from rule maf_match with matched snps for fm signals for a block
##' @param out file name to save output
##' @export
##' @return saves list with each element the TF for the named credible set
##' tf_match()

tf_match <- function(snps=SNPlocs.Hsapiens.dbSNP155.GRCh38, match_fm, out){

    message("Starting")
    traits <- fread(match_fm)
    
    ## Get rsids for matched snps, create Granges object 
    message("get rsids")
    u.match <- unique(traits$match)
    snps4gr <- paste0("chr", gsub("_", ":", u.match))
    
    g <- GRanges(seqnames=sub(":.*", "", u.match),
                 ranges=IRanges(as.numeric(sub(".*:([0-9]+)_.*", "\\1", u.match)), width=1),
                 strand="*",
                 SNP_id=snps4gr,
                 REF=gsub(".*([A-Z])_.*", "\\1", u.match),
                 ALT=gsub(".*_([A-Z]$)", "\\1", u.match)
                 )

    o <- snpsByOverlaps(snps,g)

    snps.mb <- snps.from.rsid(rsid = as.data.table(o)[,RefSNP_id],
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)
    
    snps.dt <- as.data.table(snps.mb)

    ## subset snps.mb to contain only the snps_cs

    snps.dt[, id := paste0(sub("chr", "", seqnames), ":", start, "_", REF, "_", ALT)]
    ## get rows with relevant SNPS
    rows <- snps.dt[id %in% u.match, which=T]
    
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
        if(nrow(r)){        
            r[, chrom := sub("chr", "", seqnames)][, match := paste0(chrom, ":", start, "_", REF, "_", ALT)]

            ## Add trait, block, cs and match snp annotations to r
            tf.trait <-  merge(traits, r, by="match", all.x=T, allow.cartesian=T)
            
        } else {

            tf.trait <- traits

        }
    } else {

        tf.trait <- traits
    }
           
    
    saveRDS(tf.trait, file=out)
                    
}

tf_match(snps=SNPlocs.Hsapiens.dbSNP155.GRCh38, match_fm, out)




