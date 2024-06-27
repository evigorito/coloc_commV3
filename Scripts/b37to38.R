library(data.table)
library(rtracklayer)
library(GenomicRanges)


##' Convert genomic coordinates bewteen builds
##'
##' @param bed bed file with input coordinates
##' @param bi input build, defaults 37
##' @param bo oupput build, defaults 38
##' @param out name of file to output bed file with new coordinates
##' @return saves file
##' build_convert()

build_convert(bed, bi="37", bo="38", out){

    system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")

    ## chain file from UCSC is 0-based half-open intervals
    ch = import.chain(chain)
    bed.dt <- fread(bed)
    
    grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=226061851, end=226071523))
    gr <- GRanges(
        seqnames = bed.dt$V1,
        ranges = IRanges(bed.dt$V2, end = bed.dt$V3,
                         state=bed.dt$V4,
                         
    score = 1:10,
    GC = seq(1, 0, length=10))

    results <- as.data.frame(liftOver(grObject, ch))
    }


bed="/home/ev250/rds/rds-mrc-bsu/ev250/communities/chrom18states/BSS00093_18_CALLS_segments.bed.gz"
chain="/home/ev250/rds/rds-mrc-bsu/ev250/communities/chain_file/hg19ToHg38.over.chain"
