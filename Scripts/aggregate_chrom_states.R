library(data.table)


inp <- snakemake@input[['chromstates']]
out <- snakemake@output[['out']]

## chromstates <- list.files("/home/ev250/rds/rds-mrc-bsu/ev250/communities/chrom18states_b38", "_18_CALLS_segments.bed.gz", full.names=T)


##' Aggregate chromstates from immune cells
##'
##' @param chromstates files with chrom states per cell type per individual
##' @param out file name to save output
##' @export
##' @return 
##' aggregate_chromstates()

aggregate_chromstates <- function(chromstates, out){

    ## name files for chromstates with sample name
    names(chromstates) <- sub("_.*", "", basename(chromstates))

    dt <- rbindlist(lapply(chromstates, function(i) {

        temp=fread(i)
        ## Each state in V4 is repeated V3-V2 times
        temp[, V5:=V3-V2]
        ## Count states
        temp[, sum(V5), V4]

    }), idcol="sample")


    states <- unique(dt$V4)

    ## aggregate

    act <- unlist(lapply(c("Enh[G,A,W]", "Tss[F,A]", "ZNF", "Tx"), function(i) grep(i, states, value=T)) )

    inact <- unlist(lapply(c("EnhB", "TssB","Rep"), function(i) grep(i, states, value=T)))

    quies <- unlist(lapply(c( "Het","Quies"), function(i) grep(i, states, value=T)))

    dt[, agg_state := "quiescence"][V4 %in% act, agg_state := "active"][ V4 %in% inact, agg_state := "inactive"]

    setnames(dt, "V4", "state")

    write.table(dt, out, row.names=F)
}

aggregate_chromstates(inp, out)

    
    
