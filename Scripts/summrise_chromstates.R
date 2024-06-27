library(data.table)
source("/home/ev250/communities/Functions/aux.R")

inp <- snakemake@input[['chrom_trait']]
out <- snakemake@output[['out']]

## inp <- "/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_chrom_stats_b38/traits_chromstates.rds"


##' Summarise chromatin annotations within traits. For each signal in each trait I calculate the probability of active, repressed or quiescent
##'
##' @param inp list with each element the chromatin annotations of a trait
##' @param out full file name to write output
##' @return saves a list, each element for a trait. Then for each trait a list with each element a signal. This is a matrix with rows chromatin state (active, repressed or quiescent), cols samples and values prob of each state. Cols sum to 1.
##' sum_chrom()

sum_chrom <- function(inp, out){

    l <- readRDS(inp)

    ## Number of independent signals per trait
    lapply(l, function(i) nrow(i)/84)

    ## Compute probabilities of active, inactive and quiescent chromatine state integrating info across credible set

    ## Group states by active and inactive

    states <- Reduce(function(a,b) funion(a[,.(State)], b[, .(State)]), l)$State

    act <- unlist(lapply(c("Enh[G,A,W]", "Tss[F,A]", "ZNF", "Tx"), function(i) grep(i, states, value=T)) )

    inact <- unlist(lapply(c("EnhB", "TssB","Rep"), function(i) grep(i, states, value=T)))

    quies <- unlist(lapply(c( "Het","Quies"), function(i) grep(i, states, value=T)))

    ## Summarise within trait
    ## Convert each element of cs.states to wide format and compare SNPs
    ## within the same credible set by annotation

    cs.wide <- lapply(l, function(i){
        ## sort by pp, remove duplicate entries
        d=dcast(unique(i),  id + snp + block + signal + pp ~ sample, value.var="State")
        setorder(d, id, -pp)        
    })

    ## convert each signal of cs.wide to a matrix
    samples <- unique(l[[1]]$sample)

    cs.mat <- lapply(cs.wide, function(i) {
        ## split by block/signal
        signals <- unique(i$id)
        sig_states <- setNames(lapply(signals, function(j) help_states(i[id == j,], samples, act, inact, quies)),
                               signals)
        return(sig_states)
                
    })

    saveRDS(cs.mat, out)

}

sum_chrom(inp, out)
