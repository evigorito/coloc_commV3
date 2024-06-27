library(data.table)




##' Get running time for jobs in hpc
##'
##' @param file file with running time to format
##' @return adds column with running time in minutes per job
##' runtime()

runtime <- function(file){

    dt <- fread(file)
    t <- paste(dt$V5,dt$V6, dt$V7, sep=":") 
    dt[, time := as.ITime(t)]
    ## convert dt to wide, so each file has 2 times
    dt.start <- dt[,.SD[1],V1, .SD=c("time")]
    dt.end <- dt[,.SD[.N],V1,.SD=c("time")]
    dt.w <- merge(dt.start, dt.end, by="V1", suffixes=c(".start", ".end"))
    d <- difftime(dt.w$time.end, dt.w$time.start, unit="mins")
    dt.w[, diff := d]
    write.table(dt.w, file, row.names=F)   

}

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
runtime(file)

## file="/home/ev250/rds/rds-mrc-bsu/ev250/communities/V3/time2run/map_tf2credset.txt"
