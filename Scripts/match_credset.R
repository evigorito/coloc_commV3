library(data.table)
library(magrittr)
## closest gene to snp
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(Homo.sapiens) # get TSS

## inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

fm <- args[1]
maf <- args[2]
credsize <- as.numeric(args[3])
out <- args[4]

## fm <- "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/fine-mapping/fm_chr10_block62.RData"
## credsize <-  0.95
## maf <-"~/rds/rds-basis-YRWZsjDGyaU/People/CHRIS/coloc-all/byblock/tmp_mafld_chr10_block62.RData"

(load(fm))
## get block name from fm
block <- gsub("fm_(.*)\\.RData", "\\1", basename(fm))
chrom <- gsub("_block.*", "", block)

(load(maf))
## loads MAF and LD, I use MAF

##' Define credset and match snps in credset to snps within the same chromosome by maf and distance to gene
##'
##' @param logpp list with fine mapped SNPs for IMDs within a given block
##' @param ind character vector with names traits and value the method used for fine mapping
##' @param MAF vector with maf for snps in block
##' @param credsize size for credible set, defaults to 0.95
##' @param block name of fm block
##' @param chrom chromosome to look for matches, formatted as 'chrX'
##' @param out file name to save output
##' @export
##' @return saves list with each element the TF for the named credible set
##' match_cset()

match_cset <- function(logpp, ind, credsize=0.95, block, chrom, MAF, out){

    message("Starting")
    traits <- names(ind)

    dt <- data.table(trait = sub(".tsv.gz", "", traits))

     ## QC ind: Chris said: If you have 10 or more signals for a trait then something is dodgy and I missed it.
    traits <- names(ind)


    dt <- rbindlist(lapply(names(logpp), function(i) data.table(trait=sub(".tsv.gz", "",i), signals=nrow(logpp[[i]]))))
    dt[, block := block]

    if(any(dt$signals >= 10)){
        
        write.table(dt, out, row.names=F)
        
    } else {

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

        ## For each trait, create a data table with cred set snp, block, signal and pp cols

        trait <- setNames(lapply(names(cs), function(n) data.table(trait=n,
                                                                   block=block,
                                                                   signal=rep(1:length(cs[[n]]), sapply(cs[[n]], length)),
                                                                   snp=names(unlist(cs[[n]])),

                                                                   ##get pp and add signal
                                                                   pp= unlist(cs[[n]]))[, sig := paste(block, signal, sep="_")]
                                 ),
                          names(cs))
        
        trait.dt <- rbindlist(trait)

        trait.dt[, sig := paste(block, signal, sep="_")]

        ## get allele frequency for snps (use MAF)
        snps.dt <- merge(trait.dt, data.table(snp=names(MAF), af=MAF), by="snp")
        snps.dt[, pos := as.integer(sub("[0-9]+:([0-9]+)_.*", "\\1", snp))]
        

        ## get distance to TSS of closest gene 
        aa=genes(Homo.sapiens, columns=c("ENSEMBL","SYMBOL"))
        ens=sapply(aa$ENSEMBL, paste, collapse=",")
        drop=ens=="NA"
        table(drop)
        aa=aa[!drop]

        ## I am interested in TSS, which is start or end depending on strand
        tss.dt <- as.data.table(aa)[, TSS := start][strand == "-", TSS := end][seqnames==chrom,]
        tss.gr <- GRanges(seqnames=chrom,
                          ranges=IRanges(tss.dt$TSS, width=1),
                          mcols=tss.dt[, .(ENSEMBL, SYMBOL)])

        snps.gr= GRanges(seqnames=chrom,
                         ranges=IRanges(start=snps.dt$pos, end=snps.dt$pos),
                         mcols=snps.dt[ ,-"pos", with=F])

        

        ##' Find distance to nearest gene (or TSS)
        getnearest=function(snps.gr, aa) {
            nearest=distanceToNearest(snps.gr, aa, select = "all")  %>% as.data.table()
            nearest=nearest[!duplicated(queryHits)][order(queryHits)]
            snps.gr$nearestensg=aa[ nearest$subjectHits ]$mcols.ENSEMBL  %>% sapply(., paste, collapse=",")
            snps.gr$nearestsymbol=aa[ nearest$subjectHits ]$mcols.SYMBOL  %>% sapply(., paste, collapse=",")
            snps.gr$distance=nearest$distance
            as.data.table(snps.gr@elementMetadata)
        }

        dt <- getnearest(snps.gr,tss.gr)
        ## Some fm snps for a trait in more than 1 signal, need to add an identifier
        dt[, id := paste(mcols.snp, mcols.sig, sep="-")]

        ##' Find matches avoiding repetition
        u.match <- function(snps, r.match, match.d){            
            for(s in snps){
                if(length(match.d[[s]])){
                    ## exclude used snps
                    ex <- r.match[!is.na(r.match)]
                    if(length(ex)){
                        pool <- match.d[[s]][! match.d[[s]] %in% ex]
                        if(length(pool)) r.match[[s]] <- sample(pool, 1)
                    } else {
                        r.match[[s]] <- sample(match.d[[s]], 1)
                    }
                }
            }
            return(r.match)
        }

        ##' Find matches by maf and distance
        ##' dist is selected distance to match by
        maf.d.match <- function(tmp, matches, dist2gene, dist){
            match.d <- lapply(1:nrow(tmp), function(u) {
                ## select matches for row u (maf match, and then match y dist)
                dt2 <- dist2gene[mcols.snp %in% names(matches[[u]]), ]
                dt2[between(tmp[u, distance], dt2$distance*(1-dist),
                                  dt2$distance*(1+dist)),mcols.snp]
            })
            names(match.d) <- tmp$id
            return(match.d)
        }

        ## For each snp in dt, match by maf (plus/minus 5% each side)
        ## and then look for a snp matching distance to
        ## gene. Exclude snps within cs

        matches.all <- lapply(unique(dt$mcols.trait), function(i) {

            tmp <- dt[ mcols.trait == i,]
            no.cs <- MAF[!names(MAF) %in% tmp$mcols.snp]
            matches <- lapply(tmp$mcols.af, function(u){
                ## allow for 1-af
                no.cs[between(no.cs, u*0.95, u*1.05) | between(no.cs, (1-u)*.95, (1-u)*1.05)]
            })
            names(matches) <- tmp$id
            
            ## get unique matches
            matches.u <- unique(unlist(lapply(matches, function(i) names(i))))
            
            ## get distance to gene for the potential matches
            d <- MAF[matches.u]
            ## get snp position for potential matches
            pos.d <- as.integer(gsub(".*:([0-9]+)_.*", "\\1", names(d)))
            match.gr <- GRanges(seqnames=chrom,
                                ranges=IRanges(pos.d,pos.d),
                                mcols=data.table(snp=names(d), af=d))
            
            dist2gene <- getnearest(match.gr,tss.gr)

            ## match within 5%, 10% or 20% distance (each side of target distance)
            dist <- c(.05, .1, .2)
            match.d <- maf.d.match(tmp, matches, dist2gene, dist[1])
            r.match <- setNames(rep(NA, nrow(tmp)), tmp$id)

            r.match <- u.match(tmp$id, r.match, match.d)
            
            ## try matching at 10% distance
            if(any(is.na(r.match))){

                w <- names(r.match[is.na(r.match)])
                match.d2 <- maf.d.match(tmp[id %in% w,], matches, dist2gene, dist[2])                
                r.match <-  u.match(w, r.match, match.d2)
            }

            if(any(is.na(r.match))){

                w <- names(r.match[is.na(r.match)])
                match.d3 <- maf.d.match(tmp[id %in% w,], matches, dist2gene, dist[3])                
                r.match <-  u.match(w, r.match, match.d3)
            }

            ## if still no match, select the closest snp by distance
            if(any(is.na(r.match))){
                s <- names(r.match[is.na(r.match)])
                ## go by snp to avoid duplication
                for(u in s){
                    ## select potential matches
                    pool <- dist2gene[mcols.snp %in% names(matches[[u]]),]
                    ## calculate distance from matches to index snp
                    pool[, d := abs(tmp[id == u,distance] - distance),]
                    ## exclude used snps
                    ex <- r.match[!is.na(r.match)]
                    if(length(ex)){
                        pool <- pool[!mcols.snp %in% ex,]
                        if(nrow(pool)) {
                             w <- pool[d ==min(d), which=T]
                             r.match[[u]] <- sample(pool[w, mcols.snp], 1)
                        }   
                    } else {
                         w <- pool[d ==min(d), which=T]
                         r.match[[u]] <- sample(pool[w, mcols.snp], 1)
                        
                    }

                }
                
                
                    
            }

            ## combine index snps and matches info
            match.dt <- merge(data.table(snp=names(r.match), match=r.match), dist2gene, by.x="match", by.y="mcols.snp")
            tmp <- merge(match.dt, tmp, by.x="snp", by.y="id", suffixes=c("_match", "_cs"))
            tmp[, snp:=NULL]
            setnames(tmp, names(tmp), gsub("mcols.", "", names(tmp)))
    
            return(tmp)
        })

        
        write.table(rbindlist(matches.all), out, row.names=F)    
    }

}

match_cset(logpp, ind, credsize=0.95, block, chrom, MAF, out)

       


        
    

