#' ---
#' title: Combining chromatin states from immune cells with predicted transcription factor binding site disruption for IMD fine mapped GWAS
#' author: Elena
#' output:
#'    html_document:
#'     toc: true
#'         
#' ---

#' <!--Report built by rule report from /home/ev250/communities/communitiesV2/communitiesV3 Snakefile-->
#'
#' 
##+ message=FALSE
library(data.table)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(cowplot)
library("RColorBrewer")
library(grid)
library(qusage)
library(gridExtra)



#' <!--
## block_qc=list.files("/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_chrom_stats_b38", "fm_.*_chromstates.rds", full.names=T)
## inp='/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_chrom_stats_b38/traits_prob_chromstates.rds'
## meta <- "/home/ev250/rds/rds-mrc-bsu/ev250/communities/chrom18states/main_metadata_table.tsv"
## chrom_agg <- "/home/ev250/rds/rds-mrc-bsu/ev250/communities/chrom18states_b38_aggregated/summary_chrom_states_sample.txt"
## tfs.f='/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_tf/tf_trait/traits_tfs.rds'
## paths <- c("/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/h.all.v2023.1.Hs.symbols.gmt", "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt", "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
#'"-->

##'  Inputs

args <- commandArgs(trailingOnly = TRUE)
inp <- args[2]
meta <- args[3]
chrom_agg <- args[4]
tfs.f <- args[5]
paths <- args[6:8]
block_qc <- args[9:length(args)]

##' ## Rationale
##'
##' Here, for each trait I cluster signals based on the chromatin
##' states in immnune cells. To do this, for the top SNPs within the
##' 95% of a credible set I take the average of the cell state based
##' on the chromatin annotation (active (1), inactive (-1),
##' heterochromatin (0)) weighted by the posterior probability (PP) of
##' causality. Then I calculated the average across cell types (B, T,
##' myeloid, etc). These results are visualised in a per IMD heatmap
##' with rows signals and columns cell type.
##'
##' To get some insight on the biology of those signals I then look at
##' the predicted transcription factor (TF) binding site disruption
##' across the SNPs in the credible set. I then linked the TFs to
##' pathways, here, I only show results using the HALLMARK
##' database. Using the PP of causality I calculated a pathway score
##' for each pathway within a signal. These results are presented as
##' heatmaps with rows signals and columns pathways. The rows are
##' clustered using the chromatin state annotation, to ease comparison
##' with chromatin states.
##'
##' For this analysis I only show IMDs with at least 25 independent
##' signals.


##' Start with chromstates. Make a heatmap which rows are signals and 
##'columns correspond to the average of the expected chromatin state for each 
##'cell type.

l <- readRDS(inp)
meta_samples <- fread(meta)

## select samples
m <- meta_samples[GROUP %in% c('HSC & B-cell', 'Blood & T-cell') & id !='BSS00313', ]
## group T, myeloid cells and germinal centre/lypphocyes and mononuclear as "Other"
m[ , celltype := infoline][grep("(^T|T CELL)", infoline), celltype := "T CELL"]
m[grep("CD14|MYELOID|NEUTROPHIL", infoline), celltype := "MYELOID"]
m[grep("GERMINAL|LYMPHO|MONONUCLEAR", infoline), celltype := "OTHER"]

## Select traits with 25 or more signals

t <- l[sapply(l,  length) >= 25]
## aggregate signals per trait
mat.t <- lapply(t, function(i) Reduce(rbind, i))

##' #### Cluster chromatin states by signals and samples
##'
##' Here I cluster each trait by samples and signals. To order the
##' samples and signals I first cluster all traits together. For each
##' signal in each sample I make "hard" calls for chromatin states, I
##' select the state with the highest probability.
##'


## First need to calculate distance between signals to allow
## clustering

calc_hellinger=function(d1,d2) 
    colSums((sqrt(d1)-sqrt(d2))^2)

sig.mat <- setNames(lapply(names(mat.t), function(j) {
    data <- mat.t[[j]]

    ## split by signal
    data.l <- lapply(seq(1, nrow(data), 3), function(s) data[s:(s+2),])
    todo=expand.grid(1:length(data.l), 1:length(data.l))
    
    todo=todo[ todo[,1] < todo[,2], ]

    
    DH=lapply(1:nrow(todo), function(i) calc_hellinger(data.l[[ todo[i,1] ]], data.l[[ todo[i,2] ]]))  %>%
        do.call("rbind",.)

    p=0.9
    ## celltypes with quiescent < 0.9 in at least one member of a pair, otherwise NA
    N90.use=lapply(1:nrow(todo), function(i) { data.l[[ todo[i,1] ]][2,] < p | data.l[[ todo[i,2] ]][2,] < p })  %>%
        do.call("rbind",.)
    ## counts
    N90=rowSums(N90.use)


    ## Make matrix of SNP pair-wise distance 
    MH=matrix(0,length(data.l),length(data.l), dimnames=list(names(t[[j]]), names(t[[j]])))
    todo  %<>% as.matrix()

    ## Exclude quiescent samples for both SNPs to sum distances
    ## When all samples are quiesc, rowSumns(N90.use) == 0, NAs in MH

    MH[ todo ] = rowSums(DH * N90.use) / rowSums(N90.use)
    MH[ todo[,2:1] ] = rowSums(DH * N90.use) / rowSums(N90.use)


    ## Scale to max=1
    MH <- MH/max(MH, na.rm=T)
    
   

    ## convert to adj matrix
    return(1-MH)
}
), names(mat.t))

## Select the expected value for each signal in each sample, by trait

exp.state <- setNames(lapply(names(mat.t), function(i) {

    data <- mat.t[[i]]

    ## split by signal
    data.l <- lapply(seq(1, nrow(data), 3), function(s) data[s:(s+2),])
   
    mat.exp <- Reduce(rbind, lapply(data.l, function(j) colSums(j * as.numeric(rownames(j)))))
    rownames(mat.exp) <- sapply(t[i], names)
    dt <- data.table(mat.exp, keep.rownames=T)
    names(dt)[1] <- "signal"
    return(dt)

}), names(mat.t))

## Samples
m[order(celltype),.N, .(celltype, infoline)]

exp.mean <- lapply(exp.state, function(i) {

    ## aggregate states by cell type
    tmp <- as.matrix(i[, !"signal", with=F])
    ## split tmp by cell type
    samples4cells <- setNames(lapply(unique(m$celltype), function(u) m[celltype == u, id]),
                             unique(m$celltype))
    
    cell.mean <- lapply(names(samples4cells), function(u) {
        return(matrix(rowMeans(tmp[,samples4cells[[u]]]), ncol=1, dimnames=list(i$signal, u)))
    })
    cellmean.mat <- Reduce(cbind, cell.mean)
    return(cellmean.mat)

}
)

df2 <- data.frame(row.names= unique(m$celltype), celltype= unique(m$celltype))

meantrait.p <- lapply(names(exp.mean), function(j) {
    
    data <- exp.mean[[j]]
    ## get the number of independent signals for trait j:
    n.sig <- nrow(data)

    t=gsub("([A-Z]+)_.*", "\\1", j)
    
    ## cluster by signals in each trait
    hc=hclust(as.dist(1-sig.mat[[j]]), method="complete")

    ## orderbycol <- m[order(m$celltype),.(id, infoline, celltype)]
    
    
    pheatmap(data[,rownames(df2)], annotation_col=df2, show_rownames=F, show_colnames=F,annotation_legend = T, silent=T, main=paste(t, "\n", "Signals:", n.sig), cluster_cols=F,cluster_rows=hc)[[4]] ##to save each plot into a list. note the [[4]]

    })

########################################################################
## For each signal in exp.mean, look at the predicted TFs and pathways
########################################################################

tfs <- readRDS(tfs.f)

## I need to exclude blocks with too many traits, Chris said 10 or
##  more traits per block is dodgy

exclude <- rbindlist(lapply(block_qc, function(k) {
    dt <- readRDS(k)[[1]]
    b <- gsub("fm_(.*)\\_chromstates.rds", "\\1", basename(k))    
    if(any(dt$signals >= 10)) return(dt[, block := b])
}))

    
## Identify independent signals on traits "t", exclude dodgy blocks

tfs <- lapply(tfs, function(i) {
    i <- i[trait %in% names(t) & !block %in% unique(exclude$block),][,sig := paste(block, signal, sep="_")]
    return(i)
})

## Make matrix rows signals and cols pathways, replace NAs with 0

tfs.w <- lapply(tfs, function(i) {
    DT <- dcast(i, trait + sig ~ pathway, value.var = "score")
    for (j in 3:ncol(DT)){
        set(DT,which(is.na(DT[[j]])),j,0)
    }
    return(DT)

})

##' Plot chromstates and pathways side by side ordering rows as
##' clustered by chromstates and columns by pathway similarity

## Get pathways
pth <- lapply(paths, qusage::read.gmt)
names(pth) <- names(tfs)

## Get distance between pathways

jac_idx <- function(a,b){

    length(intersect(a,b)) / length(union(a,b))
}


d.paths <- lapply(pth, function(i) {

    todo=expand.grid(1:length(i), 1:length(i))
    
    todo=todo[ todo[,1] < todo[,2], ]

    J=matrix(NA, nrow=length(i), ncol=length(i), dimnames=list(names(i), names(i)))
    
    for(j in 1:nrow(todo)) {
        J[todo[j,1], todo[j,2]] <- jac_idx(i[[todo[j,1]]], i[[todo[j,2]]])
        J[todo[j,2], todo[j,1]] <- J[todo[j,1], todo[j,2]]
    }
    diag(J) <- 1
    
    return(J)
})

## annotation Hallmark pathways can be grouped by process
## category
## (https://www.cell.com/cell-systems/pdf/S2405-4712(15)00218-5.pdf)

dt.h <- data.table(pathway=colnames(d.paths[[1]]),  process=c("signalling", "pathway", "metabolic", "proliferation", "signalling", "signalling", "signalling","DNA-damage", "proliferation", "pathway", "signalling", "development", rep("signalling",3), "development", "pathway", "immune", "immune", "cellular_component","cellular_component", "signalling", "immune", "pathway", "signalling", "signalling", "proliferation", "proliferation", "proliferation", "development","immune","metabolic", "metabolic", "metabolic", "metabolic", "pathway", "pathway", "DNA_damage","DNA_damage", "development", "metabolic", "immune", "signalling", "metabolic", "cellular_component", "immune", "development", "signalling", "signalling", "development"))


# Plot by trait, signals in the cluster order of chromstates and columns 'hallmark' pathways
cluster.path <- lapply(names(tfs.w)[1], function(i) {
    ## go through traits
    lapply(names(t), function(j) {

        dt <- tfs.w[[i]]
        
        dt <- tfs.w[[i]][trait == j,]
        t=sub("([A-Za-z0-9]+)_.*", "\\1", j)
        
        ## cluster pathways
        cols <- colnames(d.paths[[1]])[colnames(d.paths[[1]]) %in% names(dt)]
        hc.cols=hclust(as.dist(1-d.paths[[1]][cols,cols]), method='complete')

        ## get order for rows
        ## cluster chromstates by signals in each trait
        hc=hclust(as.dist(1-sig.mat[[j]]), method="complete")

        ## add signals in chromstates missing in pathways
        w <- colnames(sig.mat[[j]])[which(!colnames(sig.mat[[j]]) %in% dt$sig)]
        if(length(w)) {
            dt.add <- data.table(trait=j, sig=w)[, (cols) := lapply(seq_along(cols), function(i) 0)]

            dt <- rbind(dt, dt.add)
        }

        df3 <- data.frame(dt.h[pathway %in% cols, ])
        row.names(df3) <- df3$pathway
        df3$pathway <- NULL
        
        ## order cols as in hc.cols but dont show cluster so the plot is the same size as chromatin states
        return(pheatmap(dt[, 3:ncol(dt)][, hc.cols$order, with=F], annotation_col=df3, show_rownames=F, show_colnames=F, silent=T, main=paste(t, "\n", "Signals:", nrow(dt)), cluster_cols=F,cluster_rows=hc)[[4]])
        
        
    }
    )
})

##+ fig.height=5, fig.width=12, fig.align="center"
figs <- lapply(seq_along(meantrait.p), function(i) {
    grid.newpage()
    grid.arrange(arrangeGrob(grobs= list(meantrait.p[[i]],cluster.path[[1]][[i]]), ncol=2))
})


## QC: look at reported SNP regulating ETS2 in macrophages
## 21:39094644

pos <-39094644
chr21 <- block_qc[grep("chr21", block_qc)]

snp <- rbindlist(lapply(chr21, function(i) {
    tmp <- readRDS(i)
    if(length(tmp) == 2){
        dt <- rbindlist(tmp[[2]])
        dt[start == pos,]
    }
}))

snp[, max(pp), trait]
    
snp <- merge(snp, m, by.x="sample", by.y="id")
## same info for any trait
snp[trait == "CD_DeLange_28067908_1-hg38.tsv.gz",.N, .(State, infoline)]

## Monocyte/myeloid and CMP  are associated with enhancer activity

## look at pathways
block <- unique(snp$block)
signal <- unique(snp$signal)
sig.snp <- paste(block, signal, sep="_")

## select top 5
top <- 5
lapply(tfs, function(i) setorder(i[sig == sig.snp,], trait, -score)[, head(.SD, top), trait, .SD=c("trait", "pathway", "score")])


##' #### Conclusions
##'
##' Open chromatin regions are not enriched in annotated pathways.
##'
##' I need to go by trait and match regions that are not
##' heterochronatin/quiesent and match with pathways/tfs.
##'
##' In the report looking at predicted TF binding site disruption
##' across all IMDs, the leading pathways associated to predictions
##' were IFNg, TNFa, etc. I will look if there is any consistency
##' across traits with chromatin state.
##'
##' As QC for the input/analysis I took the SNP recently reported to
##' be causal for regulating the expression of ETS2 in macrophages and
##' affecting risk for various IMDs. In our analysis it is in the
##' credible set for CD, UC, ANS and PSC. It is in a region with
##' enhancer marks for all the monocyte samples (4), CMP and most
##' myeloid. The SNP is not within a TF binding site, it regulates a
##' nearby PU.1 motif, so we would not capture that.

