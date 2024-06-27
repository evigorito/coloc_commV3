#' ---
#' title: Comparing clustering of immune cell types based on chromatin annotations linked to IMD fine mapped GWAS
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

##'  Inputs

inp <- snakemake@input[['inp']]
meta <- snakemake@input[['meta']]
chrom_agg <- snakemake@input[['chrom_agg']]

#' <!--
## inp='/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_chrom_stats_b38/traits_prob_chromstates.rds'
## meta <- "/home/ev250/rds/rds-mrc-bsu/ev250/communities/chrom18states/main_metadata_table.tsv"
## chrom_agg <- "/home/ev250/rds/rds-mrc-bsu/ev250/communities/chrom18states_b38_aggregated/summary_chrom_states_sample.txt"
#'"-->
#'
##' ## Rationale
##'
##' For each fine-mapped signal of an IMD, I looked at the proportion
##' of quiescence/active/inactive chromatin states in 84 primary human
##' immune cells from the EpiMap repository. The selected samples were
##' B, various T and myeloid cell subsets, NK cells, lymphocytes,
##' mononuclear cells, MPP and CD34 CMP.
##'
##' I first looked at the distribution of active/quiescence/inactive
##' states across samples per IMD. Then, I clustered the samples based
##' on the similarity of the chromatin states. Here, I used the
##' Hellinger distance between samples. Each signal is summarised
##'
##'  

l <- readRDS(inp)

meta_samples <- fread(meta)

## select samples
m <- meta_samples[GROUP %in% c('HSC & B-cell', 'Blood & T-cell') & id !='BSS00313', ]


## Select traits with 25 or more signals

t <- l[sapply(l,  length) >= 25]

## Look at the distribution of chromstates across fine mapped signals, per IMD

dt.t <- setNames(lapply(names(t), function(i) {
    
    dt=as.data.table(Reduce(rbind, lapply(l[[i]], rowMeans)))
    setnames(dt, names(dt), c("Act","Q", "Inact"))
    melt(dt, measure.vars=names(dt))

}), names(t))

## Get chrom states across genome per sample

gw <- fread(chrom_agg)
gw.sample <- gw[, .(N=sum(V1)), .(agg_state, sample)]
tot.sample <- gw[,.(Nt=sum(V1)), sample]
prop.sample <- merge(gw.sample, tot.sample, by="sample")[, Proportion := N/Nt]
prop.sample[, state := "Q"][agg_state == "active", state := "Act"][agg_state == "inactive", state := "Inac"]

plots.t <- lapply(names(dt.t), function(p) {
    ggplot(dt.t[[p]], aes(x=variable, y=value)) +
        geom_boxplot() +
        theme_bw() +
        labs(x="", title=sub("([A-Z]+)_.*", "\\1", p), y="Proportion",
             subtitle=paste("N of signals:", nrow(dt.t[[p]])/3))
})



## Compare with the same samples but genome wide
p.gw <- ggplot(prop.sample, aes(x=state, y=Proportion)) +
    geom_boxplot() +
    theme_bw() +
    labs(x="", title="Genome-wide", y="Proportion")

plot_grid(plotlist=c(plots.t, list(p.gw)))

## aggregate signals per trait
mat.t <- lapply(t, function(i) Reduce(rbind, i))

##' #### Cluster by samples
## CREATE HEATMAPS PER TRAITS

## group T, myeloid cells and germinal centre/lypphocyes and mononuclear as "Other"
m[ , celltype := infoline][grep("(^T|T CELL)", infoline), celltype := "T CELL"]
m[grep("CD14|MYELOID|NEUTROPHIL", infoline), celltype := "MYELOID"]
m[grep("GERMINAL|LYMPHO|MONONUCLEAR", infoline), celltype := "OTHER"]

## To ease visualisation for the heatmaps I aggregated cell subsets,
## the number of cell types are:

 m[,.N,.(celltype, infoline)][order(celltype),]

## use celltype2 for annotation 
df <- data.frame(row.names=m$id,  celltype2=m$celltype)

## create adj matrix of distance between samples

calc_hellinger=function(d1,d2) 
    sum((sqrt(d1)-sqrt(d2))^2)


## create adj matrix of distance between samples

calc_hellinger=function(d1,d2) 
    sum((sqrt(d1)-sqrt(d2))^2)

adj.mat <- setNames(lapply(names(mat.t), function(j) {
    data <- mat.t[[j]]

    todo=expand.grid(1:ncol(data), 1:ncol(data))
    todo=todo[ todo[,1] < todo[,2], ]

    DH <- matrix(NA, nrow=ncol(data), ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
    for(i in 1:nrow(todo)) {
        DH[todo[i,1], todo[i,2]] <- calc_hellinger(data[, todo[i,1]], data[, todo[i,2]])
        DH[todo[i,2], todo[i,1]] <- DH[todo[i,1], todo[i,2]]
    }
    diag(DH) <- 0

    ## scale to 1
    DH <- DH/max(DH)
    rownames(DH) <- colnames(DH)

    ## convert to adj matrix
    return(1-DH)
}
), names(mat.t))




grid.newpage()
heat <- lapply(names(adj.mat), function(j) {
    ## get adj matrix
    data <- adj.mat[[j]]
    ## get the number of independent signals for trait j:
    n.sig <- nrow(mat.t[[j]])/3

    hc=hclust(as.dist(1-data), method="complete")
    t=gsub("([A-Z]+)_.*", "\\1", j)
    
    p <- pheatmap(data, annotation_col=df, show_rownames=F, show_colnames=F,annotation_legend = T, silent=T, main=paste(t, "\n", "Signals:", n.sig), cluster_cols=hc,cluster_rows=hc)
    grid.draw(p)
    grid.newpage()
    ## Look at the order of samples
    m2 <- m[match(colnames(data)[p$tree_col$order], id),.(id, celltype, infoline)]
    return(m2[, trait := j])
}
)

##' #### Cluster by signals
##'
##'
##' Here I calculate Hellinger distance as above.
##'
##' Also, I require quiescence <0.9 in at least one vector, otherwise
##' the distance is NA. This is to avoid close proximity based only on
##' quiescence.
##'
##' 
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

    ## ## swap meaning of 1 and -1 (this is to allow close distance when there is high correlation between inactive vs active states across samples
    ## DHrev=lapply(1:nrow(todo), function(i) calc_hellinger(data.l[[ todo[i,1] ]], data.l[[ todo[i,2] ]][3:1,]))  %>%
    ##     do.call("rbind",.)
    
    ## ## keep minimum distance
    ## DH = pmin(DH,DHrev)

    p=0.9
    ## celltypes with quiescent < 0.9 in at least one member of a pair, otherwise NA
    N90.use=lapply(1:nrow(todo), function(i) { data.l[[ todo[i,1] ]][2,] < p | data.l[[ todo[i,2] ]][2,] < p })  %>%
        do.call("rbind",.)
    ## counts
    N90=rowSums(N90.use)


    ## Make matrix of SNP pair-wise distance 
    MH=matrix(0,length(data.l),length(data.l))
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

sig.plots <- lapply(names(sig.mat), function(i) {

    hc=hclust(as.dist(1-sig.mat[[i]]), method="complete")
    t=gsub("([A-Z]+)_.*", "\\1", i)
    p <- pheatmap(sig.mat[[i]],  show_rownames=F, show_colnames=F, main=t, annotation_legend = F, cluster_cols=hc,cluster_rows=hc)
    grid.draw(p)
    grid.newpage()
})

##' #### Cluster by signals and samples
##'
##' Here I cluster each trait by samples and signals. To order the
##' samples and signals I first cluster all traits together. For each
##' signal in each sample I make "hard" calls for chromatin states, I
##' select the state with the highest probability.
##'
##' When a signal is shared across many traits I take the average
##'

## Look at the number of total signals across traits and unique signals

## total
length(Reduce("c", sapply(t, names)))

## unique
length(unique(Reduce("c", sapply(t, names))))


## Select the expected value for each signal in each sample, by trait

exp.state <- setNames(lapply(names(mat.t), function(i) {

    data <- mat.t[[i]]

    ## split by signal
    data.l <- lapply(seq(1, nrow(data), 3), function(s) data[s:(s+2),])
    ## mat.hard <- Reduce(rbind, lapply(data.l, function(j) apply(j, 2, function(k) as.numeric(rownames(j)[k==mean(k)]))))
    mat.exp <- Reduce(rbind, lapply(data.l, function(j) colSums(j * as.numeric(rownames(j)))))
    rownames(mat.exp) <- sapply(t[i], names)
    dt <- data.table(mat.exp, keep.rownames=T)
    names(dt)[1] <- "signal"
    return(dt)

}), names(mat.t))

## Aggregate exp.state across traits

agg.exp <-  rbindlist(exp.state)

## Collapse signals
collapse.exp <- agg.exp[,lapply(.SD,mean),by=signal]

agg.p <- pheatmap(collapse.exp[, !"signal", with=F], annotation_col=df, annotation_legend = T, show_rownames=F, show_colnames=F, silent=T, main="All traits")

rows <- collapse.exp$signal[agg.p$tree_row$order]
cols <- names(collapse.exp[, !"signal", with=F])[agg.p$tree_col$order]

grid.draw(agg.p)
grid.newpage()

##  SAMPLE ORDER WITH FULL ANNOTATION
m[order(match(id, cols)), .(id, celltype, infoline)]

## Plot by trait, cluster by signal and take the mean of the expected
## chrom states values across cell types

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

    ## mat <- as.matrix(data[, !"signal", with=F], rownames=data$signal)
    
    ## cluster by signals in each trait
    hc=hclust(as.dist(1-sig.mat[[j]]), method="complete")

    ## orderbycol <- m[order(m$celltype),.(id, infoline, celltype)]
    
    
    p <- pheatmap(data[,rownames(df2)], annotation_col=df2, show_rownames=F, show_colnames=F,annotation_legend = T, silent=T, main=paste(t, "\n", "Signals:", n.sig), cluster_cols=F,cluster_rows=hc)
    grid.draw(p)
    grid.newpage()

    })

## Plot by trait, cluster by signal and look at the distribution
## across cell types

exp.distribution <- lapply(exp.state, function(i) {

    ## aggregate states by cell type
    tmp <- as.matrix(i[, !"signal", with=F])
    ## split tmp by cell type
    samples4cells <- setNames(lapply(unique(m$celltype), function(u) m[celltype == u, id]),
                             unique(m$celltype))
    
    sum <- lapply(names(samples4cells), function(u) {
        tmp2 <- t(apply(tmp[,samples4cells[[u]]], 1, summary))
        colnames(tmp2) <- paste(u, colnames(tmp2), sep="_")
        return(tmp2)
    })
    sum.mat <- Reduce(cbind, sum)
    rownames(sum.mat) <- i$signal
    return(sum.mat)

}
)

df2 <- data.frame(row.names=colnames(exp.distribution[[1]]), celltype=sub("(.*)_.*", "\\1",
                                                                          colnames(exp.distribution[[1]])))


trait.p <- lapply(names(exp.distribution), function(j) {
    
    data <- exp.distribution[[j]]
    ## get the number of independent signals for trait j:
    n.sig <- nrow(data)

    t=gsub("([A-Z]+)_.*", "\\1", j)

    ## mat <- as.matrix(data[, !"signal", with=F], rownames=data$signal)
    
    ## cluster by signals in each trait
    hc=hclust(as.dist(1-sig.mat[[j]]), method="complete")

    ## orderbycol <- m[order(m$celltype),.(id, infoline, celltype)]
    
    
    p <- pheatmap(data[,rownames(df2)], annotation_col=df2, show_rownames=F, show_colnames=F,annotation_legend = T, silent=T, main=paste(t, "\n", "Signals:", n.sig), cluster_cols=F,cluster_rows=hc)
    grid.draw(p)
    grid.newpage()

    })




##' ## Conclusions
##'
##' The distribution of chromatin states across samples per IMD showed
##' that active and quiescence states tend to be the most
##' abundant. For some traits active was more prominent but for others
##' quiescence. Inactive was always a minority. Compared to the
##' distribution of chromatin states genome wide in the same samples
##' showed an enrichment on active states, or reduction of quiescent,
##' which suggest we are enriching on active regions.
##'
##' Using the GWAS fine-mapped SNPs in IMDs cluster cell types to some
##' extent, but for most IMDs not very convincing. I think the most
##' clear are RA and CEL, and perhaps AITD to some extent.
##'
##' RA seemed to cluster most T cells together, although with
##' some of the B cells in between. NK are close to T cells and
##' 'Other' can be expected within T cells as they are not very well
##' lineage defined.
##'
##' I expected that the IMDs with more signals or with higher
##' distribution of active states to give better definition of the
##' cell populations, but was not a clear trend. RA has 50 signals and
##' CEL 29. Chris also made the point that the (un)certainty on the
##' mapping of the regions needs to be considered.
##'
##' Anyway, what I havent tested is whether using the genome wide
##' annotations will cluster this set of samples by cell type (positive
##' control). I am just making the assumption that chromatin states
##' should group samples by cell types.
##'
##' Clustering traits by signals showed connected regions in all
##' traits. However, it is not easy to interpret what the clustering
##' of specfic regions of the genome mean. I can look at the
##' transciption factor binding site predicted disruption and see if
##' it helps interpretation.
##'
##' Clustering the traits by signals and samples showed that signals
##' are the main driver for clustering. 
##'
##' Next, I will look at clustering traits by shared pathways based on
##' transcription binding site disruption by credible set of the fine
##' mapped SNPs.
##'
##' 
