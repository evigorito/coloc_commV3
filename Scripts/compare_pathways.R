#' ---
#' title: Comparing clustering of IMDs based on the predicted effect of fine mapped GWAS SNPs on disrupting transcription factor  binding sites
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
library(magrittr)
library(qusage)

##'  Inputs

args <- commandArgs(trailingOnly = TRUE)
inp <- args[2]
paths <- args[3:length(args)]

#' <!--
## inp='/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_tf/tf_trait/traits_tfs.rds'
## paths <- c("/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/h.all.v2023.1.Hs.symbols.gmt", "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt", "/home/ev250/rds/rds-cew54-wallace-share/People/Chris/Elena/pathways/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
#'"-->
#'
##' ## Rationale
##'
##' For each fine-mapped signal of an IMD, I looked at the predicted
##' effect of each SNP of a credible set for the disruption of a
##' transcription factor binding site. 
##'
##' Then, I looked at the annotated pathways for those TFs using 3
##' databases: hallmark (h), kegg and reactome.
##'
##' Weigheted by the posterior probability of being causal I assign
##' scores to each pathway within a credible set. Then, I calculated
##' the aggregated score for each pathway across all signals, to see
##' the most prevalent pathways. I also cluster signals by shared
##' pathways. Here, I only plotted the results using "Hallmark"
##' because it is a relatively small number of pathways covering
##' signalling, immune function and some general metabolic pathways.
##' 
##' 

tfs <- readRDS(inp)

## Identify independent signals

tfs <- lapply(tfs, function(i) i[, sig := paste(block, signal, sep="_")])

##' For each trait rank pathways aggregating posterior probabilities
##' across signals

p <- lapply(tfs, function(i) i[,.N, .(pathway, trait)])

path.t <- lapply(tfs, function(i) i[, .(score=sum(score)), .(pathway,trait)])

## Look at the number of signals per pathway

sigpath <- lapply(tfs, function(i) {
    sig.path <- i[,.(Nsig_trait=.N),.(pathway, trait)]
    sig <- i[,.N, .(sig, trait)][, .(Nsig=.N), trait]
    merge(sig.path, sig, by="trait")
})

## Look at the number of genes per pathway

## Get pathways
pth <- lapply(paths, qusage::read.gmt)
names(pth) <- names(tfs)

genes.path <- lapply(pth, function(p){
    ## get number of genes per pathway
    x <- sapply(p, length)
    data.table(pathway=names(x), N=x)
})

names(genes.path) <- names(path.t)


## Plot pathways per trait, by pathway dataset

path.p <- lapply(names(path.t), function(i) {

    dt <- path.t[[i]]
    dt <- merge(dt, genes.path[[i]], by="pathway")
    ## add number of signals per pathway
    dt <- merge(dt, sigpath[[i]], by=c("trait", "pathway"))
    lapply(unique(dt$trait), function(u) {
        dt2 <- dt[trait == u,]
        setorder(dt2, score)
        dt2[ , pathway := factor(pathway, levels=unique(pathway))]
        t=gsub("([A-Z]+)_.*", "\\1", u)
        
        ## Reactome has too many pathways, take top 20
        if(i == "reactome"){
            if(nrow(dt2) > 20){
                dt2 <- tail(dt2, n=20)
            }
        }        
        ggplot(dt2, aes(score, pathway)) +
            geom_point() +
            theme_bw() +
            xlim(NA, max(dt2$score)*1.2) + 
            geom_text(aes(x=max(dt2$score)*1.2, y=pathway, label=N)) +
            geom_text(aes(x=max(dt2$score)*1.1, y=pathway, label=Nsig_trait)) +
            labs(title=paste(t, "with", i))
    })
})


##' Here I plot the aggregated score for each pathway across all
##' signals. The numbers on the first column on the right hand side of
##' the plot correspond to the number of signals that share the
##' pathway, while the second column to the number of genes annotated
##' in each pathway.

path.p[[1]]




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

## Plot by trait, cluster by signal using 'hallmark' pathways
## annotation Hallmark pathways can be grouped by process
## category
## (https://www.cell.com/cell-systems/pdf/S2405-4712(15)00218-5.pdf)

df.h <- data.frame(row.names=colnames(d.paths[[1]]),  process=c("signalling", "pathway", "metabolic", "proliferation", "signalling", "signalling", "signalling","DNA-damage", "proliferation", "pathway", "signalling", "development", rep("signalling",3), "development", "pathway", "immune", "immune", "cellular_component","cellular_component", "signalling", "immune", "pathway", "signalling", "signalling", "proliferation", "proliferation", "proliferation", "development","immune","metabolic", "metabolic", "metabolic", "metabolic", "pathway", "pathway", "DNA_damage","DNA_damage", "development", "metabolic", "immune", "signalling", "metabolic", "cellular_component", "immune", "development", "signalling", "signalling", "development"))


##' Below is a heatmap combining all diseases together for Hallmark I
##' normalised the pathway score by the number of signals per trait so
##' the heatmap is not dominated by traits with most signals
##' 
path <- 'h'
traits.path <- path.t[[path]]
## add number of signals per pathway per trait
traits.path <- merge(traits.path, sigpath[[path]], by=c("trait", "pathway"))

## Add column of normalised score, score/Nsig, 
traits.path[, norm_score := score/Nsig]

## Transform to wide, replace NAs with 0
traits.path.w <- dcast(traits.path, trait ~ pathway, value.var="norm_score", fill=0)

## Get distance between traits by using cosine distance to compare scores 

## Make an adjacency matrix connecting traits via pathway score
cos.d <- function(a,b) {
    num=sum(a*b)
    den=sqrt(sum(a^2)*sum(b^2))
    return(num/den)

}

traits.adj <- matrix(NA, nrow(traits.path.w), nrow(traits.path.w), dimnames=list(traits.path.w$trait, traits.path.w$trait))
for (x in 1:(nrow(traits.adj)-1)){   
    for (j in (x+1):nrow(traits.adj)){
        
        traits.adj[x,j] <- traits.adj[j,x] <- cos.d(traits.path.w[x,-"trait"], traits.path.w[j,-"trait"])
        
    }       
}
diag(traits.adj) <- 1


## Heatmap

## Cluster traits
hc.t <- hclust(as.dist(1-traits.adj), method="complete")
## Cluster pathways
path.cols <- colnames(d.paths[[path]])[colnames(d.paths[[path]]) %in% names(traits.path.w)]
adj.path <- 1-d.paths[[path]]
colnames(adj.path) <- rownames(adj.path) <- gsub("HALLMARK_", "", rownames(adj.path))

hc.path <- hclust(as.dist(adj.path[gsub("HALLMARK_", "",path.cols),  gsub("HALLMARK_", "",path.cols)]), method='complete')
## simplify trait names

path.mat <- as.matrix(traits.path.w[,path.cols, with=F ])
rownames(path.mat) <- sub("([A-Za-z0-9]+)_.*", "\\1", traits.path.w$trait)
colnames(path.mat) <- gsub("HALLMARK_", "", colnames(path.mat))

df.annot <- df.h[row.names(df.h) %in% path.cols,, drop=F]
rownames(df.annot) <- gsub("HALLMARK_", "",rownames(df.annot))
##+ fig.height=12, fig.width=12, fig.align="center"
pheatmap(path.mat, annotation_col=df.annot,
         show_rownames=T, show_colnames=T, silent=F, main="Clustering traits by pathways",
         cluster_cols=hc.path,cluster_rows=hc.t)



## Plot by trait, cluster by signal using 'hallmark' pathways
## annotation Hallmark pathways can be grouped by process
## category
## (https://www.cell.com/cell-systems/pdf/S2405-4712(15)00218-5.pdf)
hall <- pheatmap(d.paths[[1]], annotation_col=df.h, show_rownames=F, show_colnames=F, main="Hallmark pathway", silent =T)
hall.order <- colnames(d.paths[[1]])[hall$tree_row$order]



##' Similarity between pathways

grid.newpage()
grid.draw(hall)

##' #### Cluster signals by pathway distance
##'
##' Here I cluster signals by pathways

## Make matrix rows signals and cols pathways, replace NAs with 0

tfs.w <- lapply(tfs, function(i) {
    DT <- dcast(i, trait + sig ~ pathway, value.var = "score")
    for (j in 3:ncol(DT)){
        set(DT,which(is.na(DT[[j]])),j,0)
    }
    return(DT)

})


sig.mat <- lapply(tfs.w[1], function(i){
    ## go by trait
    setNames(lapply(unique(i$trait), function(u){
        ## calculate cosine distance between signals
        dt <- i[trait == u,]
        s <- unique(dt$sig)
        dt <- dt[, -c("trait","sig")]
        if(nrow(dt) <2) return(NULL)
        d <- matrix(NA, length(s), length(s), dimnames=list(s,s))
        for (x in 1:(length(s)-1)){   
            for (j in (x+1):length(s)){
                
                d[x,j] <- d[j,x] <- cos.d(dt[x,], dt[j,])
                
            }       
        }
        diag(d) <- 1
        return(d)
    }),
    unique(i$trait)
    )
})





cluster.path <- lapply(names(sig.mat), function(i) {
    ## go through traits
    lapply(names(sig.mat[[i]]), function(j) {

        mat <- sig.mat[[i]][[j]]
        if(is.matrix(mat)){
        dt <- tfs.w[[i]][trait == j,]
        t=sub("([A-Za-z0-9]+)_.*", "\\1", j)
        ## cluster by signals in each trait
        hc=hclust(as.dist(1-mat), method="complete")
        ## cluster pathways
        cols <- colnames(d.paths[[1]])[colnames(d.paths[[1]]) %in% names(dt)]
        hc.cols=hclust(as.dist(1-d.paths[[1]][cols,cols]), method='complete')

        df3 <- df.h[row.names(df.h) %in% cols,,drop=F]
        
        p <- pheatmap(dt[, 3:ncol(dt)],annotation_col=df3, show_rownames=F, show_colnames=F, silent=T, main=paste(t, "\n", "Signals:", nrow(dt)), cluster_cols=hc.cols,cluster_rows=hc)
        grid.draw(p)
        grid.newpage()
        }
    }
    )
})



##' ## Conclusions
##'
##' This is a very exploratory analysis as it doesnt include any
##' statistical test for enrichment.
##'
##' Across traits the most popular pathways were TNFa, IFNG, KRAS and
##' other immune mediated and metabolic.  Myogenesis and pancreas beta
##' cells were appear in some traits, but I need to check which are
##' the TFs driving those pathways.
##'
##' Clustering signals by shared pathways showed that in general few
##' signals concentrate most of the annotated pathways.
##'
##' Next, I will combine the annotated pathways with chromatin states
##' in immune cells.

    







