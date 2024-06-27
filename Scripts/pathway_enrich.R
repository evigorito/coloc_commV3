#' ---
#' title: Assess significance for pathway enrichment based on predicted tf-binding site disruption by fine-mapped IMD SNPs
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
library(ggplot2)
library(cowplot)

##'  Inputs

args <- commandArgs(trailingOnly = TRUE)
fm <- args[2]
match <- args[3]
cs_tf <- grep("credset_fm", args, value=T)
match_tf <- grep("match_fm", args, value=T)

#' <!--
## fm='/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_tf/tf_trait/traits_tfs.rds'
## match='/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_match_snps/tf_trait/traits_matched_tfs.rds'
## cs_tf <- list.files('/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_traits_tf/tf_sites', full.names=T)
## match_tf <- list.files('/home/ev250/rds/rds-mrc-bsu/ev250/communities/fm_match_snps/tf_sites', full.names=T)
#'"-->
#'
##' ## Rationale
##'
##' In previous analysis I identified transcription factors predicted
##' to be disrupted by snps in credible sets. The I computed scores
##' adding posterior probabilities for transcription factors within
##' pathways. I have noticed that some pathways related to TNF
##' signalling, IFN among others tended to have the highest scores
##' across IMDs.
##'
##' Here, I implement an statistical test to see whether it is just by
##' chance or if there is anything interesting. The steps are as
##' follows:
##'
##' 1- For each snp in a credible set, find a matching snp
##' within the same block. The matching is performed by maf (5%
##' tolerance) and distance to proximal gene TSS.
##'
##' 2- Predict transcription binding site disruption with strong
##' effect using motifbreakR.
##'
##' 3- Using annotated pathways (HALLMARK, KEGG and REACTOME) group
##' transcription factors and compute an score by adding posterior
##' probabilities (PP) for all TFs within a pathway by credible set
##' and trait.
##'
##' 4- Use a paired t-test to assess mean difference in scores by
##' pathway and trait.
##'
##' 5- Perform silimar analysis without grouping tfs.

fm.path <- readRDS(fm)
fm.match <- readRDS(match)

##  Combine in long format
l <- mapply(function(a,b) {
    ## compare scores using t-test by trait and pathway
    rbind(a[, sig := paste(block, signal, sep="_")][,snp := "fm"], b[,snp:="matched"],fill=T)
    ## merge(a[, sig := paste(block, signal, sep="_")],b, by=c("trait", "sig", "pathway"), all=T, suffixes= c(".fm", ".match"))
},
a=fm.path,
b=fm.match,
SIMPLIFY=F)

## boxplots by trait and pathways
box.p <- lapply(l, function(i) {
    ## go by trait
    setNames(lapply(unique(i$trait), function(t) {
        ## go by pathway
        setNames(lapply(unique(i[trait == t,pathway]), function(p) {
            ggplot(i[trait == t & pathway == p,], aes(x=snp, y=score)) +
                geom_boxplot() +
                labs(title=t, subtitle=p)
        }), unique(i[trait == t,pathway]))
    }), unique(i$trait))

})

##  Combine in wide format
l.wide <- mapply(function(a,b) {
    m=merge(a[, sig := paste(block, signal, sep="_")],b, by=c("trait", "sig", "pathway"), all=T, suffixes= c(".fm", ".match"))
    ## replace NA with 0 for score.fm or score.match
    m[is.na(score.fm), score.fm := 0][is.na(score.match), score.match := 0,]
    return(m)

},
a=fm.path,
b=fm.match,
SIMPLIFY=F)

## Look at the distribution of observations (credible sets) across IMD and pathway

obs <- lapply(l.wide, function(i) i[,.N, .(trait, pathway)][, .(N)])

obs.p <- lapply(names(obs), function(i) ggplot(obs[[i]], aes(x=N)) +
                                        geom_histogram(color="black", fill="white", bins=20) +
                                        theme_bw() +
                                        labs(x="Number of obs", title=paste("Pathway",i)))

plot_grid(plotlist=obs.p, ncol=1)
 

## t-test
t.test <- lapply(l.wide, function(i) {
    ## go by trait
    rbindlist(lapply(unique(i$trait), function(t) {
        
        ## go by pathway
        rbindlist(lapply(unique(i[trait == t,pathway]), function(p) {
            
            ## check if enough observations, at least 20
            sub <- i[trait==t & pathway == p,]
            if(nrow(sub) >= 20){
                test=t.test(sub$score.fm, sub$score.match, paired=T)
                return(data.table(trait=t, pathway=p, obs=nrow(sub), pvalue=test$p.value))
            } else {
                return(sub[,.N, .(trait, pathway)])
            }
        }),fill=T)
    }), fill=T)
    
})

t.test.tf.p <- lapply(names(t.test), function(i) ggplot(t.test[[i]][!is.na(pvalue),], aes(x=pvalue))+
                                        geom_histogram(color="black", fill="white", bins=20) +
                                        theme_bw() +
                                        labs(x="pvalue", title=paste("Pathway",i)))
plot_grid(plotlist=t.test.tf.p, ncol=1)

## Number of tests
sapply(t.test, nrow)

## Distribution of p-values
lapply(t.test, function(i) summary(i$pvalue))


##' Try with individual tfs
##'
##' 


tfs_cs <- rbindlist(lapply(cs_tf, function(i) rbindlist(readRDS(i),fill=T)), fill=T)[effect == 'strong',]

tfs_cs[, sig:=paste(block, signal, sep="_")]

tfs_m <- rbindlist(lapply(match_tf, function(i) {
        dt <- readRDS(i)
        if("effect" %in% names(dt)){
            return(dt[effect == 'strong',])
        }
        }))

## Sum pp by tf, signal (cred set) and trait (IMD)

cs.score <- tfs_cs[, .(score=sum(pp)), by=.(geneSymbol, sig,trait)]
match.score <- tfs_m[,.(score=sum(pp)), by=.(geneSymbol, sig,trait)]

## merge by sig and run t.test

tf.score <- merge(cs.score, match.score, by=c("geneSymbol", "sig", "trait"), suffixes=c(".cs", ".match"), all=T)

## Replace NAs with 0s

for (j in seq_len(ncol(tf.score))) set(tf.score,which(is.na(tf.score[[j]])),j,0)

## Number of obs by trait and TF

ggplot(tf.score[, .N, .(trait, geneSymbol)] , aes(x=N)) +
    geom_histogram(color="black", fill="white", bins=50) +
    theme_bw() +
    labs(title="Observations per IMD and TF")

## Get p-values paired t-test

tf.t.test <- rbindlist(lapply(unique(tf.score$trait), function(t){
    ## go by TF if enough info
    rbindlist(lapply(unique(tf.score[trait==t, geneSymbol]), function(g){
        sub=tf.score[trait==t &geneSymbol == g,]
        ## 10 obs min
        if(nrow(sub) >= 10){
            test=t.test(sub$score.cs, sub$score.match, paired=T)
            return(data.table(trait=t, TF=g, obs=nrow(sub), pvalue=test$p.value))
        }
    }))
}))


ggplot(tf.t.test, aes(x=pvalue))+
    geom_histogram(color="black", fill="white", bins=20) +
    theme_bw() +
    labs(x="pvalue", title="Distribution of pvalues across TFs and IMDs")

## Number of tests
nrow(tf.t.test)

## Distribution of pvalues
summary(tf.t.test$pvalue)

##' ## Conclusions
##'
##' I didn't find much evidence of significant pathway enrichment
##' using any of the databases, for the number of tests performed.
##'
##' The distribution of p-values look pretty uniform.
##'
##' When tested significance by individual TFs, no evidence of
##' enrichment was found either.
