2###################################################################################
## Snakefile for clustering coloc signals V3
# Here we start again, this time with a descriptive analysis of the IMD fine mapping
# and improving in the communities clustering as most communities were difficult to
# define when multiple signals in the same region. Chris is looking into multi-coloc
# to coloc all traits together instead of 2by2 comparison.
###################################################################################


shell.prefix("source ~/.bashrc; ")

configfile: "../../config.yaml"

localrules: all

import numpy as np
import pandas as pd
import os

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

# Functions

def select_samples_chromstates(File, sep="\t", search_values = ['HSC & B-cell', 'Blood & T-cell'], exclude=['BSS00313']):
    """Selects values in column GROUP and creates a dictionary with keys unique values for infoline column and values id column. The default File is the metadata file. """
    
    data=pd.read_csv(File, sep=sep)
    data=data[data.GROUP.str.contains('|'.join(search_values))]
    dic={}

    for k in data.infoline.unique():
        
        dic[k]=[ x for x in data.loc[data['infoline'] == k, 'id'].tolist() if x not in exclude]
    return(dic)

# variables and wildcards
dic_18states = select_samples_chromstates(File=config['communities'] + "/chrom18states/main_metadata_table.tsv")

chrom_samples_id = [item for l in list(dic_18states.values()) for item in l]
# fine mapped basename removing ext
fm_blocks = [os.path.splitext(os.path.basename(x))[0] for x in os.listdir(config['communities'] + "/fine-mapping")]

blocks=[x.replace("fm_", "") for x in fm_blocks]  


paths=["h.all.v2023.1.Hs.symbols.gmt", "c2.cp.kegg.v2023.1.Hs.symbols.gmt", "c2.cp.reactome.v2023.1.Hs.symbols.gmt", "c2.cp.biocarta.v2023.1.Hs.symbols.gmt", "c2.cp.pid.v2023.1.Hs.symbols.gmt","c2.cp.wikipathways.v2023.1.Hs.symbols.gmt","c7.immunesigdb.v2023.1.Hs.symbols.gmt"]

chroms=[x+1 for x in range(22)]

rule all:
    input:
        # expand(config['out_dir'] + "/chrom18states_b38/{id}_18_CALLS_segments.bed.gz", id=chrom_samples_id),
        # expand(config['out_dir'] + "/fm_traits_chrom_stats_b38/{fm_block}_chromstates.rds", fm_block=fm_blocks),
        # config['out_dir'] + "/fm_traits_chrom_stats_b38/traits_chromstates.rds",
        # config['out_dir'] + "/fm_traits_chrom_stats_b38/traits_prob_chromstates.rds",
        # config['out_dir'] + "/chrom18states_b38_aggregated/summary_chrom_states_sample.txt",
        # "Scripts/compare_chromstates.html"
        # expand(config['out_dir'] + "/fm_traits_tf/tf_sites/tf_prediction_credset_{fm_block}.rds", fm_block=fm_blocks),
        # config['out_dir'] + "/fm_traits_tf/tf_entrez/entrez_traits_tfs.txt",
        # config['out_dir'] + "/fm_traits_tf/tf_trait/traits_tfs.rds",
        # "Scripts/compare_pathways.html",
        #  "Scripts/tfs_chromstates.html",
        # config['out_dir'] + "/V3/time2run/map_tf2credset.txt"
        # expand(config['out_dir'] + "/fm_match_snps/allel_freq/ALL.chr{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.txt", chrom=chroms)
        # expand(config['out_dir'] + "/fm_match_snps/matched_snps/match4credset_fm_{block}.txt",block=blocks)
        # expand(config['out_dir'] + "/fm_match_snps/tf_sites/tf_prediction_match_{fm_block}.rds", fm_block=fm_blocks)
        # config['out_dir'] + "/fm_match_snps/tf_trait/traits_matched_tfs.rds",
        "Scripts/pathway_enrich.html"
        
       

rule down_samples_chromstates_b38:
    """ Download samples for immune cell types in built38"""
    input:
        HTTP.remote("https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/{id}_18_CALLS_segments.bed.gz", keep_local=True)
    output:
        config['out_dir'] + "/chrom18states_b38/{id}_18_CALLS_segments.bed.gz"
    run:
        outputName=os.path.join(config['out_dir'] + "/chrom18states_b38/", os.path.basename(input[0]))
        shell("mv {input} {outputName} ")

        # CHECK PIPELINE, I MADE A MISTAKE, I EXCLUDED BLOCKS WITH 10 OR MORE TRAITS, INSTEAD OF EXCLUDING A TRAIT WITHIN A BLOCK WITH 10 OR MORE SIGNALS (APPLIES TO CHROM_STATES AND DOWNSTREAM HTML FOR TFS). CORRECTED CODE IN match_credset.R (18/6)
        
rule fm_credset_chromstate:
    """Compute probabilities of active, inactive and quiescent chromatine state integrating info across 95% credible sets for fine mapped IMDs. Output is a list with first element a table with traits per block and second element a list with each element a trait with it matched chromstates."""
    input:
        chromstates=expand(config['out_dir'] + "/chrom18states_b38/{id}_18_CALLS_segments.bed.gz", id=chrom_samples_id),
        fm=config['communities'] + "/fine-mapping/{fm_block}.RData"
    params:
        credsize=0.95
    output:
        out=config['out_dir'] + "/fm_traits_chrom_stats_b38/{fm_block}_chromstates.rds"
    script:
        "Scripts/fm_credset_chromstate.R"

        
rule chromstate_by_trait:
    """In this rule I combine the chromstate annotations by trait. In previous rule where combined by block."""
    input:
        chrom_credsets=expand(config['out_dir'] + "/fm_traits_chrom_stats_b38/{fm_block}_chromstates.rds", fm_block=fm_blocks)
    output:
        out=config['out_dir'] + "/fm_traits_chrom_stats_b38/traits_chromstates.rds"
    script:
        "Scripts/traits_chromstates.R"


rule summarise_chromstates:
    """Here I summarise within an between chromatin annotations by trait"""
    input:
        chrom_trait=config['out_dir'] + "/fm_traits_chrom_stats_b38/traits_chromstates.rds"
    output:
        out=config['out_dir'] + "/fm_traits_chrom_stats_b38/traits_prob_chromstates.rds"
    script:
        "Scripts/summrise_chromstates.R"

rule aggregate_chrom_states_in_immune_samples:
    """Here I calculate the proportion of active/quiescence/inactive marks genome wide across the epimap samples. This will be useful to compare with the profile obtained when subsetting the GWAS fine-maaped SNPs"""
    input:
        chromstates=expand(config['out_dir'] + "/chrom18states_b38/{id}_18_CALLS_segments.bed.gz", id=chrom_samples_id),
    output:
        out=config['out_dir'] + "/chrom18states_b38_aggregated/summary_chrom_states_sample.txt"
    script:
        "Scripts/aggregate_chrom_states.R"
    
rule compare_chromstates:
    """Using the aggregated chromatin states per trait I look for clusters of chromstates samples within traits"""
    input:
        inp=config['out_dir'] + "/fm_traits_chrom_stats_b38/traits_prob_chromstates.rds",
        meta=config['out_dir'] + "/chrom18states/main_metadata_table.tsv",
        chrom_agg=config['out_dir'] + "/chrom18states_b38_aggregated/summary_chrom_states_sample.txt",
        script="Scripts/compare_chromstates.R"
    output:
        out="Scripts/compare_chromstates.html"
    script:
        "../../Scripts/RenderReport.R"

rule map_tf2credset:
    """Here I annotate fine mapped SNPs to the most likely TF binding site disruption. Based on rule tf_annot from communitiesV2"""
    input:
        fm=config['communities'] + "/fine-mapping/{fm_block}.RData",
        r=config['r43'],
        script="Scripts/tf_credset.R"
    resources:
        tmpdir=config['out_dir']
    threads:
        4
    params:
        credsize=0.95
    output:
        out=config['out_dir'] + "/fm_traits_tf/tf_sites/tf_prediction_credset_{fm_block}.rds",
    shell:
        "conda activate r43;"
        "{input.r}  {input.script} {input.fm} {params.credsize} {output.out} "

rule tf_symbol:
    """Select TFs with binding sites predicted to be disrupted by fm SNPs with strong effect and get entrez id"""
    input:
        fm=expand(config['out_dir'] + "/fm_traits_tf/tf_sites/tf_prediction_credset_{fm_block}.rds", fm_block=fm_blocks),
        r=config['r43'],
        script="Scripts/tf_entrez.R"
    params:
        inp_dir=config['out_dir'] + "/fm_traits_tf/tf_sites"
    output:
        out=config['out_dir'] + "/fm_traits_tf/tf_entrez/entrez_traits_tfs.txt"
    shell:
        "conda activate r43;"
        "{input.r}  {input.script} {params.inp_dir} {output.out} "

        
rule combine_tf:
    """Here I group predicted TF binding sites across IMDs and look for pathways within disease"""
    input:
        entrez=config['out_dir'] + "/fm_traits_tf/tf_entrez/entrez_traits_tfs.txt",
        paths=expand(config['out_dir'] + "/objects/pathways/{path}", path=paths[0:3]),
        fm=expand(config['out_dir'] + "/fm_traits_tf/tf_sites/tf_prediction_credset_{fm_block}.rds", fm_block=fm_blocks),
        r=config['r43'],
        script="Scripts/traits_tfs.R"
    params:
        n=["h", "kegg", "reactome"],
    output:
        out=config['out_dir'] + "/fm_traits_tf/tf_trait/traits_tfs.rds"
    shell:
        "conda activate r43;"
        "{input.r}  {input.script} {input.entrez} {input.paths} {params.n}  "
        "{output.out} {input.fm} "
        

rule compare_pathways:
    """Similar to rule compare_chromstates, I look if pathways can cluster IMDs"""
    input:
        inp=config['out_dir'] + "/fm_traits_tf/tf_trait/traits_tfs.rds",
        paths=expand(config['out_dir'] + "/objects/pathways/{path}", path=paths[0:3]),
        r=config['r43'],
        script2call="Scripts/compare_pathways.R",
        script="Scripts/RenderReportconda.R"
    output:
        out="Scripts/compare_pathways.html"
    shell:
        "conda activate r43;"
        "{input.r} {input.script} {input.script2call} {input.inp} "
        " {input.paths} "

rule time2run_tfs:
    """The analysis in compare_pathways does not have an statistical test. Would be good to design one. To do that I first calculate how much time it takes to run all steps in the pipeline in order to design a permutation test. Here I extract running time from log files"""
    input:
        # map2tf2credset=expand(config['out_dir'] + "/V3/logs/map_tf2credset.fm_block\={fm_block}.out", fm.block=fm_blocks),
        r=config['r43'],
        script="Scripts/run_time.R"
    params:
        pattern=config['out_dir'] + "/V3/logs/map_tf2credset*"
    output:
        out=config['out_dir'] + "/V3/time2run/map_tf2credset.txt",
    shell:
        "grep 2024 {params.pattern} | awk -F':' '{{print $1,$2,$3,$4}}' > {output.out}; "
        "conda activate r43;"
        "{input.r} {input.script} {output.out} "
        

rule get_maf_from_vcf:
    """Here I use 1000G built38 to extract allele frequency info. At the end I didnt use these files. I have matched using snps within the same block"""
    input:
        bcf=config["ref38"] + "/ALL.chr{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
    output:
        out=config['out_dir'] + "/fm_match_snps/allel_freq/ALL.chr{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.txt"
    shell:
        "bcftools query -f '%CHROM %POS %REF %ALT %AF\\n' {input.bcf} > {output.out} "
        
rule maf_match:
    """To design stat test I need to match snps in cred set with snps within the same block by maf and if possible by distance to TSS. I used the files that Chris generated containing MAF for SNPs wihin a fine mapping block extracted from EUR from UKBB. Then I also match for distance to nearest gene TSS. When a snp is in cred sest for many triats, I match that snp with another snp for all traits."""
    input:
        fm=config['communities'] + "/fine-mapping/fm_{block}.RData",
        maf=config['coloc'] + "/tmp_mafld_{block}.RData",
        r=config['r43'],
        script="Scripts/match_credset.R"
    params:
        credsize=0.95
    output:
        out=config['out_dir'] + "/fm_match_snps/matched_snps/match4credset_fm_{block}.txt",
    shell:
        "conda activate r43;"
        "{input.r}  {input.script} {input.fm} {input.maf} {params.credsize} {output.out} "
        
rule map_tf2match:
    """Here I annotate fine mapped SNPs to the most likely TF binding site disruption. Based on rule tf_annot from communitiesV2"""
    input:
        match_fm=config['out_dir'] + "/fm_match_snps/matched_snps/match4credset_{fm_block}.txt",
        r=config['r43'],
        script="Scripts/tf_match.R"
    resources:
        tmpdir=config['out_dir']
    threads:
        4
    output:
        out=config['out_dir'] + "/fm_match_snps/tf_sites/tf_prediction_match_{fm_block}.rds",
    shell:
        "conda activate r43;"
        "{input.r}  {input.script} {input.match_fm}  {output.out} "
        
rule combine_match_tf:
    """Here I group predicted TF binding sites across matched snps for fine mapped IMDs and look for pathways within disease"""
    input:
        paths=expand(config['out_dir'] + "/objects/pathways/{path}", path=paths[0:3]),
        tf_match=expand(config['out_dir'] + "/fm_match_snps/tf_sites/tf_prediction_match_{fm_block}.rds", fm_block=fm_blocks),
        r=config['r43'],
        script="Scripts/traits_matched_tfs.R"
    params:
        n=["h", "kegg", "reactome"],
    output:
        out=config['out_dir'] + "/fm_match_snps/tf_trait/traits_matched_tfs.rds"
    shell:
        "conda activate r43;"
        "{input.r}  {input.script} {input.paths} {params.n}  "
        "{output.out} {input.tf_match} "

rule assess_pathway_enrichment:
    """Here I test for significant enrichment of a pathway within an IMD comparing scores of fm-cred sets vs matched snps, both for pathways or individual tfs"""
    input:
        fm=config['out_dir'] + "/fm_traits_tf/tf_trait/traits_tfs.rds",
        control=config['out_dir'] + "/fm_match_snps/tf_trait/traits_matched_tfs.rds",
        match_tf=expand(config['out_dir'] + "/fm_match_snps/tf_sites/tf_prediction_match_{fm_block}.rds",
        fm_block =fm_blocks),
        cs_tf=expand(config['out_dir'] + "/fm_traits_tf/tf_sites/tf_prediction_credset_{fm_block}.rds", fm_block=fm_blocks),
        r=config['r43'],
        script2call="Scripts/pathway_enrich.R",
        script="Scripts/RenderReportconda.R"
        # script="Scripts/pathway_enrich.R"
    output:
        out="Scripts/pathway_enrich.html"
    shell:
        "conda activate r43;"
        "{input.r} {input.script} "
        "{input.script2call} "
        "{input.fm} {input.control} "
        "{input.match_tf} {input.cs_tf} "
        
    



        
rule combine_chromstates_tfs:
    """Here I look whether the intersect of active regions and TF pathways"""
    input:
        inp=config['out_dir'] + "/fm_traits_chrom_stats_b38/traits_prob_chromstates.rds",
        meta=config['out_dir'] + "/chrom18states/main_metadata_table.tsv",
        chrom_agg=config['out_dir'] + "/chrom18states_b38_aggregated/summary_chrom_states_sample.txt",
        tfs=config['out_dir'] + "/fm_traits_tf/tf_trait/traits_tfs.rds",
        paths=expand(config['out_dir'] + "/objects/pathways/{path}", path=paths[0:3]),
        qcblocks=expand(config['out_dir'] + "/fm_traits_chrom_stats_b38/{fm_block}_chromstates.rds", fm_block=fm_blocks),
        r=config['r43'],
        script2call="Scripts/tfs_chromstates.R",
        script="Scripts/RenderReportconda.R"
        # script="Scripts/tfs_chromstates.R"
    output:
        out="Scripts/tfs_chromstates.html"
    shell:
        "conda activate r43;"
        "{input.r} {input.script} "
        "{input.script2call} "
        "{input.inp} {input.meta} "
        "{input.chrom_agg} {input.tfs} {input.paths} {input.qcblocks} "
        
                

 
       
# snakemake  -k -j 500 --cluster-config cpu.json --cluster "sbatch -A {cluster.account} -p {cluster.partition}  -c {cluster.cpus-per-task}   -t {cluster.time} --output {cluster.error} -J {cluster.job} "
