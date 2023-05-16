#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.loadJson = 'params.json'
params.compare = 'compare.csv'
params.GO_anno = 'go_gene.txt'
params.GO_description = 'kegg_description.txt'
params.KEGG_anno = 'kegg_gene.txt'
params.KEGG_description = 'kegg_description.txt'

process RNASEQSTAT2 { 
    time '1d'
    publishDir = [path: {"${params.outdir}/../2.RNASeqStat2"},
                  mode: params.publish_dir_mode]
    label "process_downstream"
    
    input:
    path 'rsem.merged.gene_counts.tsv'
    
    output:
    path "*_*_*_*"        
    path "Rplots.pdf"
    path "stat_2_txt.csv"
    path "*_vs_*"    ,emit: stat  
    
    script:
    """
    Rscript ${projectDir}/bin/rnaseqstat2.r \
    ${params.compare} \
    ./rsem.merged.gene_counts.tsv \
    ${params.GO_anno} \
    ${params.GO_description} \
    ${params.KEGG_anno} \
    ${params.KEGG_description}
    python3 ${projectDir}/bin/stat2txt.py
    python3 ${projectDir}/bin/stat3csv.py
    python3 ${projectDir}/bin/stat4lollipop.py
    """
}
