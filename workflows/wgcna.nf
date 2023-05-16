#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process WGCNA 
{ 
    publishDir = [path: {"${params.outdir}/../3.WGCNA/"},
                  mode: params.publish_dir_mode]
    label "wgcna"
    
    input:
    path "rsem.merged.gene_tpm.tsv"
    
    output:
    path '*'

    script:
    """
    mv ./rsem.merged.gene_tpm.tsv ./rsem.merged.gene_tpm.txt
    Rscript ${projectDir}/bin/wgcna.r ./rsem.merged.gene_tpm.txt
    """
    
}

/*
workflow {
    WGCNA ()
}
*/