#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LOLLIPOP 
{ 
    publishDir = [path: {"${params.outdir}/../2.RNASeqStat2"},
                  mode: params.publish_dir_mode]
    label "lollipop"
    
    input:
    path "*"

    output:
    path "*/*.pdf"

    script:
    """
    Rscript ${projectDir}/bin/lollipop.r
    """
    
}

/*
workflow {
    (params.outdir)
}
*/
