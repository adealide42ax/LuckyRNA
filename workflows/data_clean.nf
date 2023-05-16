#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = './rnaseq_result'

process DATA_CLEAN 
{ 

    label "data_clean"
    
    input:
    path outdir
    
    output:
    path 'rsem.merged.gene_counts.tsv' ,emit: counts
    path 'rsem.merged.gene_tpm.tsv' ,emit: tpm

    when :
    file("${params.outdir}/star_rsem/rsem.merged.gene_counts.tsv").exists()&&file("${params.outdir}/star_rsem/rsem.merged.gene_tpm.tsv").exists()
    
    script:
    """
    cp ${params.outdir}/star_rsem/rsem.merged.gene_counts.tsv ./
    cp ${params.outdir}/star_rsem/rsem.merged.gene_tpm.tsv ./
    python3 ${projectDir}/bin/del_transcript_id.py ./rsem.merged.gene_counts.tsv
    python3 ${projectDir}/bin/del_transcript_id.py ./rsem.merged.gene_tpm.tsv
    """
    
}

/*
workflow 
{
    DATA_CLEAN(outdir)

}
*/