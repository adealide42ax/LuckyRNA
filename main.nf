#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.loadJson      = 'params.json'
params.compare       = 'compare.csv'
params.GO_anno       = 'go_gene.txt'
params.GO_description = 'go_description.txt'
params.KEGG_anno     = 'kegg_gene.txt'
params.KEGG_description = 'kegg_description.txt'

params.fasta         = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf           = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff           = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed      = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.bbsplit_index = WorkflowMain.getGenomeAttribute(params, 'bbsplit')
params.star_index    = WorkflowMain.getGenomeAttribute(params, 'star')
params.hisat2_index  = WorkflowMain.getGenomeAttribute(params, 'hisat2')
params.rsem_index    = WorkflowMain.getGenomeAttribute(params, 'rsem')
params.salmon_index  = WorkflowMain.getGenomeAttribute(params, 'salmon')

WorkflowMain.initialise(workflow, params, log)

include { RNASEQ } from './workflows/rnaseq'
include { DATA_CLEAN } from './workflows/data_clean'
include { RNASEQSTAT2 } from './workflows/rnaseqstat'
include { WGCNA } from './workflows/wgcna'
include { STATISTICS } from './workflows/lollipop'

workflow NFCORE_RNASEQ {
    RNASEQ ()
    DATA_CLEAN(params.outdir)
    WGCNA(DATA_CLEAN.out.tpm)
    RNASEQSTAT2(DATA_CLEAN.out.counts)
    STATISTICS(RNASEQSTAT2.out.stat)
}

workflow 
{
    NFCORE_RNASEQ ()
}
