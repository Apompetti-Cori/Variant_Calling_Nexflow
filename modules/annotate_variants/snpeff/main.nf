#!/usr/bin/env nextflow

/*
Coriell Institute for Medical Research
Whole Exome Seq Pipeline. Started August 2023.

Contributors:
Anthony Pompetti <apompetti@coriell.org>
*/

/*
Enable Nextflow DSL2
*/
nextflow.enable.dsl=2

/*
Define local params 
*/
params.outdir = "./results"
params.pubdir = "annotate_variants"

process SNPEFF {
    maxForks 4
    memory '40 GB'
    cpus 4
    
    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*ann.g.vcf"), emit: vcf

    script:
    def file_id = params.sample_table ? meta.id : meta

    """
    snpEff -Xmx8g \
    ann -noStats -geneId \
    -filterInterval \
    ./scripts/loci.bed hg38 ${reads} | SnpSift  annotate ${params.annotations} > ${file_id}.ann.g.vcf
    """
}