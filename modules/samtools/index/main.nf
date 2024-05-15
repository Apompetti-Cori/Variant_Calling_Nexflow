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
params.pubdir = "samtools_index"

process SAMTOOLS_INDEX {
    maxForks 4
    memory '40 GB'
    
    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path("*.bam.bai"), emit: bam_index

    script:
    def file_id = params.sample_table ? meta.id : meta
    
    """
    samtools index \
    ${reads} \
    "${file_id}.bam.bai"
    """
}