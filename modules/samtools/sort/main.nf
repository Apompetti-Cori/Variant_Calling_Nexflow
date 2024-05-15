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
params.pubdir = "samtools_sort"

process SAMTOOLS_SORT {
    maxForks 4
    memory '40 GB'
    cpus 4
    
    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.id == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path("*.bam.bai"), emit: bam_index

    script:
    def file_id = params.sample_table ? meta.id : meta
    
    """
    samtools sort \
    --threads ${task.cpus} \
    ${reads} \
    -O BAM \
    -o "${file_id}.bam"
    
    samtools index \
    ${file_id}.bam \
    "${file_id}.bam.bai"
    """

}