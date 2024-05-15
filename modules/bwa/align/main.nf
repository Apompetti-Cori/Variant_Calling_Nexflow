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
params.pubdir = "bwa_align"

process BWA_ALIGN {
    maxForks 4
    memory '40 GB'
    cpus 4
    
    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.id == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)
    val(state)

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
    def single_end = params.sample_table ? (meta.single_end ? true : false) : (params.single_end ? true : false)
    def file_id = params.sample_table ? meta.id : meta
    
    if(single_end){
        """
        bwa mem \
        -t ${task.cpus} \
        ${params.index} \
        ${reads} \
        -R "@RG\\tID:${file_id}\\tSM:${file_id}\\tLB:lib1\\tPL:illumina\\tPU:group1" > ${file_id}.sam
        """
    }
    else{
        """
        bwa mem \
        -t ${task.cpus} \
        ${params.index} \
        ${reads[0]} ${reads[1]} \
        -R "@RG\\tID:${file_id}\\tSM:${file_id}\\tLB:lib1\\tPL:illumina\\tPU:group1" > ${file_id}.sam
        """
    }
}