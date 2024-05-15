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
params.pubdir = "funcotator"

process FUNCOTATOR {
    maxForks 4
    memory '40 GB'
    cpus 4
    
    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path(funcotator_source)
    
    output:
    tuple val(meta), path("*.funcotated.maf"), emit: maf

    script:
    def file_id = params.sample_table ? meta.id : meta

    """
    gatk t \
    -I ${reads}
    
    gatk Funcotator \
    --variant ${reads} \
    --reference ${funcotator.db} \
    --ref-version ${params.genome} \
    --data-sources-path ${funcotator_source} \
    --output ${file_id}.funcotated.maf \
    --output-file-format MAF
    """
}