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
params.pubdir = "trim_galore"

/*
Run trim_galore on each read stored within the reads_ch channel
*/
process TRIM_GALORE {
    maxForks 4
    memory '8 GB'
    cpus 4

    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.id == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*{_val_1,_val_2,_trimmed}.fq.gz"), emit: reads
    path("*trimming_report.txt"), emit: report

    script:
    def single_end = params.sample_table ? (meta.single_end ? '' : '--paired') : (params.single_end ? '' : '--paired')
    def file_id = params.sample_table ? meta.id : meta
    
    """
    trim_galore \
    ${single_end} \
    --retain_unpaired \
    --illumina \
    --quality 30 \
    --phred33 \
    --clip_R1 25 --clip_R2 25 \
    --cores ${task.cpus} \
    ${reads}
    """
}