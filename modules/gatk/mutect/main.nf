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
params.pubdir = "mutect/raw"
params.chip_wb = "${workflow.projectDir}/tables/interval_chip_wb.bed"

process MUTECT {
    maxForks 4
    memory '40 GB'
    cpus 4
    
    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.id == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)
    tuple path(fasta), path(index), path(dict)

    output:
    tuple val(meta), path("*.raw.vcf.gz"), path("*.gz.tbi"), path("*.gz.stats"), emit: vcf
    
    shell:
    '''
    gatk Mutect2 \
    -R !{fasta} \
    -I !{reads} \
    --max-reads-per-alignment-start 10000 \
    -L !{params.chip_wb} \
    -O !{meta.id}.raw.vcf.gz
    '''
}