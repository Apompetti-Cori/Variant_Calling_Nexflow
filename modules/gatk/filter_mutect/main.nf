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
params.pubdir = "mutect/filtered"
params.chip_wb = "${workflow.projectDir}/tables/interval_chip_wb.bed"

process FILTER_MUTECT {
    maxForks 4
    memory '40 GB'
    cpus 4
    
    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.id == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(vcf_stats)
    tuple path(fasta), path(index), path(dict)

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf
    path("*.filteringStats.tsv"), emit: stats
    
    shell:
    '''
    gatk FilterMutectCalls \
    -R !{fasta} \
    -V !{vcf} \
    -O !{meta.id}.filtered.vcf.gz
    '''
}