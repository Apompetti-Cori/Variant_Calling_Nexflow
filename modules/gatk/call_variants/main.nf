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
params.pubdir = "call_variants"
params.chip_wb = "${workflow.projectDir}/tables/interval_chip_wb.bed"

process CALL_VARIANTS {
    maxForks 4
    memory '40 GB'
    cpus 4
    
    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path(reads_index)
    tuple path(fasta), path(index), path(dict)

    output:
    tuple val(meta), path("*.g.vcf.gz"), emit: vcf
    
    shell:
    '''
    gatk \
    --java-options "-Xmx8g" HaplotypeCaller \
    --max-reads-per-alignment-start 10000 \
    -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
    -L !{params.chip_wb} \
    --reference !{fasta} \
    --sample-name !{meta.id} \
    --input !{reads} --output !{meta.id}.g.vcf.gz \
    -ERC GVCF
    '''
}