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
params.pubdir = "multiqc"
params.multiqc_config = "${workflow.projectDir}/multiqc_config.yaml"
params.multiqc_report_title = "MultiQC Report"

process MULTIQC {
    memory '32 GB'
    cpus 1

    conda '/opt/miniconda3/envs/multiqc'

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'
    
    input:
    path('multiqc_input/*')

    output:
    path("*.html")

    script:
    """
    multiqc multiqc_input/ --config ${params.multiqc_config} --title ${params.multiqc_report_title}
    """
}