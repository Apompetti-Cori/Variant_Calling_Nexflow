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

process PREPROCESS {
    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("*{_mL,_sL,_sL_1,_sL_2,_mL_1,_mL_2}.fq.gz")

    script:
    if ( meta.multi_lane ){
        if (meta.single_end){
            """
            cat ${reads_1} > ${meta.id}_mL.fq.gz
            """
        }
        else {
            """
            cat ${reads_1} > ${meta.id}_mL_1.fq.gz
            cat ${reads_2} > ${meta.id}_mL_2.fq.gz
            """
        }
    }
    else {
        if (meta.single_end){
            """
            mv ${reads_1} ${meta.id}_sL.fq.gz
            """
        }
        else {
            """
            mv ${reads_1} ${meta.id}_sL_1.fq.gz
            mv ${reads_2} ${meta.id}_sL_2.fq.gz
            """
        }
    }
}