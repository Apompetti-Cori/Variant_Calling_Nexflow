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

//Configurable variables for pipeline
params.fastq_folder = "${workflow.projectDir}/fastq"
params.reads = "${params.fastq_folder}/*{_L00,}{1,2,3,4}{_R,_}{1,2}*.{fastq,fq}.gz"
params.single_end = false
params.multi_lane = false
params.genome = false
params.sample_table = false
params.db = params.genomes ? params.genomes[ params.genome ].db ?:false : false
params.index = params.genomes ? params.genomes[ params.genome ].index ?:false : false
params.annotations = params.genomes ? params.genomes[ params.genome ].annotations ?:false : false
params.funcotator = params.genomes ? params.genomes[ params.genome ].funcotator ?:false : false

//Include modules to main pipeline
include { PREPROCESS } from './modules/preprocess/main'
include { FASTQC as PRETRIM_FASTQC } from './modules/fastqc/main' addParams(pubdir: 'pretrim_fastqc')
include { TRIM_GALORE } from './modules/trim_galore/main'
include { FASTQC as POSTTRIM_FASTQC } from './modules/fastqc/main' addParams(pubdir: 'posttrim_fastqc')
include { BWA_ALIGN } from './modules/bwa/align/main'
include { SAMTOOLS_SORT } from './modules/samtools/sort/main'
include { SAMTOOLS_INDEX } from './modules/samtools/index/main'
include { CALL_VARIANTS } from './modules/gatk/call_variants/main'
include { MUTECT } from './modules/gatk/mutect/main'
include { FILTER_MUTECT } from './modules/gatk/filter_mutect/main'
include { FUNCOTATOR } from './modules/annotate_variants/funcotator/main'
include { SNPEFF } from './modules/annotate_variants/snpeff/main'
include { MULTIQC } from './modules/multiqc/main'

//Provide sample table in csv format to have pipeline process samples via sample table
if ( params.sample_table ){
    // Channel for the samplesheet
    ch_samplesheet = Channel.fromPath(params.sample_table)
    
    // Parse it line by line
    reads_ch = ch_samplesheet.splitCsv(header:true).filter{
        // Filter out observations without files
        it['r1_L1'] != ''
    }.map{

        // This is the read1 and read2 entry
        r1_L1 = it['r1_L1']
        r1_L2 = it['r1_L2']
        r1_L3 = it['r1_L3']
        r1_L4 = it['r1_L4']
        r2_L1 = it['r2_L1']
        r2_L2 = it['r2_L2']
        r2_L3 = it['r2_L3']
        r2_L4 = it['r2_L4']

        // Detect wiether single-end or paired-end
        is_singleEnd = r2_L1.toString()=='' ? true : false
        is_multiLane = r1_L2.toString()=='' ? false : true
        is_emptyObsv = r1_L1.toString()=='' ? true : false
        
        // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and a single_end logical entry
        meta = [id: it['sample'], batch: it['batch'], single_end: is_singleEnd, multi_lane: is_multiLane, empty_obsv: is_emptyObsv]
        
        // We return a nested map, the first entry is the meta map, the second one is the read(s)
        if ( is_singleEnd ){
            if ( r1_L4.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3, r1_L4], []]
            }
            else if ( r1_L3.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3], []]
            }
            else if ( r1_L2.toString()!='' ){
                [meta, [r1_L1, r1_L2], []]
            }
            else {
                [meta, [r1_L1], []]
            }
        }
        else{
            if (r2_L4.toString()!=''){
                [meta, [r1_L1, r1_L2, r1_L3, r1_L4], [r2_L1, r2_L2, r2_L3, r2_L4]]
            }
            else if ( r2_L3.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3], [r2_L1, r2_L2, r2_L3]]
            }
            else if ( r2_L2.toString()!='' ){
                [meta, [r1_L1, r1_L2], [r2_L1, r2_L2]]
            }
            else {
                [meta, [r1_L1], [r2_L1]]
            }
        }
    }
}
else {
    exit 1, "No sample table provided please submit path of sample table using --sample_table flag"
}

// Create channel for fasta, index, and dictionary
db_ch = Channel.fromPath(params.db.fasta).concat(Channel.fromPath(params.db.index), Channel.fromPath(params.db.dict)).collect(flat: false)

workflow somatic_variants {
    
    //PREPROCESS the sample table to change the files listed inside. Concatenates any multilane files.
    PREPROCESS(reads_ch)
    reads_ch = PREPROCESS.out

    //Run trim_galore on raw reads
    TRIM_GALORE(reads_ch)

    //Run fastqc on trimmed reads
    POSTTRIM_FASTQC(
        TRIM_GALORE.out.reads
    )
    
    //Run bismark_align on trimmed reads
    //State Dependency: Wait until TRIM_GALORE is done to run 
    //State Dependency: Wait until POSTTRIM_FASTQC is done to run 
    state = POSTTRIM_FASTQC.out.collect()
    BWA_ALIGN(
        TRIM_GALORE.out.reads.collect(flat: false).flatMap(),
        state
    )
    
    //Run samtools_sort on aligned reads
    SAMTOOLS_SORT(
        BWA_ALIGN.out.sam.collect(flat: false).flatMap()
    )
    
    //Run gatk Mutect on sorted reads
    MUTECT(
      SAMTOOLS_SORT.out.bam,
      db_ch
    )
    
    //Run gatk FilterMutectCalls on mutect reads
    FILTER_MUTECT(
      MUTECT.out.vcf,
      db_ch
    )
    
    //Run snpEff on vcf file
    //ANNOTATE_VARIANTS(CALL_VARIANTS.out.vcf.collect())

    //Run Funcotator on vcf file
    //FUNCOTATOR(
    //    FILTER_MUTECT.out.vcf,
    //    params.funcotator.source
    //)

    
    //Run multiqc on pretrim fastqc output, trim_galore trimming report, posttrim fastqc output, bismark conversion output
    //MULTIQC(PRETRIM_FASTQC.out.collect().combine(POSTTRIM_FASTQC.out.collect()).combine(TRIM_GALORE.out.report.collect()))
}


workflow {
    somatic_variants()
}
