#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Bulk RNA-seq Analysis Pipeline
========================================================================================
*/

// Parameter validation
def validateParameters() {
    if (!params.input) {
        error "Input samplesheet not specified with e.g. '--input samplesheets/GSE135251_samplesheet.csv'"
    }
    if (!params.genome) {
        error "Reference genome not specified with '--genome'"
    }
    if (!params.gtf) {
        error "GTF annotation file not specified with '--gtf'"
    }
}

// Include modules and subworkflows
include { FASTQC } from './modules/fastqc'
include { TRIMGALORE } from './modules/trimgalore'
include { STAR_ALIGN } from './modules/star'
include { FEATURECOUNTS } from './modules/featurecounts'
include { DESEQ2 } from './modules/deseq2'
include { MULTIQC } from './modules/multiqc'

// Main workflow
workflow {
    // Validate parameters
    validateParameters()
    
    // Read samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> tuple(
            row.sample,
            file(row.fastq_1),
            file(row.fastq_2),
            row.group_in_paper,
            row.disease,
            row.Fibrosis_stage,
            row.nas_score,
            row.Stage
        )}
        .set { input_samples }
    
    // Quality control
    FASTQC(input_samples)
    
    // Trimming if enabled
    if (params.trim) {
        TRIMGALORE(input_samples)
        reads_for_alignment = TRIMGALORE.out.trimmed_reads
    } else {
        reads_for_alignment = input_samples
    }
    
    // Alignment
    STAR_ALIGN(
        reads_for_alignment,
        file(params.genome),
        file(params.gtf),
        params.star_index
    )
    
    // Quantification
    FEATURECOUNTS(
        STAR_ALIGN.out.bam,
        file(params.gtf)
    )
    
    // Differential expression
    DESEQ2(
        FEATURECOUNTS.out.counts,
        input_samples.map { it -> tuple(it[0], it[3], it[4], it[5], it[6], it[7]) }
    )
    
    // Quality report
    MULTIQC(
        Channel.empty().mix(
            FASTQC.out.zip,
            TRIMGALORE.out.reports,
            STAR_ALIGN.out.logs,
            FEATURECOUNTS.out.summary
        ).collect()
    )
}

// Workflow completion notification
workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    """
} 