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

// Load modules
include { FASTQC } from './modules/fastqc'
include { STAR_ALIGN } from './modules/star'
include { FEATURECOUNTS } from './modules/featurecounts'
include { DESEQ2 } from './modules/deseq2'
include { MULTIQC } from './modules/multiqc'

// Main workflow
workflow {
    // Validate parameters
    validateParameters()
    
    // Create channels for input samples from the sample manifest
    ch_samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample
            def group = row.group_in_paper
            def disease = row.disease
            def fibrosis = row.Fibrosis_stage
            def nas = row.nas_score
            def stage = row.Stage
            def r1 = file(row.fastq_1)
            def r2 = file(row.fastq_2)
            return [sample_id, group, disease, fibrosis, nas, stage, r1, r2]
        }

    // Run FASTQC on raw reads
    FASTQC(
        ch_samples.map { sample_id, group, disease, fibrosis, nas, stage, r1, r2 -> 
            [sample_id, r1, r2, group, disease, fibrosis, nas, stage]
        }
    )
    
    // Run STAR alignment
    STAR_ALIGN(
        ch_samples.map { sample_id, group, disease, fibrosis, nas, stage, r1, r2 -> 
            [sample_id, r1, r2, group, disease, fibrosis, nas, stage]
        },
        file(params.genome),
        file(params.gtf),
        file(params.star_index)
    )
    
    // Run featureCounts
    FEATURECOUNTS(STAR_ALIGN.out.bam, file(params.gtf))
    
    // Collect all count files and sample metadata for DESeq2
    ch_counts_files = FEATURECOUNTS.out.counts.collect()
    
    // Extract metadata from ch_samples and collect it
    ch_sample_metadata = ch_samples.map { sample_id, group, disease, fibrosis, nas, stage, r1, r2 -> 
        return [sample_id, group, disease, fibrosis, nas, stage]
    }.collect()
    
    // Run DESeq2 analysis
    DESEQ2(ch_counts_files, ch_sample_metadata)
    
    // Collect QC reports for MultiQC
    reports = Channel.empty()
    reports = reports.mix(FASTQC.out.zip.collect())
    reports = reports.mix(STAR_ALIGN.out.logs.collect())
    reports = reports.mix(FEATURECOUNTS.out.summary.collect())
    
    // Run MultiQC to summarize results
    MULTIQC(reports.collect())
}

// Workflow completion notification
workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    """
} 