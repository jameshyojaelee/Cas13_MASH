process STAR_ALIGN {
    tag "$sample"
    label 'high_memory'
    
    container 'quay.io/biocontainers/star:2.7.10b--h9ee0642_0'
    
    input:
    tuple val(sample), path(read1), path(read2), val(group), val(disease), val(fibrosis), val(nas), val(stage)
    path genome
    path gtf
    path star_index
    
    output:
    tuple val(sample), path("${sample}Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample}Log.final.out", emit: logs
    
    script:
    def index_cmd = star_index ? "" : """
        mkdir star_index
        STAR --runMode genomeGenerate \
            --genomeDir star_index \
            --genomeFastaFiles ${genome} \
            --sjdbGTFfile ${gtf} \
            --runThreadN ${task.cpus}
        """
    def index_path = star_index ?: "star_index"
    """
    ${index_cmd}
    
    STAR --genomeDir ${index_path} \
        --readFilesIn ${read1} ${read2} \
        --runThreadN ${task.cpus} \
        --outFileNamePrefix ${sample} \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --quantMode GeneCounts
    """
} 