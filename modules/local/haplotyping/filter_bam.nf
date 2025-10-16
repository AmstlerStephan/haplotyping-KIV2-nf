process FILTER_BAM {
  tag "${sample}"
  publishDir "${params.output}/${sample}/filtered_bam/", mode: 'copy', enabled: "${params.verbose}"

  input:
  tuple val(sample), path(bam_file), path(bam_file_index), path(cluster_stats)
  path filter_bam_py

  output:
  tuple val("${sample}"), path("*.bam"), path("*.bam.bai"), emit: filtered_bam

  script:
  """
    python ${filter_bam_py} \\
        --bam_file ${bam_file} \\
        --cluster_stats ${cluster_stats} \\
        --min_reads_per_cluster ${params.min_reads_per_cluster} \\
        --max_reads_per_cluster ${params.max_reads_per_cluster} \\
        -o ./

    samtools sort filtered_bam.bam -o filtered_bam.bam && samtools index filtered_bam.bam
    """
}
