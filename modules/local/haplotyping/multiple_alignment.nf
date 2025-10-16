process MULTIPLE_ALIGNMENT {
  tag "${sample}"
  publishDir "${params.output}/${sample}/mafft/tree", pattern: "*.tree", mode: 'copy'
  publishDir "${params.output}/${sample}/mafft", pattern: "merged_haplotypes_ma.fasta", mode: 'copy'

  input:
  tuple val(sample), path(merged_haplotypes)

  output:
  path "*.tree"
  path "merged_haplotypes_ma.fasta"

  script:
  """
    ginsi \\
        --thread ${task.cpus} \\
        --treeout \\
        ${merged_haplotypes} > merged_haplotypes_ma.fasta
    """
}
