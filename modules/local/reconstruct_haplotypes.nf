process RECONSTRUCT_HAPLOTYPES {
  tag "${sample}-${region}"
  publishDir "${params.output}/${sample}/${region}/haplotyping/", pattern: "*.fasta", mode: 'copy'

  input:
  tuple val(sample), val(region), path(merged_haplotypes), path(positions), path(reference)
  path reconstruct_haplotypes_py

  output:
  tuple val("${sample}"), val("${region}"), path("reconstructed_haplotypes.fasta"), emit: reconstructed_haplotypes

  script:
  """
    python ${reconstruct_haplotypes_py} \\
        --merged_haplotypes ${merged_haplotypes} \\
        --positions ${positions} \\
        --reference ${reference} \\
        --reference_start ${params.reference_start} \\
        -o ./
    """
}
