process MERGE_AND_RECONSTRUCT_HAPLOTYPES {
  tag "${sample}-${region}"
  publishDir "${params.output}/${sample}/${region}/haplotyping/", pattern: "*.fasta", mode: 'copy'
  publishDir "${params.output}/${sample}/${region}/stats/", pattern: "*.tsv", mode: 'copy'

  input:
  tuple val(sample), val(region), path(fastx_file), path(positions), path(reference)
  path merge_and_reconstruct_py

  output:
  tuple val("${sample}"), val("${region}"), path("merged_haplotypes.fasta"), emit: merged_haplotypes
  tuple val("${sample}"), val("${region}"), path("reconstructed_haplotypes.fasta"), emit: reconstructed_haplotypes
  path "*.fasta"
  path "*.tsv"

  script:
  """
    python ${merge_and_reconstruct_py} \\
        --fastx_file ${fastx_file} \\
        --positions ${positions} \\
        --reference ${reference} \\
        --reference_start ${params.reference_start} \\
        --variant_cutoff ${params.variant_cutoff} \\
        --max_edit_distance ${params.max_edit_distance} \\
        --output_format ${params.output_format} \\
        -o ./
    """
}
