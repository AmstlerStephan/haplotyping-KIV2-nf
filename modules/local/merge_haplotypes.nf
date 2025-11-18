process MERGE_HAPLOTYPES {
  tag "${sample}-${region}"
  publishDir "${params.output}/${sample}/${region}/haplotyping/", pattern: "*fasta", mode: 'copy'
  publishDir "${params.output}/${sample}/${region}/stats/", pattern: "*tsv", mode: 'copy'

  input:
  tuple val(sample), val(region), path(fastx_file)
  path merge_haplotypes_py

  output:
  tuple val("${sample}"), val("${region}"), path("*merged_haplotypes.fasta"), emit: merged_haplotypes
  path "*fasta"
  path "*.tsv"

  script:
  """
    python ${merge_haplotypes_py} \\
        --fastx_file ${fastx_file} \\
        --variant_cutoff ${params.variant_cutoff} \\
        --output_format ${params.output_format} \\
        --max_edit_distance ${params.max_edit_distance} \\
        -o ./
    """
}
