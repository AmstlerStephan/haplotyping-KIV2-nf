process EXTRACT_HAPLOTYPES {
  tag "${sample}-${region}"
  publishDir "${params.output}/${sample}/${region}/haplotyping/", mode: 'copy', pattern: "*${params.output_format}"
  publishDir "${params.output}/${sample}/${region}/stats/", mode: 'copy', pattern: "*stats.tsv"

  input:
  tuple val(sample), val(region), path(bam_file), path(bam_file_index)
  path variant_calling_positions
  path extract_haplotypes_py

  output:
  tuple val("${sample}"), val("${region}"), path("haplotypes_filtered.${params.output_format}"), emit: extracted_haplotypes
  path "haplotypes.${params.output_format}"
  path "*stats.tsv"

  script:
  def hardmask = params.hardmask ? "--hardmask" : ""
  def use_variant_calling_positions = params.use_variant_calling_positions ? "--use_variant_calling_positions" : ""
  def variant_calling_positions_arg = params.use_variant_calling_positions ? "--variant_calling_positions ${variant_calling_positions}" : ""
  
  // Get region-specific exclusion ranges
  def exclusion_ranges = params.region_exclusion_ranges?.get(region) ?: ""
  def ranges_arg = exclusion_ranges && !exclusion_ranges.isEmpty() ? "--ranges_to_exclude ${exclusion_ranges}" : ""
  """
    python ${extract_haplotypes_py} \\
        --bam_file ${bam_file} \\
        ${ranges_arg} \\
        --min_qscore ${params.min_qscore} \\
        ${hardmask} \\
        ${use_variant_calling_positions} \\
        ${variant_calling_positions_arg} \\
        --output_format ${params.output_format} \\
        -o ./
    """
}
