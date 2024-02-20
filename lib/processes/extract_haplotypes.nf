process EXTRACT_HAPLOTYPES {
    publishDir "${params.output}/${sample}/haplotyping/", mode: 'copy', pattern: "*${params.output_format}"
    publishDir "${params.output}/${sample}/stats/", mode: 'copy', pattern: "haplotype_stats.tsv"

  input:
    tuple val( sample ), path( bam_file ), path( bam_file_index )
    path variant_calling_positions
    path extract_haplotypes_py
  output:
    tuple val( "${sample}" ), path( "haplotypes_filtered.${params.output_format}" ), emit: extracted_haplotypes
    path "haplotypes.${params.output_format}"
    path "haplotype_stats.tsv"
  script:
    def hardmask = params.hardmask ? "--hardmask" : ""
    def use_variant_calling_positions = params.use_variant_calling_positions ? "--use_variant_calling_positions" : ""
    def variant_calling_positions = params.use_variant_calling_positions ? "--variant_calling_positions ${variant_calling_positions}" : ""
  """
    python ${extract_haplotypes_py} \
      --bam_file ${bam_file} \
      --ranges_to_exclude ${params.ranges_to_exclude} \
      --min_qscore ${params.min_qscore} \
      $hardmask \
      $use_variant_calling_positions \
      $variant_calling_positions \
      --output_format ${params.output_format} \
      -o ./
  """
}