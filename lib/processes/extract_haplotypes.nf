process EXTRACT_HAPLOTYPES {
    publishDir "${params.output}/${run}/${sample}/haplotyping/", mode: 'copy', pattern: "*${params.output_format}"
    publishDir "${params.output}/${run}/${sample}/stats/", mode: 'copy', pattern: "haplotype_stats.tsv"

  input:
    tuple val( sample ), path( bam_file )
    path extract_haplotypes_py
  output:
    tuple val( "${sample}" ), path( "haplotypes_filtered.${params.output_format}" ), emit: extracted_haplotypes
    path "*${params.output_format}"
    path "haplotype_stats.tsv"
  script:
    def hardmask = params.hardmask ? "--hardmask" : ""
  """
    python ${extract_haplotypes_py} \
      --bam_file ${bam_file} \
      --ranges_to_exclude ${params.ranges_to_exclude} \
      --min_qscore ${params.min_qscore} \
      $hardmask \
      --output_format ${params.output_format} \
      -o ./
  """
}