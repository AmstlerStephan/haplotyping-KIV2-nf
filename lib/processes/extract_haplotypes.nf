process EXTRACT_HAPLOTYPES {
    publishDir "${params.output}/${run}/${barcode}/stats/", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output}/${run}/${barcode}/haplotyping/", mode: 'copy', pattern: "*${params.output_format}"

  input:
    tuple val( run ), val( barcode ), path( bam_file )
    path extract_haplotypes_py
  output:
    tuple val( "${run}" ), val( "${barcode}" ), path( "haplotypes_filtered.${params.output_format}" ), emit: extracted_haplotypes
    path "*tsv"
    path "${params.output_format}"
  script:
    def filter_haplotypes = params.filter_haplotypes ? "--filter_haplotypes" : ""
    def hardmask = params.hardmask ? "--hardmask" : ""
  """
    python ${extract_haplotypes_py} \
      --bam_file ${bam_file} \
      --ranges_to_exclude ${params.ranges_to_exclude} \
      --min_qscore ${params.min_qscore} \
      $filter_haplotypes \
      $hardmask \
      --output_format ${params.output_format} \
      -o ./
  """
}