process FILTER_BAM {
    publishDir "./", mode: 'copy'


  input:
    tuple val( sample ), path( bam_file ), path( bam_file_index ), path( cluster_stats )
    path filter_bam_py


  output:
    tuple val( "${sample}" ), path( "*.bam" ), path( "*.bam.bai" ), emit: filtered_bam
    path "*test.txt"

  script:
  """
  echo $sample and $bam_file and $filter_bam_py > ${sample}_test.txt
    cp ${bam_file} ./test.bam
    cp ${bam_file_index} ./test.bam.bai
  """
}