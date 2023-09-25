process EXTRACT_HAPLOTYPES {
    publishDir "${params.output}/${run}/${sample}/haplotype/", pattern: "*fasta", mode: 'copy'
    publishDir "${params.output}/${run}/${sample}/haplotype/stats", pattern: "*tsv", mode: 'copy'
  input:
    tuple val( sample ), path( fastx_file )
    path merge_haplotypes_py
  output:
    tuple val( "${run}" ), val( "${sample}" ), path( "*haplotype.fasta" ), emit: haplotype
    path "*haplotype.tsv"
  script:
  """
    Rscript ${merge_haplotypes_py} \
    --fastx_file ${fastx_file}
      """
}