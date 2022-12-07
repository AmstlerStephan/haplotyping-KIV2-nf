process EXTRACT_HAPLOTYPES {
    publishDir "${params.output}/${run}/${barcode}/haplotype/", pattern:"fasta", mode: 'copy'
    publishDir "${params.output}/${run}/${barcode}/haplotype/stats", pattern:"tsv", mode: 'copy'
  input:
    tuple val( run ), val( barcode ), path( fasta_aligned ), path ( sample_sheet )
    path extract_haplotypes_R
  output:
    tuple val( "${run}" ), val( "${barcode}" ), path( "*haplotype.fasta" ), emit: haplotype
    path "*haplotype.tsv"
  script:
  """
    Rscript ${extract_haplotypes_R} \
    --sample_sheet ${sample_sheet} \
    --run ${run} \
    --barcode ${barcode} \
    --umi_cutoff_R9 ${params.umi_cutoff_R9} \
    --umi_cutoff_V14 ${params.umi_cutoff_V14} \
    --aligned_fasta ${fasta_aligned}
      """
}