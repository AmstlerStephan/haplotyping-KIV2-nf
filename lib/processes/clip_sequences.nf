process CLIP_SEQUENCES {
    publishDir "${params.output}/${run}/${barcode}/fasta/", mode: 'copy'
  input:
    tuple val( run ), val( barcode ), path( bam_file )
    path clip_sequences_R
  output:
    tuple val( "${run}" ), val( "${barcode}" ), path( "*clipped.fasta" ), emit: fasta_clipped
    path "*sequences.fasta"
  script:
  """
    Rscript ${clip_sequences_R} --bam_file ${bam_file}
  """
}