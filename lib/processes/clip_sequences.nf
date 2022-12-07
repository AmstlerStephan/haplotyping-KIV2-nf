process CLIP_SEQUENCES {
    publishDir "${params.output}/${run}/${barcode}/fasta/", mode: 'copy'
  input:
    tuple val( run ), val( barcode ), path( bam_file )
    path sample_sheet
    path clip_sequences_R
  output:
    tuple val( "${run}" ), val( "${barcode}" ), path( "*clipped.fasta" ), emit: clipped_fasta
    path "*sequences.fasta"
  script:
  """
    Rscript ${clip_sequences_R} --bam_file ${bam_file}
  """
}