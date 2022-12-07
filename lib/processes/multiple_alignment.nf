process MULTIPLE_ALIGNMENT {
    publishDir "${params.output}/${run}/${barcode}/fasta/", mode: 'copy'
  input:
    tuple val( run ), val( barcode ), path( fasta_clipped )
  output:
    tuple val( "${run}" ), val( "${barcode}" ), path( "multiple_alignment.fasta" ), emit: fasta_aligned
  script:
  """
    mafft \
    --thread ${params.threads} \
    --anysymbol \
    --bl 62 \
    --op 1.53 \
    --ep 0.123 \
    --reorder \
    --retree 2 \
    --treeout \
    --maxiterate 2 \
    ${fasta_clipped} > multiple_alignment.fasta
  """
}