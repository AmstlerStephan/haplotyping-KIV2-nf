mafft_fasta="merged_haplotypes_ma.fasta"
process MULTIPLE_ALIGNMENT {
    publishDir "${params.output}/${sample}/haplotyping/mafft/tree", pattern: "*.tree", mode: 'copy'
    publishDir "${params.output}/${sample}/haplotyping/mafft", mode: 'copy'
  input:
    tuple val( sample ), path( merged_haplotypes )
  output:
    tuple val( "${sample}" ), path( "*.tree" )
    path "${mafft_fasta}"
  script:
  """
    ginsi \
    --thread ${params.threads} \
    --treeout \
    ${merged_haplotypes} > ${mafft_fasta}
  """
}