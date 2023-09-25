process MERGE_HAPLOTYPES {
    publishDir "${params.output}/${sample}/haplotyping/", mode: 'copy'
  input:
    tuple val( sample ), path( fastx_file )
    path merge_haplotypes_py
  output:
    tuple val( "${sample}" ), path( "*merged_haplotypes.${params.output_format}" ), emit: merged_haplotypes
    path "*${params.output_format}"
    
  script:
  """
    python ${merge_haplotypes_py} \
      --fastx_file ${fastx_file} \
      --variant_cutoff ${params.variant_cutoff} \
      --output_format ${params.output_format} \
      -o ./
    """
}