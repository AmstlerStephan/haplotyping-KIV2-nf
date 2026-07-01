nextflow.enable.dsl = 2

include { EXTRACT_HAPLOTYPES               } from '../modules/local/extract_haplotypes.nf'
include { MERGE_AND_RECONSTRUCT_HAPLOTYPES } from '../modules/local/merge_and_reconstruct_haplotypes.nf'
include { MULTIPLE_ALIGNMENT               } from '../modules/local/multiple_alignment.nf'
include { FILTER_BAM                       } from '../modules/local/filter_bam.nf'

workflow EXTRACT_HAPLOTYPES_WF {

  // Validate required parameters (temporarily disabled)
  // WorkflowMain.validate(params)

  // scripts
  extract_haplotypes_py = file("${projectDir}/bin/extract_haplotypes.py", checkIfExists: true)
  merge_and_reconstruct_haplotypes_py = file("${projectDir}/bin/merge_and_reconstruct_haplotypes.py", checkIfExists: true)
  filter_bam_py = file("${projectDir}/bin/filter_bam.py", checkIfExists: true)


  // Input channels - process all regions
  bam_files = Channel
    .fromPath("${params.input}/barcode*/*/align/consensus/${params.bam_pattern}", type: "file")
    .map { file ->
      def barcode = file.parent.parent.parent.parent.name
      def region = file.parent.parent.parent.name
      tuple(barcode, region, file)
    }

  bam_file_indexes = Channel
    .fromPath("${params.input}/barcode*/*/align/consensus/${params.bam_pattern}.bai", type: "file")
    .map { file ->
      def barcode = file.parent.parent.parent.parent.name
      def region = file.parent.parent.parent.name
      tuple(barcode, region, file)
    }

  cluster_stats = Channel
    .fromPath("${params.input}/barcode*/*/stats/raw/${params.cluster_stats_pattern}")
    .map { file ->
      def barcode = file.parent.parent.parent.parent.name
      def region = file.parent.parent.parent.name
      tuple(barcode, region, file)
    }

  bam_stats_tuples = bam_files
    .join(bam_file_indexes, by: [0, 1])
    .join(cluster_stats, by: [0, 1])
    .map { barcode, region, bam, index, stats ->
      tuple(barcode, region, bam, index, stats)
    }

  // Process workflow
  FILTER_BAM(bam_stats_tuples, filter_bam_py)

  // Add optional region-specific variant positions file to filtered BAM samples.
  // When use_variant_calling_positions=false, pass an empty string and let
  // EXTRACT_HAPLOTYPES discover polymorphic sites directly from the BAM.
  bam_with_variant_positions = FILTER_BAM.out.filtered_bam.map { sample, region, bam, bai ->
    def region_variant_file = params.region_variant_calling_positions?.get(region) ?: ""
    def variant_file = params.use_variant_calling_positions ? region_variant_file : ""

    if (params.use_variant_calling_positions && (!variant_file || variant_file.isEmpty())) {
      throw new IllegalArgumentException("Missing region_variant_calling_positions entry for region '${region}' while use_variant_calling_positions=true")
    }

    tuple(sample, region, bam, bai, variant_file)
  }

  EXTRACT_HAPLOTYPES(bam_with_variant_positions, extract_haplotypes_py)

  extracted_haplotypes_filtered = EXTRACT_HAPLOTYPES.out.extracted_haplotypes.filter { _barcode, _region, fasta_file -> fasta_file.countFasta() >= 1 }

  // Prepare inputs for fused merge+reconstruct process:
  // Join extracted haplotypes with their polymorphic positions and region reference.
  haplotypes_with_positions_and_ref = extracted_haplotypes_filtered
    .join(EXTRACT_HAPLOTYPES.out.positions, by: [0, 1])
    .map { sample, region, haplotypes, positions ->
      def region_ref = params.region_references?.get(region) ?: ""
      if (!region_ref || region_ref.isEmpty()) {
        throw new IllegalArgumentException("Missing region_references entry for region '${region}' required for full sequence reconstruction")
      }
      def ref_file = file(region_ref, checkIfExists: true)
      tuple(sample, region, haplotypes, positions, ref_file)
    }

  // Fused process: merge haplotypes AND reconstruct full-length sequences in one pass
  MERGE_AND_RECONSTRUCT_HAPLOTYPES(haplotypes_with_positions_and_ref, merge_and_reconstruct_haplotypes_py)

  // Feed reconstructed sequences into MULTIPLE_ALIGNMENT

  if (params.perform_multiple_alignment) {
    MULTIPLE_ALIGNMENT(MERGE_AND_RECONSTRUCT_HAPLOTYPES.out.reconstructed_haplotypes)
  }
}
