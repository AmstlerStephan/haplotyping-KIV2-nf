nextflow.enable.dsl = 2

include { EXTRACT_HAPLOTYPES } from '../modules/local/extract_haplotypes.nf'
include { MERGE_HAPLOTYPES   } from '../modules/local/merge_haplotypes.nf'
include { MULTIPLE_ALIGNMENT } from '../modules/local/multiple_alignment.nf'
include { FILTER_BAM         } from '../modules/local/filter_bam.nf'

workflow EXTRACT_HAPLOTYPES_WF {

  // Validate required parameters (temporarily disabled)
  // WorkflowMain.validate(params)

  // scripts
  extract_haplotypes_py = file("${projectDir}/bin/extract_haplotypes.py", checkIfExists: true)
  merge_haplotypes_py = file("${projectDir}/bin/merge_haplotypes.py", checkIfExists: true)
  filter_bam_py = file("${projectDir}/bin/filter_bam.py", checkIfExists: true)

  // Set up variant calling positions file
  def no_file = file("${projectDir}/data/variant_calling/NO_FILE/NO_FILE.txt", checkIfExists: true)

  // Input channels - process all regions
  bam_files = Channel.fromPath("${params.input}/barcode*/*/align/consensus/${params.bam_pattern}", type: "file")
    .map { file ->
      def barcode = file.parent.parent.parent.parent.name
      def region = file.parent.parent.parent.name
      tuple(barcode, region, file)
    }

  bam_file_indexes = Channel.fromPath("${params.input}/barcode*/*/align/consensus/${params.bam_pattern}.bai", type: "file")
    .map { file ->
      def barcode = file.parent.parent.parent.parent.name
      def region = file.parent.parent.parent.name
      tuple(barcode, region, file)
    }

  cluster_stats = Channel.fromPath("${params.input}/barcode*/*/stats/raw/${params.cluster_stats_pattern}")
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

  // Add region-specific variant positions file to all filtered BAM samples
  bam_with_variant_positions = FILTER_BAM.out.filtered_bam.map { sample, region, bam, bai ->
    def region_variant_file = params.region_variant_calling_positions?.get(region) ?: ""
    def variant_file = region_variant_file && !region_variant_file.isEmpty() && params.use_variant_calling_positions
      ? file(region_variant_file, checkIfExists: true)
      : no_file
    tuple(sample, region, bam, bai, variant_file)
  }

  EXTRACT_HAPLOTYPES(bam_with_variant_positions, extract_haplotypes_py)

  extracted_haplotypes_filtered = EXTRACT_HAPLOTYPES.out.extracted_haplotypes.filter { _barcode, _region, fasta_file -> fasta_file.countFasta() >= 1 }

  MERGE_HAPLOTYPES(extracted_haplotypes_filtered, merge_haplotypes_py)
}
