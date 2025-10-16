nextflow.enable.dsl = 2

include { EXTRACT_HAPLOTYPES } from '../modules/local/haplotyping/extract_haplotypes.nf'
include { MERGE_HAPLOTYPES   } from '../modules/local/haplotyping/merge_haplotypes.nf'
include { MULTIPLE_ALIGNMENT } from '../modules/local/haplotyping/multiple_alignment.nf'
include { FILTER_BAM         } from '../modules/local/haplotyping/filter_bam.nf'

workflow EXTRACT_HAPLOTYPES_WF {

  // Validate required parameters
  WorkflowMain.validate(params)

  // scripts
  extract_haplotypes_py = file("${projectDir}/bin/extract_haplotypes.py", checkIfExists: true)
  merge_haplotypes_py = file("${projectDir}/bin/merge_haplotypes.py", checkIfExists: true)
  filter_bam_py = file("${projectDir}/bin/filter_bam.py", checkIfExists: true)

  // files
  if (params.use_variant_calling_positions) {
    variant_calling_positions = file("${params.variant_calling_positions}", checkIfExists: true)
  }
  else {
    variant_calling_positions = file("${projectDir}/data/variant_calling/NO_FILE.txt", checkIfExists: true)
  }

  // Input channels
  bam_files = Channel.fromPath("${params.input}/barcode*/${params.region}/align/consensus/${params.bam_pattern}", type: "file")
    .map { file ->
      def barcode = file.parent.parent.parent.parent.name
      def region = file.parent.parent.parent.name
      tuple(barcode, region, file)
    }

  bam_file_indexes = Channel.fromPath("${params.input}/barcode*/${params.region}/align/consensus/${params.bam_pattern}.bai", type: "file")
    .map { file ->
      def barcode = file.parent.parent.parent.parent.name
      def region = file.parent.parent.parent.name
      tuple(barcode, region, file)
    }

  cluster_stats = Channel.fromPath("${params.input}/barcode*/${params.region}/stats/raw/${params.cluster_stats_pattern}")
    .map { file ->
      def barcode = file.parent.parent.parent.parent.name
      def region = file.parent.parent.parent.name
      tuple(barcode, region, file)
    }

  bam_stats_tuples = bam_files
    .join(bam_file_indexes)
    .join(cluster_stats)
    .map { barcode, region, bam, index, stats ->
      tuple(barcode, region, bam, index, stats)
    }

  // Process workflow
  FILTER_BAM(bam_stats_tuples, filter_bam_py)

  EXTRACT_HAPLOTYPES(FILTER_BAM.out.filtered_bam, variant_calling_positions, extract_haplotypes_py)

  extracted_haplotypes_filtered = EXTRACT_HAPLOTYPES.out.extracted_haplotypes.filter { _barcode, _region, fasta_file -> fasta_file.countFasta() >= 1 }

  MERGE_HAPLOTYPES(extracted_haplotypes_filtered, merge_haplotypes_py)
}
