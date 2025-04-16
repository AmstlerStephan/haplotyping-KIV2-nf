nextflow.enable.dsl = 2

include {EXTRACT_HAPLOTYPES} from '../processes/extract_haplotypes.nf'
include {MERGE_HAPLOTYPES} from '../processes/merge_haplotypes.nf'
include {MULTIPLE_ALIGNMENT} from '../processes/multiple_alignment.nf'
include {FILTER_BAM} from '../processes/filter_bam.nf'

workflow EXTRACT_HAPLOTYPES_WF {

    requiredParams = [
        'input', 'output'
    ]

    requiredParams.each { param -> 
        if (!params[param]) {
        exit 1, "Parameter ${param} is required."
        }
    }

    // scripts
    extract_haplotypes_py = file( "${projectDir}/bin/extract_haplotypes.py", checkIfExists: true)
    merge_haplotypes_py = file( "${projectDir}/bin/merge_haplotypes.py", checkIfExists: true)
    filter_bam_py = file( "${projectDir}/bin/filter_bam.py", checkIfExists: true)

    // files
    if (params.use_variant_calling_positions) {
        variant_calling_positions = file( "${params.variant_calling_positions}", checkIfExists: true)
    } else {
        variant_calling_positions = file( "${projectDir}/data/variant_calling/NO_FILE.txt", checkIfExists: true)
    }

    // STAGE CHANNELS
    bam_files = 
        Channel.fromPath("${params.input}/barcode*/*/align/consensus/${params.bam_pattern}", type: "file")
        .map { file -> tuple(file.parent.parent.parent.name, file) }

    bam_file_indexes = 
        Channel.fromPath("${params.input}/barcode*/*/align/consensus/${params.bam_pattern}.bai", type: "file")
        .map { file -> tuple(file.parent.parent.parent.name, file) }

    cluster_stats = Channel.fromPath("${params.input}/barcode*/*/stats/raw/${params.cluster_stats_pattern}")
        .map { file -> tuple(file.parent.parent.parent.name, file) }

    bam_stats_tuples = bam_files
        .join(bam_file_indexes)
        .join(cluster_stats)
        .map { barcode, bam, index, stats ->
            tuple(barcode, bam, index, stats)
        }

    FILTER_BAM(bam_stats_tuples, filter_bam_py)
    
    EXTRACT_HAPLOTYPES(FILTER_BAM.out.filtered_bam, variant_calling_positions, extract_haplotypes_py)

    extracted_haplotypes_filtered = EXTRACT_HAPLOTYPES.out.extracted_haplotypes
    filter{ barcode, fasta_file -> fasta_file.countFasta() >= 1 }
    
    MERGE_HAPLOTYPES(extracted_haplotypes_filtered, merge_haplotypes_py)

    //MULTIPLE_ALIGNMENT(MERGE_HAPLOTYPES.out.merged_haplotypes)
    
}
