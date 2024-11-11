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
    bam_file_paths = 
    Channel.fromPath("${params.input}/barcode*/align/consensus/${params.bam_pattern}", type: "file")

    bam_file_index_paths = 
    Channel.fromPath("${params.input}/barcode*/align/consensus/${params.bam_pattern}.bai", type: "file")

    cluster_stats_paths = 
    Channel.fromPath("${params.input}/barcode*/stats/raw/${params.cluster_stats_pattern}", type: "file")

    // wait for all files to be found
    bam_file_paths.await()
    bam_file_index_paths.await()
    cluster_stats_paths.await()

    bam_file_paths
    .map { 
        bam_file_path -> 
            def barcode = (bam_file_path =~ /barcode\d*/)[0]
            tuple ( barcode, bam_file_path)
    }
    .set { bam_files }

    cluster_stats_paths
    .map { 
        cluster_stats_path -> 
            def barcode = (cluster_stats_path =~ /barcode\d*/)[0]
            tuple ( barcode, cluster_stats_path)
    }
    .set { cluster_stats }

    bam_file_index_paths
    .map { 
        bam_file_index_path -> 
            def barcode = (bam_file_index_path =~ /barcode\d*/)[0]
            tuple ( barcode, bam_file_index_path)
    }
    .set { bam_file_indexes }

    bam_files
    .join(bam_file_indexes, remainder: false)
    .join(cluster_stats, remainder: false)
    .set { bam_stats_tuples }

    FILTER_BAM(bam_stats_tuples, filter_bam_py)
    
    EXTRACT_HAPLOTYPES(FILTER_BAM.out.filtered_bam, variant_calling_positions, extract_haplotypes_py)

    extracted_haplotypes_filtered = EXTRACT_HAPLOTYPES.out.extracted_haplotypes
    filter{ barcode, fasta_file -> fasta_file.countFasta() >= 1 }
    
    MERGE_HAPLOTYPES(extracted_haplotypes_filtered, merge_haplotypes_py)

    //MULTIPLE_ALIGNMENT(MERGE_HAPLOTYPES.out.merged_haplotypes)
    
}
