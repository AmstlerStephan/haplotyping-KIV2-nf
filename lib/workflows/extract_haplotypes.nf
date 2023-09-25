nextflow.enable.dsl = 2

requiredParams = [
    'input', 'output'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}

// scripts
extract_haplotypes_py = file( "${projectDir}/bin/extract_haplotypes.py", checkIfExists: true)
merge_haplotypes_py = file( "${projectDir}/bin/merge_haplotypes.py", checkIfExists: true)

// STAGE CHANNELS
bam_file_paths = Channel.fromPath("${params.input}/barcode*/align/**${params.bam_pattern}", type: "file")
bam_file_index_paths = Channel.fromPath("${params.input}/barcode*/align/**${params.bam_pattern}.bai", type: "file")

bam_files = bam_file_paths
.map { 
    bam_file_path -> 
        barcode = (bam_file_path =~ /barcode\d*/)[0]
        tuple ( barcode, bam_file_path)
}

bam_file_indexes = bam_file_index_paths
.map { 
    bam_file_index_path -> 
        barcode = (bam_file_index_path =~ /barcode\d*/)[0]
        tuple ( barcode, bam_file_index_path)
}

bam_file_tuples = 
bam_files
.join(bam_file_indexes)

include {EXTRACT_HAPLOTYPES} from '../processes/extract_haplotypes.nf'
include {MERGE_HAPLOTYPES} from '../processes/merge_haplotypes.nf'
include {MULTIPLE_ALIGNMENT} from '../processes/multiple_alignment.nf'

workflow EXTRACT_HAPLOTYPES_WF {

    EXTRACT_HAPLOTYPES(bam_file_tuples, extract_haplotypes_py)

    extracted_haplotypes_filtered = EXTRACT_HAPLOTYPES.out.extracted_haplotypes
    filter{ barcode, fasta_file -> fasta_file.countFasta() >= 1 }
    
    MERGE_HAPLOTYPES(extracted_haplotypes_filtered, merge_haplotypes_py)

    MULTIPLE_ALIGNMENT(MERGE_HAPLOTYPES.out.merged_haplotypes)
}
