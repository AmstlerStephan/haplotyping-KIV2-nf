nextflow.enable.dsl = 2

requiredParams = [
    'input', 'all_runs', 'output'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}

// scripts
clip_sequences = file( "${projectDir}/bin/clip_sequences.R", checkIfExists: true)
extract_haplotype = file( "${projectDir}/bin/extract_haplotypes.R", checkIfExists: true)

// STAGE CHANNELS
if (params.all_runs) {
    bam_files = Channel.fromPath("${params.input}/run*/ont_pl/barcode*/**${params.bam_pattern}", type: 'file')
    .view()

    sample_sheets = [:]
    Channel.fromPath("${params.input}/run*/lib/*${params.sample_sheet}", type: 'file')
    .view()
    .map { 
        sample_sheet_path ->
            run = ( sample_sheet_path =~ /run\d*_*V*\d*/)[0]
            sample_sheets.put("$run", sample_sheet_path)
    }
    .view()
} else {
    bam_files = Channel.fromPath("${params.input}/ont_pl/**${params.bam_pattern}", type: 'file')
    sample_sheets = [:]
    Channel.fromPath("${params.input}/lib/*${params.sample_sheet}", type: 'file')
    .view()
    .map { 
        sample_sheet_path ->
            run = ( sample_sheet_path =~ /run\d*_*V*\d*/)[0]
            sample_sheets.put("$run", sample_sheet_path)
    }
    .view()
}



bam_files
.map { 
    bam_file_path -> 
        run = (bam_file_path =~ /run\d*_*V*\d*/)[0]
        barcode = (bam_file_path =~ /barcode\d*/)[0]
        tuple ( run, barcode, bam_file_path) 
}
.view()
.set { bam_file_tuples }


include {CLIP_SEQUENCES} from '../processes/clip_sequences.nf'
include {MULTIPLE_ALIGNMENT} from '../processes/multiple_alignment.nf'
include {EXTRACT_HAPLOTYPES} from '../processes/extract_haplotypes.nf'

workflow EXTRACT_HAPLOTYPES_WF {

    CLIP_SEQUENCES(bam_file_tuples, clip_sequences)

    CLIP_SEQUENCES.out.fasta_clipped.view()
    MULTIPLE_ALIGNMENT(CLIP_SEQUENCES.out.fasta_clipped)

    MULTIPLE_ALIGNMENT.out.fasta_aligned.view()

    MULTIPLE_ALIGNMENT.out.fasta_aligned
    .map{ run, barcode, fasta_aligned ->
        sample_sheet = sample_sheets.get("$run")
        tuple (run, barcode, fasta_aligned, sample_sheet )
    }
    .view()
    .set{ fasta_aligned_sample_sheet }

    EXTRACT_HAPLOTYPES(fasta_aligned_sample_sheet, extract_haplotype)
}