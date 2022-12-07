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
    barcodes = Channel.fromPath("${params.input}/run*/ont_pl/barcode*/**${params.bam_pattern}", type: 'file') 
} else {
    barcodes = Channel.fromPath("${params.input}/ont_pl/barcode*/**${params.bam_pattern}", type = 'file')
}
barcodes_tuple = barcodes
.map { 
    barcode_path -> 
        run = (barcode_path =~ /run\d*_*V*\d*/)[0]
        barcode = barcode_path.baseName
        tuple( run, barcode, barcode_path ) 
}

sample_sheets = [:]
Channel.fromPath("${params.input}/**/${params.sample_sheet}", type: 'file')
.map { 
    sample_sheet_path ->
        run = ( sample_sheet_path =~ /run\d*_*V*\d*/)[0]
        sample_sheets.put("$run", sample_sheet_path)
}

include {CLIP_SEQUENCES} from '../processes/clip_sequences.nf'
include {EXTRACT_HAPLOTYPES} from '../processes/extract_haplotypes.nf'
include {MULTIPLE_ALIGNMENT} from '../processes/multiple_alignment.nf'

workflow EXTRACT_HAPLOTYPES_WF {


}