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
clip_sequences = file( "${projectDir}/bin/parse_nanostat.R", checkIfExists: true)
merge_parsed_run = file( "${projectDir}/bin/merge_parsed_run.R", checkIfExists: true)
parse_metrics = file( "${projectDir}/bin/parse_run_metrics.R", checkIfExists: true)

// STAGE CHANNELS
if (params.all_runs) {
    barcodes = Channel.fromPath("${params.input}/run*/fastq_pass/barcode*", type: 'dir') 
} else {
    barcodes = Channel.fromPath("${params.input}/fastq_pass/barcode*", type = 'dir')
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

run_metrics = Channel.fromPath("${params.input}/**/*.md", type: 'file')
.map { 
    run_metrics_path ->
        run = ( run_metrics_path =~ /run\d*_*V*\d*/)[0]
        tuple( run, run_metrics_path )
}


workflow EXTRACT_HAPLOTYPES {


}