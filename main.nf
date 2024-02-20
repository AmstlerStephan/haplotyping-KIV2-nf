#!/usr/bin/env nextflow

// ENABLE DSL2
nextflow.enable.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         ===============================================
          E C S E Q - t e m p l a t e   P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run ecseq/dnaseq [OPTIONS]...

         Options: GENERAL
              --input [path/to/input/dir]     [REQUIRED] Provide the directory containing fastq file(s) in "*{1,2}.fastq.gz" format

              --reference [path/to/ref.fa]    [REQUIRED] Provide the path to the reference genome in fasta format

              --output [STR]                  A string that can be given to name the output directory. [default: "."]


         Options: MODIFIERS
              --SE                            Indicate to the pipeline whether fastq files are SE reads in "*.fastq.gz" format. [default: off]

              --FastQC                        Generate FastQC report of trimmed reads. [default: off]

              --bamQC                         Generate bamQC report of alignments. [default: off]

              --keepReads                     Keep trimmed fastq reads. [default: off]


         Options: TRIMMING
              --forward                       Forward adapter sequence. [default: "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"]

              --reverse                       Reverse adapter sequence. [default: "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"]

              --minQual                       Minimum base quality threshold. [default: 20]

              --minLeng                       Minimum read length threshold. [default: 25]

              --minOver                       Minimum overlap threshold. [default: 3]


         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    


         Example: 
              nextflow run ecseq/dnaseq \
              --input /path/to/input/dir --reference /path/to/genome.fa \
              --FastQC --bamQC --keepReads

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\

         ===============================================
          Nextflow UMI haplotyping Pipeline
         ===============================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         ==============================================="
log.info "          U M I - H A P L O T Y P I N G"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ===============================================" }
else {
log.info "         ===============================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir    : ${params.input}"
log.info "         output dir   : ${params.output}"
log.info "         ==============================================="
log.info "         RUN NAME: ${workflow.runName}"
log.info ""

include {EXTRACT_HAPLOTYPES_WF} from './lib/workflows/extract_haplotypes.nf'

// MAIN WORKFLOW 
workflow {
        EXTRACT_HAPLOTYPES_WF()
}


//////////////////
// END PIPELINE //
//////////////////

// WORKFLOW TRACING # what to display when the pipeline finishes
// eg. with errors
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// eg. in general
workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${!params.debug && workflow.success ? "(cleared)" : "" }"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    // run a small clean-up script to remove "work" directory after successful completion 
    if (!params.debug && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}
