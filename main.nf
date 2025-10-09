#!/usr/bin/env nextflow

// ENABLE DSL2
nextflow.enable.dsl=2

include {EXTRACT_HAPLOTYPES_WF} from './lib/workflows/extract_haplotypes.nf'


// Main Workflow 
workflow {
    // Print help and exit
    if (params.help) {
        println """\
        ===============================================
        U M I   H A P L O T Y P I N G   P I P E L I N E
        ===============================================
        ~ version ${workflow.manifest.version}

        Usage:
        nextflow run AmstlerStephan/haplotyping-KIV2-nf [OPTIONS]...

        Options:
        --input [path/to/input/dir]                       [REQUIRED] Provide the directory containing input BAM files.
        --ont_pl_dir [path/to/ont_pl_dir]                 [REQUIRED] Provide the directory associated with ONT pipeline.
        --output [path/to/output/dir]                     [REQUIRED] Specify the directory where the pipeline will write its output.
        --variant_calling_positions [path/to/positions]   Specify the path to the file specifying variant calling positions.
        --region [region name]                            Specify the region name to be processed.
        --bam_pattern [pattern]                           Pattern to match BAM files within the input directory. Default: "masked_consensus.bam"
        --cluster_stats_pattern [pattern]                 Pattern to match cluster statistics files. Default: "split_cluster_stats.tsv"
        --min_reads_per_cluster [number]                  Minimum number of reads per cluster. Default: 10
        --max_reads_per_cluster [number]                  Maximum number of reads per cluster. Default: 200
        --max_edit_distance [number]                      Maximum edit distance allowed. Default: 2
        --use_variant_calling_positions                   Use variant calling positions. Default: false
        --ranges_to_exclude [ranges]                      Comma-separated list of ranges to exclude. Default: "2472,2506"
        --min_qscore [number]                             Minimum quality score required. Default: 45
        --output_format [format]                          Output format for haplotype results. Default: "fasta"

        --threads [number]                                Number of threads used during execution. Default: (availableProcessors() - 1)
        """
        exit 0
    }

    // Print version and exit
    if (params.version) {
        println """\
        ===============================================
        U M I   H A P L O T Y P I N G   P I P E L I N E
        ===============================================
        ~ version ${workflow.manifest.version}
        """
        exit 0
    }

    // Print standard logging info
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

    EXTRACT_HAPLOTYPES_WF()


    // Workflow tracing - what to display when the pipeline finishes
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

        // Run a small clean-up script to remove "work" directory after successful completion 
        if (!params.debug && workflow.success) {
            ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
        }
    }
}

