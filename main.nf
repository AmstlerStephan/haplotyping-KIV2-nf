#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { validateParameters; paramsHelp; paramsSummaryLog  } from 'plugin/nf-validation'
include { EXTRACT_HAPLOTYPES_WF                             } from './workflows/extract_haplotypes.nf'

workflow {

     if (params.help) {
          def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
          def String command = "nextflow run AmstlerStephan/${workflow.manifest.name} -r v0.2.0 -profile test,docker"
          log.info paramsHelp(command) + citation
          exit 0
     }

     if (params.version) {
          log.info WorkflowMain.version(workflow)
          exit 0
     }

     // Validate input parameters (temporarily disabled for map parameter issues)
     // if (params.validate_params & !params.version) {
     //      validateParameters()
     // }
     // Print summary of supplied parameters
     log.info paramsSummaryLog(workflow)

     // Run the workflow
     EXTRACT_HAPLOTYPES_WF()
}

workflow.onError {
     log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
     log.info ""
     log.info "Pipeline completed at: $workflow.complete"
     log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
     WorkflowMain.onComplete(workflow, baseDir, params)
}

