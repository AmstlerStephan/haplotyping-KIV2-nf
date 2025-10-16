import nextflow.Nextflow

class WorkflowMain {

    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis, please cite:\n\n" +
            "Nanopore sequencing with unique molecular identifiers enables accurate mutation analysis " +
            "and haplotyping in the complex lipoprotein(a) KIV-2 VNTR. Genome Med 16, 117 (2024). " +
            "https://doi.org/10.1186/s13073-024-01391-8";
    }

    public static void validate(params) {
        def requiredParams = [
            'input', 'output'
        ]

        requiredParams.each { param ->
            if (params[param] == null) {
                Nextflow.error("Parameter ${param} is required to run the pipeline.")
            }
        }

        // Validate input directory exists
        if (params.input && !file(params.input).exists()) {
            Nextflow.error("Input directory does not exist: ${params.input}")
        }

        // Validate region parameter if provided
        if (params.region && params.region.trim().isEmpty()) {
            Nextflow.error("Parameter 'region' cannot be empty if specified.")
        }

        // Validate numeric parameters
        if (params.min_reads_per_cluster < 1) {
            Nextflow.error("Parameter 'min_reads_per_cluster' must be >= 1")
        }

        if (params.max_reads_per_cluster < params.min_reads_per_cluster) {
            Nextflow.error("Parameter 'max_reads_per_cluster' must be >= min_reads_per_cluster")
        }

        if (params.min_qscore < 0 || params.min_qscore > 60) {
            Nextflow.error("Parameter 'min_qscore' must be between 0 and 60")
        }

        // Validate output format
        def validFormats = ['fasta', 'fastq']
        if (!(params.output_format in validFormats)) {
            Nextflow.error("Parameter 'output_format' must be one of: ${validFormats.join(', ')}")
        }
    }

    public static String version(workflow) {
        return "" +
            "          ===================                          \n" +
            "            ${workflow.manifest.name.toUpperCase()}    \n" +
            "          ===================                          \n" +
            "            version ${workflow.manifest.version}       \n";
    }

    public static String onComplete(workflow, baseDir, params) {
        // run a small clean-up script to remove "work" directory after successful completion
        if (workflow.success && !params.verbose) {
            ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
        }
    }
}