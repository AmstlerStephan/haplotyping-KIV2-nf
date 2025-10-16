[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-brightgreen.svg)](https://www.nextflow.io/)

haplotyping-KIV2-nf Pipeline
======================

**haplotyping-KIV2-nf** is a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

## Overview
`haplotyping-KIV2-nf` is designed for haplotyping using a series of Python scripts. The workflow involves multiple stages such as filtering BAM files, extracting haplotypes, and merging haplotype information. The pipeline extracts the polymorphic positions of aligned reads provided in bam-file format. It removes the soft- and hard-clips of the sequences and extracts the polymorphic positions (haplotype) per sample either per-sample or using user defined positions.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/)

To run the haplotyping workflow, execute the following command:

```bash
nextflow run AmstlerStephan/haplotyping-KIV2-nf -r v0.1.1 -c <custom.config> -profile <docker/conda> 
```

Replace `<path_to_input_directory>` and `<path_to_output_directory>` with the actual paths.

## Scripts and Files
- The workflow utilizes several Python scripts:
  - `extract_haplotypes.py`: Used for haplotype extraction.
  - `merge_haplotypes.py`: Used for merging haplotype information.
  - `filter_bam.py`: Used for filtering BAM files.
- Additionally, a file named `variant_calling_positions` is used for haplotype extraction. Its location is determined by the `use_variant_calling_positions` parameter.

## Input Channels
- Input data is read from specific directories (`input/barcode*/align/consensus/`). BAM files, BAM file indexes, and cluster statistics are organized into tuples based on barcode information.

## Workflow Stages
1. **Filter BAM Files**
   - Filters BAM files using the `filter_bam.py` script based on specified criteria.

2. **Extract Haplotypes**
   - Uses the `extract_haplotypes.py` script to extract haplotypes from the filtered BAM files. Haplotypes are filtered based on the provided `variant_calling_positions` file.

3. **Merge Haplotypes**
   - Merges the extracted haplotypes using the `merge_haplotypes.py` script.

## Configuration
- The workflow may have configurable parameters in the `nextflow.config` file. Check for customization options there.

### Basic Parameters

- **help**: A boolean flag indicating whether to display help information. Default is `false`.
  
- **version**: A boolean flag indicating whether to display the workflow version. Default is `false`.
  
- **debug**: A boolean flag enabling or disabling debug mode. When set to `true`, additional debugging information may be provided during workflow execution. Default is `false`.

### Input/Output Parameters

- **input**: The directory containing input data for the workflow. This parameter is required for the workflow to locate and process input files.

- **ont_pl_dir**: The directory associated with the consensus reads obtained from the https://github.com/genepi/umi-pipeline-nf analysis workflow. Default is `null`.

- **output**: The directory where the workflow will write its output. This parameter is required for storing the results of the haplotyping workflow.

- **variant_calling_positions**: A file specifying variant calling positions. If provided, the workflow uses this file during haplotype extraction.

- **bam_pattern**: The pattern used to match BAM files within the input directory. Default is `"masked_consensus.bam"`.

- **cluster_stats_pattern**: The pattern used to match cluster statistics files within the input directory. Default is `"split_cluster_stats.tsv"`.

- **min_reads_per_cluster**: Minimum number of reads per cluster to be considered during processing. Default is `10`.

- **max_reads_per_cluster**: Maximum number of reads per cluster to be considered during processing. Default is `200`.

- **max_edit_distance**: Maximum edit distance allowed during merging of haplotype clusters. Default is `2`.

- **use_variant_calling_positions**: A boolean flag indicating whether to use variant calling positions. If `true`, the workflow considers the `variant_calling_positions` file.

- **region_exclusion_ranges**: A map defining region-specific exclusion ranges. Each region can have its own exclusion ranges or none at all. Example:
  ```groovy
  region_exclusion_ranges = [
      "lpa2645": "",              // No positions to exclude
      "lpa5104": "2472,2506"      // Exclude positions 2472-2506
  ]
  ```

- **min_qscore**: The minimum quality score required during processing. Default is `45`.

- **output_format**: The output format for haplotype results. Default is `"fasta"`.

## Other Parameters

- **threads**: The number of threads used during workflow execution. It is set to `(Runtime.runtime.availableProcessors() - 1)` by default.

## Output
- The workflow generates output files, including filtered BAM files, extracted haplotypes, and merged haplotype information.

## Credits

These scripts were originally written for use by [GENEPI](https://genepi.i-med.ac.at/), by ([@StephanAmstler](https://github.com/AmstlerStephan)).
