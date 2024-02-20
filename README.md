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
nextflow run AmstlerStephan/haplotyping-KIV2-nf -r main -c <custom.config> -profile <docker/conda> 
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

## Output
- The workflow generates output files, including filtered BAM files, extracted haplotypes, and merged haplotype information.

## Credits

These scripts were originally written for use by [GENEPI](https://genepi.i-med.ac.at/), by ([@StephanAmstler](https://github.com/AmstlerStephan)).
