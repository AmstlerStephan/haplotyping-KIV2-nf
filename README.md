[<img width="200" align="right" src="docs/images/ecseq.jpg">](https://www.ecseq.com)
[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/ecseq/dnaseq.svg)](https://hub.docker.com/r/ecseq/dnaseq)

umi-pipeline-nf Pipeline
======================

**haplotyping-KIV2-nf** is a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

## Overview
`haplotyping-KIV2-nf` extracts the polymorphic positions of aligned reads provided in bam-file format. It removes the soft- and hard-clips of the sequences, performs a multiple sequence alignment and extracts the polymorphic positions (haplotype) per sample.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/)

2. Start running your own analysis!
2.1 Download and adapt the config/custom.config with paths to your data (relative and absolute paths possible)

```bash
nextflow run AmstlerStephan/haplotyping-KIV2-nf -r main -c <custom.config> -profile <docker/conda> 
```

### Credits

These scripts were originally written for use by [GENEPI](https://genepi.i-med.ac.at/), by ([@StephanAmstler](https://github.com/AmstlerStephan)).
