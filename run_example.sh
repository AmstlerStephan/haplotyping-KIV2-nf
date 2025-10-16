#!/bin/bash

# Example command to run the modernized haplotyping pipeline for lpa5104 region

nextflow run main.nf \
  --input data \
  --region lpa5104 \
  --output results_lpa5104 \
  --use_variant_calling_positions false \
  --min_reads_per_cluster 10 \
  --max_reads_per_cluster 200 \
  --max_edit_distance 2 \
  --min_qscore 45 \
  --threads 8 \
  -profile docker

# Note: The pipeline will automatically use the region-specific exclusion ranges
# for lpa5104: positions 2472-2506 will be excluded as configured in nextflow.config