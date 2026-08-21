#!/bin/bash

vcftools \
  --vcf /path/to/snps_beforeLD.vcf \
  --SNPdensity 50000 \
  --out /path/to/output/directory/<file_name>
  