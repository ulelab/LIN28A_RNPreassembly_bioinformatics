#!/bin/bash

# Running in JupyterAnalysis conda environment, using pyranges (https://github.com/biocore-NTNU/pyranges) and https://github.com/MuhammedHasan/gencode_utr_fix

# Full annotation
gencode_utr_fix \
--input_gtf '../../data/genomes/Goodwright_m39/gencode.vM28.primary_assembly.annotation.gtf.gz' \
--output_gtf '../../data/genomes/Goodwright_m39/gencode.vM28.primary_assembly.annotation.gencode_utr_fix.gtf.gz'