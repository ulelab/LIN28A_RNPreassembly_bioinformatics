#!/bin/bash

#SBATCH --job-name=cv_coverage
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --mem=16GB


# Weighted approach
# Unmasked genome
python3 ../../../src/cv_coverage.py \
$file \
$g \
UTR3 \
5 \
/camp/home/kuretk/home/data/genomes/MusMusculus/Goodwright_m39/GRCm39.primary_assembly.genome.fa \
/camp/home/kuretk/home/data/genomes/MusMusculus/Goodwright_m39/GRCm39.primary_assembly.genome.fa.fai \
/camp/home/kuretk/home/data/genomes/MusMusculus/Goodwright_m39/sorted.annotated.regions.gtf \
6 \
None \
300 \
True \
8 \
20000 \
20
