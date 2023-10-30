#!/bin/bash

#SBATCH --job-name=bs_assign
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --mem=32GB


python3 ../../src/bs_assign/bs_assignement.py \
$g \
$file \
$prtxn \
/camp/home/kuretk/home/data/genomes/MusMusculus/Goodwright_m39/GRCm39.primary_assembly.genome.fa \
/camp/home/kuretk/home/data/genomes/MusMusculus/Goodwright_m39/GRCm39.primary_assembly.genome.fa.fai \
20 \
5 \
8 \
10000 \
$outDir \
$consensus

