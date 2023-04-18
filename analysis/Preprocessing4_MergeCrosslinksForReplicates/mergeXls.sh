#!/bin/bash

# PABPC data

fpath="../../data/Pabpc1Pabpc4Iclip_2022/Crosslinks/"
file_ids="DOX_C1_Crick DOX_C1_Lj DOX_C4_Benthyl DOX_C4_Proteintech KO_C1_Crick KO_C1_Lj KO_C4_Benthyl KO_C4_Proteintech"
outpath="../../data/Pabpc1Pabpc4Iclip_2022/Crosslinks/mergedXls/"

# Running in bedtools conda env

for id in $file_ids
do
    echo $id
    files=$(ls ${fpath}${id}*.bed.gz)
    echo $files
    zcat $files | \
    sort -k 1,1 -k2,2n -k3,3n -k6,6 | \
    bedtools groupby -i stdin -g 1,2,3,6 -c 5 -o sum | \
    awk '{OFS="\t"}{print $1,$2,$3,".",$5,$4}' | \
    pigz > ${outpath}${id}_merged.bed.gz
    echo "Merge completed."
done 2>&1 | tee ${outpath}stdout.txt

# LIN28 data

fpath="../../data/LIN28_220626_results/Crosslinks/"
file_ids="LIN28A-WT_ESCiLIF LIN28A-WT_ESC_LIF-CHIR"
outpath="../../data/LIN28_220626_results/Crosslinks/mergedXls/"


for id in $file_ids
do
    echo $id
    files=$(ls ${fpath}${id}*.bed.gz)
    echo $files
    zcat $files | \
    sort -k 1,1 -k2,2n -k3,3n -k6,6 | \
    bedtools groupby -i stdin -g 1,2,3,6 -c 5 -o sum | \
    awk '{OFS="\t"}{print $1,$2,$3,".",$5,$4}' | \
    pigz > ${outpath}${id}_merged.bed.gz
    echo "Merge completed."
done 2>&1 | tee ${outpath}stdout.txt
