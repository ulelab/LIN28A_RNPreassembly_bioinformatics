#!/bin/bash

files=$(ls ../../data/mm39/FlagIp/Crosslinks/AllSamples_merged.bed*)

# Define binding sites based on motif groups

prtxn="../../data/mm39/FlagIp/prtxn/PrtxnFile.tsv"
dir=results
mkdir $dir
while read line
do
    g=$(echo $line | cut -d" " -f2)
    consensus=$(echo $line | cut -d" " -f1)
    echo $consensus
    echo $g
    for f in $files
    do
        echo $f
        # Run bs assignment
        prtxn=$prtxn file=$f g=$g consensus=$consensus outDir=$dir sbatch --output=$dir/${consensus}_FILENAME_${fname}.out --error=$dir/${consensus}_FILENAME_${fname}.err run_bsAssign.sbatch
    done
done < CLUSTER_motifGroupsFromAllKmers.txt
