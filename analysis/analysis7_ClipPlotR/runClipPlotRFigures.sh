#!/bin/bash

# Running in clipplotr environment
pXl='../../data/LIN28_220626_results/Crosslinks'
pPABC1='../../data/Pabpc1Pabpc4Iclip_2022/Crosslinks'
bsPath='../../data/BindingSites/merged0'
bsPathPa='../../data/BindingSites/AAA_sitesPabpcCrosslinks'
AllPaSites='../../results/FindMotifs/motifs_AAA_merged.bed'
pathPaSignal='../../results/FindMotifs/motifs_AATAAA.bed'
annotation='../../data/genomes/Goodwright_m39/post_filtering.gencode.vM28.primary_assembly.annotation.gtf'

xls=${pXl}/LIN28A-WT_ESC_LIF-CHIR-FGF0220626_MM_1.bed.gz,\
${pXl}/LIN28A-WT_ESC_LIF-CHIR-FGF0220626_MM_2.bed.gz,\
${pXl}/LIN28A-WT_ESCiLIF0220626_MM_1.bed.gz,\
${pXl}/LIN28A-WT_ESCiLIF0220626_MM_2.bed.gz,\
${pXl}/LIN28A-WT_ESCiLIF-OLD0220626_MM.bed.gz,\
${pXl}/LIN28A-S200A_ESC_LIF-CHIR-FGF0220626_MM_1.bed.gz,\
\
${pPABC1}/DOX_C1_Crick1.bed.gz,\
${pPABC1}/DOX_C1_Crick2.bed.gz,\
${pPABC1}/KO_C1_Crick1.bed.gz,\
${pPABC1}/KO_C1_Crick2.bed.gz,\
\
${pPABC1}/DOX_C4_Proteintech_1.bed.gz,\
${pPABC1}/DOX_C4_Proteintech_2.bed.gz,\
${pPABC1}/KO_C4_Proteintech_1.bed.gz,\
${pPABC1}/KO_C4_Proteintech_2.bed.gz


labels=LIN28A-WT_ESC_LIF-CHIR_1,\
LIN28A-WT_ESC_LIF-CHIR_2,\
LIN28A-WT_ESC_iLIF_1,\
LIN28A-WT_ESC_iLIF_2,\
LIN28A-WT_ESC_iLIF_3,\
LIN28A-S200A_ESC_LIF-CHIR_1,\
\
PABPC1-DOX_1,\
PABPC1-DOX_2,\
PABPC1-KO_1,\
PABPC1-KO_2,\
\
PABPC4-DOX_1,\
PABPC4-DOX_2,\
PABPC4-KO_1,\
PABPC4-KO_2


auxilliaryFiles=\
${bsPath}/merged0_AllSamples_merged_Unfiltered_AGGG-UGGG_Cluster_1.0.bed,\
${bsPath}/merged0_AllSamples_merged_Unfiltered_UGAU_Cluster_2.0.bed,\
${bsPath}/merged0_AllSamples_merged_Unfiltered_UUUU_Cluster_3.0.bed,\
${pathPaSignal},\
${AllPaSites}


auxilliaryLabels='GGGA/GGGU,UGAU,UUUU,AAUAAA,AAA'

dir=figures_libnorm
mkdir $dir


while read line
do
    name=$(echo $line | cut -d, -f1)
    locus=$(echo $line | cut -d, -f2)
    rbp=$(echo $line | cut -d_ -f2)

    echo ${name}
    echo ${locus}
    echo ${rbp}

    ~/dev/clipplotr/clipplotr \
    -x $xls \
    -l $labels \
    -c '#2a9d8f,#2a9d8f,#f4a261,#f4a261,#f4a261,#e76f51,#2a9d8f,#2a9d8f,#264653,#264653,#48cae4,#48cae4,#023e8a,#023e8a' \
    --groups 'LIN28A,LIN28A,LIN28A,LIN28A,LIN28A,LIN28A,PABPC1,PABPC1,PABPC1,PABPC1,PABPC4,PABPC4,PABPC4,PABPC4' \
    -n libsize \
    -s rollmean \
    -w 50 \
    -g $annotation \
    -r $locus \
    -a transcript \
    -y $auxilliaryFiles \
    --auxiliary_labels $auxilliaryLabels \
    --ratios '1,0.3,0.2,0.2' \
    --size_x 310 \
    --size_y 220 \
    --scale_y \
    -o ${dir}/combined_LIN28ATracks/${name}_${locus}_fullBindingSites.pdf


    ~/dev/clipplotr/clipplotr \
    -x $xls \
    -l $labels \
    -c '#2a9d8f,#2a9d8f,#f4a261,#f4a261,#f4a261,#e76f51,#2a9d8f,#2a9d8f,#264653,#264653,#48cae4,#48cae4,#023e8a,#023e8a' \
    --groups 'LIN28A_FCL,LIN28A_FCL,LIN28A_2iL,LIN28A_2iL,LIN28A_2iL,LIN28A_S200A,PABPC1,PABPC1,PABPC1,PABPC1,PABPC4,PABPC4,PABPC4,PABPC4' \
    -n libsize \
    -s rollmean \
    -w 50 \
    -g $annotation \
    -r $locus \
    -a transcript \
    -y $auxilliaryFiles \
    --auxiliary_labels $auxilliaryLabels \
    --ratios '1,0.3,0.2,0.2' \
    --size_x 310 \
    --size_y 220 \
    --scale_y \
    -o ${dir}/separateTracks/${name}_${locus}_separateTracks_fullBindingSites.pdf
done < FigureSelection1.txt
