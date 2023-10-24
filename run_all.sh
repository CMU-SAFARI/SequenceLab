#!/bin/bash

run_sl(){
    TSS_FILE=$1
    READ_COUNT=$2
    EXPERIMENT_NAME=$3
    MIN_EDIT_DISTANCE=$4
    MAX_EDIT_DISTANCE=$5
    MAXLENGTH=$6
    ./main --debug-mode 0 \
        --grid-size 4 \
        --kmer-length 5 \
        --maxlength $MAXLENGTH \
	--file $TSS_FILE \
        --read-count $READ_COUNT \
        --iteration-count 1 \
        --grid-file-prefix "${EXPERIMENT_NAME}_grid" \
        --histogram-file-prefix "${EXPERIMENT_NAME}_hist" \
        --runtime-file "${EXPERIMENT_NAME}_runtime.csv" \
        --min-edit-distance $MIN_EDIT_DISTANCE \
        --max-edit-distance $MAX_EDIT_DISTANCE \
	--activate-individually 0000010000000
}
#1111111111111

run_sl testdata/ont/ONT_mapped_10k_top.tss 1000 ONT_mapped_top 0 20 450000
run_sl testdata/ont/ONT_mapped_10k_bottom.tss 9000 ONT_mapped_bottom 0 20 95000
run_sl testdata/ont/ONT_chained_10k_top.tss 1000 ONT_chained_top 0 20 140000
run_sl testdata/ont/ONT_chained_10k_bottom.tss 9000 ONT_chained_bottom 0 20 4000
run_sl testdata/hifi/HiFi_mapped_100k_top.tss 10000 HiFi_mapped_top 0 10 40000
run_sl testdata/hifi/HiFi_mapped_100k_bottom.tss 90000 HiFi_mapped_bottom 0 10 2200
run_sl testdata/hifi/HiFi_chained_100k_top.tss 9999 HiFi_chained_top 0 10 35000
run_sl testdata/hifi/HiFi_chained_100k_bottom.tss 90001 HiFi_chained_bottom 0 10 2500
run_sl testdata/illumina/Illumina_mapped_10M_top.tss 874877 Illumina_mapped_top 0 10 400
run_sl testdata/illumina/Illumina_mapped_10M_bottom.tss 9125123 Illumina_mapped_bottom 0 10 185
run_sl testdata/illumina/Illumina_chained_10M_top.tss 926785 Illumina_chained_top 0 10 400
run_sl testdata/illumina/Illumina_chained_10M_bottom.tss 9073215 Illumina_chained_bottom 0 10 180
