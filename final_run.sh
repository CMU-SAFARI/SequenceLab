#!/bin/bash
srun bash -c "./main --debug-mode 0 \
        --grid-size 4 \
        --kmer-length 5 \
        --maxlength $6 \
	--file $1 \
        --read-count $2 \
        --iteration-count 1 \
        --grid-file-prefix $3_grid \
        --histogram-file-prefix $3_hist \
        --runtime-file $3_runtime.csv \
        --min-edit-distance $4 \
        --max-edit-distance $5 \
	--activate-individually 00000000010000 \
        --triplex-pattern-length 6 \
        --triplex-random \
        --triplex-random-pattern-count 256"
