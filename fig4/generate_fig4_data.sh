#!/bin/bash

############################################################
# compute multi-MUMs/MEMs for chr8 and potato chr3
# using mumemto version 1.1.1
############################################################

DATA="../../../data"

mkdir -p output

# Step 1: run mumemto on chr8, find all multi-MUMs
/usr/bin/time -v mumemto ${DATA}/fasta/hg_chr8/*.fa -o output/hg_chr8
# finds all MEMs
/usr/bin/time -v mumemto ${DATA}/fasta/hg_chr8/*.fa -o output/hg_chr8_mems -f 0 -k 2

# Step 2: run mumemto on potato chr3
/usr/bin/time -v mumemto ${DATA}/fasta/potato_chr3/*.fa -o output/potato_chr3
# finds all MEMs
/usr/bin/time -v mumemto ${DATA}/fasta/potato_chr3/*.fa -o output/potato_chr3_mems -f 0 -k 2


# Compute mem density
python mem_density.py -m output/hg_chr8_mems -n $(wc -l < output/hg_chr8.lengths) -l $(awk 'BEGIN{max=0} {if($2>max) max=$2} END{print max}' output/hg_chr8.lengths)

python mem_density.py -m output/potato_chr3_mems -n $(wc -l < output/potato_chr3.lengths) -l $(awk 'BEGIN{max=0} {if($2>max) max=$2} END{print max}' output/potato_chr3.lengths)
