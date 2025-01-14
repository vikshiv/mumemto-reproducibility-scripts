#!/bin/bash

############################################################
# compute multi-MUMs of chr17 and 19 to recreate figure 3
# using mumemto version 1.1.1
############################################################

DATA="../../../data"

mkdir -p output

# Step 1: run mumemto on chr17, find partial mums
/usr/bin/time -v mumemto ${DATA}/fasta/hg_chr17/*.fa -o output/chr17_partial -k -1

# Step 2: run mumemto on chr19
/usr/bin/time -v mumemto ${DATA}/fasta/hg_chr19/*.fa -o output/chr19