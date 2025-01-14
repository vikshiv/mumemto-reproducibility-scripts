#!/bin/bash

############################################################
# recreate minigraph cactus graph of chr19 haplotypes
# using cactus-pangenome version 7.0.0, cactus version 2.9.1
############################################################

DATA="../../../data"

# Generate seqs.txt from fasta files
for fa in ${DATA}/fasta/hg_chr19/*.fa; do
    name=$(basename "$fa" _chr19.fa)
    echo -e "${name}\t${fa}" >> seqs.txt
done

# Step 1: run minigraph-cactus
/usr/bin/time -v cactus-pangenome ./js seqs.txt --outDir output --outName chr19_mgc --reference chm13 --binariesMode local

# Step 2: extract missing contigs from the graph
cat output/chr19_mgc.gfa | grep ^W | awk '{print $4 "\t" $5 "\t" $6}' | bedtools sort > output/chr19_mgc.contigs.bed
bedtools subtract -a output/chr19_lengths.bed -b output/chr19_mgc.contigs.bed > output/chr19_mgc.missing.bed
