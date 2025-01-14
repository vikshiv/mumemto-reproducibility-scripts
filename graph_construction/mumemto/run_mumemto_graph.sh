#!/bin/bash

############################################################
# build mumemto-based graph of chr19 haplotypes
# using mumemto version 1.1.1
############################################################

DATA="../../../data"

mkdir -p output

# Generate seqs.txt from fasta files
for fa in ${DATA}/fasta/hg_chr19/*.fa; do
    name=$(basename "$fa" _chr19.fa)
    echo -e "${name}\t${fa}" >> output/chr19_seqs.txt
done

# Step 1: run mumemto
/usr/bin/time -v mumemto ../../../fasta/hg_chr19/*.fa -o output/chr19

# Step 2: build mumemto graph
python mum_to_graph.py -i output/chr19 -c 19 -o output

# Step 3: make SV-only graph
python graph_to_sv_only.py -g output/chr19_full.gfa -f output/chr19_filelist.txt -o output -c 19

# Step 4: run minigraph-cactus on the SV-only graph
/usr/bin/time -v cactus-graphmap ./js output/chr19_seqs.txt output/chr19_sv_only.gfa chr19.paf --reference chm13 --outputFasta chr19.sv.gfa.fa
/usr/bin/time -v cactus-align ./js output/chr19_seqs.txt chr19.paf chr19.hal --reference chm13 --pangenome --outVG
/usr/bin/time -v cactus-graphmap-join ./js --hal chr19.hal --outDir output --outName chr19 --reference chm13 --vg chr19.vg

# Step 5: extract missing contigs from the mumemto + MC graph
cat output/chr19_full.gfa | grep ^W | awk '{print $4 "\t" $5 "\t" $6}' | bedtools sort > output/chr19_full.contigs.bed
bedtools subtract -a output/chr19_lengths.bed -b output/chr19_mgc.contigs.bed > output/chr19_mgc.missing.bed

