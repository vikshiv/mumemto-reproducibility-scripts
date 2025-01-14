import numpy as np
from tqdm.auto import tqdm
import os
from Bio import SeqIO
import argparse
import re

def parse_arguments():    
    parser = argparse.ArgumentParser(description="builds an SV-only minigraph-like graph from a mumemto-based GFA file")
    parser.add_argument('--gfa','-g', dest='gfa', help='gfa file')
    parser.add_argument('--chr', '-c', dest='chr', help='chromosome', required=True)
    parser.add_argument('--filelist', '-f', dest='filelist', help='list of paths for fasta files')
    parser.add_argument('--output', '-o', dest='output', help='output dir')
    args = parser.parse_args()
    return args

def write_rfga(graphs, out, id_map):
    ### graphs is a list of oblocks as graphs
    ### format: dict of node_name -> node sequence
    ### format: list of links: node1, strand, node2, strand
    nodes, links = graphs
    for n, (s, sn, so) in nodes.items():
        idm = id_map[sn.split('|')[-1]]
        out.write(f'S\t{n}\t{s}\tLN:i:{len(s)}\tSN:Z:id={sn}\tSO:i:{so}\tSR:i:{idm}\n')
    for l in links:
        l = '\t'.join(l)
        out.write(f'L\t{l}\t0M\n')
        
def make_sv_graph(path_id, decorated_mums, gap_names, all_mum_nodes, all_mum_links, chm13_nodes, allele_count, allele_length, len_range, id_map):
    candidate_gaps = [k for k, x in gap_names.items() if len(x) <= allele_count and 
                                                        all([len(all_mum_nodes[y]) < allele_length[1] and len(all_mum_nodes[y]) > allele_length[0] for y in x]) and 
                                                        (max([len(all_mum_nodes[y]) for y in x]) - min([len(all_mum_nodes[y]) for y in x]) > len_range)]
    sv_nodes = chm13_nodes + sum([gap_names[k] for k in candidate_gaps], [])
    sv_nodes_set = set(sv_nodes)
    sv_links = [l for l in all_mum_links if l[0] in sv_nodes_set and l[2] in sv_nodes_set]
    ## write it the minigraph way with nodes
    os.makedirs(args.output, exist_ok=True)
    node_name_map = {n : 's'+str(idx) for idx, n in enumerate(sv_nodes)}
    out = open(os.path.join(args.output, 'mumemto_sv_only.gfa'), 'w')
    write_rfga(({node_name_map[n] : decorated_mums[n] for n in sv_nodes}, [[node_name_map[l[0]], l[1], node_name_map[l[2]], l[3]] for l in sv_links]), out, id_map)
    out.close()
    print(f'length of graph chr{path_id}: ', sum([len(all_mum_nodes[n]) for n in sv_nodes_set]))

def main(args):    
    seqs = [list(SeqIO.parse(f.split()[0], 'fasta'))[0] for f in open(args.filelist, 'r')]
    ids = [s.id for s in seqs]
    id_map = {i : idx for idx, i in enumerate(ids)}
    ### build make shift "SV" graph
    graph = open(args.gfa, 'r')
    all_mum_nodes, all_mum_walks, all_mum_links = [], [], []
    for l in graph:
        l = l.strip().split()
        if l[0] == 'S':
            all_mum_nodes.append(l[1:])
        elif l[0] == 'W':
            # l[6] = [re.findall(r'([<>])([^<>]*)', l[6])]
            all_mum_walks.append(l[1:])
        elif l[0] == 'L':
            all_mum_links.append(l[1:])
    chm13_walks = [w for w in all_mum_walks if w[0] == 'CHM13']
    all_mum_nodes = {k : v for k,v in all_mum_nodes}
    for w in tqdm(all_mum_walks):
        w[5] = re.findall(r'([<>])([^<>]*)', w[5])
    
    ### look at graph and find the "SVs"
    gap_names = {}
    interblock_names = {}
    for k in all_mum_nodes.keys():
        name = k
        k = name.split('_')
        if (len(k) == 3 and k[2][0] == 'a'):
            pref = '_'.join(k[:2])
            if pref in gap_names:
                gap_names[pref].append(name)
            else:
                gap_names[pref] = [name]
        elif (len(k) == 2 and k[1][0] == 'a'):
            pref = k[0]
            if pref in interblock_names:
                interblock_names[pref].append(name)
            else:
                interblock_names[pref] = [name]
              
    ### try reformatting what exists to fit minigraph
    decorated_mums = {}
    for w in chm13_walks:
        offset = int(w[3])
        for _, node in w[5]:
            decorated_mums[node] = (all_mum_nodes[node], 'chm13|chr'+args.chr, offset)
            offset += len(all_mum_nodes[node])
    for w in all_mum_walks:
        if w[0] == 'CHM13':
            continue
        offset = int(w[3])
        for _, node in w[5]:
            if node not in decorated_mums:
                decorated_mums[node] = (all_mum_nodes[node], f'{w[0]}.{w[1]}|chr{args.chr}_{w[0]}.{w[1]}', offset)
            offset += len(all_mum_nodes[node])
    chm13_nodes = sum([[n for _, n in w[5]] for w in chm13_walks], [])
    make_sv_graph(args.chr, decorated_mums, gap_names, all_mum_nodes, all_mum_links, chm13_nodes, 45, (50, 100000), 1000, id_map)
    
if __name__ == '__main__':
    args = parse_arguments()
    main(args)