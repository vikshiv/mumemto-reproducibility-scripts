import numpy as np
from tqdm.auto import tqdm
import os
from Bio import SeqIO
import argparse

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    parser.add_argument('--filelist', '-f', dest='filelist', help='if the filelist is provided, then FASTA filenames are used as labels')
    parser.add_argument('--chr', '-c', dest='chr', help='chromosome', required=True)
    parser.add_argument('--output', '-o', dest='output', help='output dir')
    
    args = parser.parse_args()
    if args.mumfile:
        args.prefix = os.path.splitext(args.mumfile)[0]
    else:
        args.mumfile = args.prefix + '.mums'
        
    args.lens = args.prefix + '.lengths'
    args.filelist = args.prefix + '_filelist.txt'
    if not args.chr.startswith('chr'):
        args.chr = 'chr' + args.chr
    return args

def reverse_strand(mum, seq_lengths):
    l, starts, strands = mum
    new_starts = np.array([p if s == '+' else seq_lengths[idx] - p - l for idx, (p, s) in enumerate(zip(starts, strands))])
    return [l, new_starts, strands]

def get_hap_ids(args):
    hap_ids = [os.path.basename(l.split()[0]).split('_' + args.chr)[0] for l in open(args.lens, 'r').read().splitlines()]
    hap_ids = [tuple(h.split('.') + [args.chr]) if h != 'chm13' else ('CHM13', '1', args.chr) for h in hap_ids]
    return hap_ids
def get_mums(args, seq_lengths):
    mums = [l.split() for l in open(args.mumfile, 'r').read().splitlines()]
    mums = [[int(l[0]), np.fromstring(l[1], dtype=int, sep=','), tuple(l[2].split(','))] for l in mums]
    # mums = [reverse_strand(mum, seq_lengths) for mum in mums]
    mums = sorted(mums, key=lambda x: x[1][0])
    return mums

def find_common_mum_gaps(mums):
    starts = np.array([m[1] for m in mums])
    mum_orders = starts.transpose().argsort()
    mum_gaps = []
    flips = set([])
    for i in tqdm(range(mum_orders.shape[0])):
        cur = []
        for l in range(1, mum_orders.shape[1]):
            left, right = mum_orders[i][l-1], mum_orders[i][l]
            if mums[left][2][i] == mums[right][2][i]:
                if mums[left][2][i] == '+':
                    cur.append((left, right))
                else:
                    cur.append((right, left))
                    flips.add((right, left))
        mum_gaps.append(cur)
    common_gaps = set.intersection(*map(set, mum_gaps))
    left, right = zip(*common_gaps)
    return common_gaps, mum_orders

def get_gap_lengths(gap_order, mums):
    gap_lengths = []
    for l, r in gap_order:
        lens = np.full(len(mums[l][1]), mums[l][0])
        lens[(mums[r][1] < mums[l][1])] = mums[r][0] 
        gap_lengths.append(np.abs(mums[l][1] - mums[r][1]) - lens)
    gap_lengths = np.array(gap_lengths)
    return gap_lengths

def trim_mums(gap_lengths, common_gap_order, mums):
    ### trim mums, returns void
    to_trim = np.unique(np.where(gap_lengths < 0)[0])
    for idx in to_trim:
        l,r = common_gap_order[idx]
        gap_lens = gap_lengths[idx]
        trim_bases = -gap_lens.min()
        if mums[l][0] > mums[r][0]:
            mums[l][0] -= trim_bases
            if '-' in mums[l][2]:
                for i, s in enumerate(mums[l][2]):
                    if s=='-':
                        mums[l][1][i] += trim_bases
        else:
            mums[r][0] -= trim_bases
            for i, s in enumerate(mums[l][2]):
                if s=='+':
                    mums[r][1][i] += trim_bases

def get_oblocks(common_gap_order, mums):
    ### first "oriented blocks"
    oblocks = []
    last_l, last_r = common_gap_order[0]
    cur_block = [last_l, last_r]
    for i in range(1, len(common_gap_order)):
        cur_l, cur_r = common_gap_order[i]
        if cur_l == last_r and mums[cur_l][2] == mums[last_r][2]:
            cur_block.append(cur_r)
        else:
            oblocks.append(cur_block)
            cur_block = [cur_l, cur_r]
        last_l, last_r = cur_l, cur_r
    return oblocks

def trim_oblocks(oblocks, mums):
    ### look at oblock coverage and oblock gaps
    # trim the oblock overlaps
    for i in range(1, len(oblocks)):
        l,r = oblocks[i-1][-1], oblocks[i][0]
        lens = np.full(len(mums[l][1]), mums[l][0])
        lens[(mums[r][1] < mums[l][1])] = mums[r][0] 
        gap_lens = np.abs(mums[l][1] - mums[r][1]) - lens
        if not np.any(gap_lens < 0):
            continue
        trim_bases = -gap_lens.min()
        if mums[l][0] > mums[r][0]:
            mums[l][0] -= trim_bases
            if '-' in mums[l][2]:
                for i, s in enumerate(mums[l][2]):
                    if s=='-':
                        mums[l][1][i] += trim_bases
        else:
            mums[r][0] -= trim_bases
            for i, s in enumerate(mums[l][2]):
                if s=='+':
                    mums[r][1][i] += trim_bases

### graph aux functions
def get_node_seq(seqs, i, start, end, forward = True):
        return seqs[i][start:end] if forward else seqs[i][start:end].reverse_complement()
    
def link_oriented_block(sequences, oblock_id, oblock, mums):
    def add_mum(l):
        return get_node_seq(sequences, 0, mums[l][1][0], mums[l][1][0] + mums[l][0], True)
    def add_intercol_mum(l, r, oblock_id, gap_id):
        cur_seqs = {}
        hap_links = []
        empties = set()
        for idx in range(NUM_SEQS):
            if mums[l][2][idx] == '+':
                seq = get_node_seq(sequences, idx, mums[l][1][idx] + mums[l][0], mums[r][1][idx], True)
            else:
                seq = get_node_seq(sequences, idx, mums[r][1][idx] + mums[r][0], mums[l][1][idx], False)
            if seq not in cur_seqs.keys():
                cur_seqs[seq] = "oblock%d_g%d_a%d"%(oblock_id, gap_id, len(cur_seqs))
            if seq:
                hap_links.append(cur_seqs[seq])
            else:
                hap_links.append(None)
                empties.add(cur_seqs[seq])
        return cur_seqs, hap_links, empties
    strands = mums[oblock[0]][2]
    NUM_SEQS = len(strands)
    hap_walks = [[] for _ in range(len(strands))]
    nodes = {}
    node_order = []
    all_empties = set()
    for m in range(len(oblock) - 1):
        mum_name = "oblock%d_m%d"%(oblock_id, oblock[m])
        node_order.append(mum_name)
        nodes[mum_name] = add_mum(oblock[m])
        for i in range(len(hap_walks)):
            hap_walks[i].append(mum_name)
        im_seqs, hap_path, empties = add_intercol_mum(oblock[m], oblock[m] + 1, oblock_id, m)
        for s, n in im_seqs.items():
            if n not in empties:
                nodes[n] = s
        node_order.append(list(im_seqs.values()))
        for idx, h in enumerate(hap_path):
            if h != None:
                hap_walks[idx].append(h)
        all_empties.update(empties)
    mum_name = "oblock%d_m%d"%(oblock_id, oblock[-1])
    node_order.append(mum_name)
    nodes[mum_name] = add_mum(oblock[-1])
    for i in range(len(hap_walks)):
        hap_walks[i].append(mum_name)
    ### with nodes, connect them in order to build links
    links = []
    for i in range(0, (len(node_order) // 2) - 1, 2):
        mum1 = node_order[i]
        ims = node_order[i+1]
        mum2 = node_order[i+2]
        for im_node in ims:
            if im_node in all_empties:
                links.append((mum1, '+', mum2, '+'))
            else:
                links.append((mum1, '+', im_node, '+'))
                links.append((im_node, '+', mum2, '+'))
            
    ### if one of the haplotypes is inverted
    if '-' in strands:
        # first add the reverse links
        for i in range(0, (len(node_order) // 2) - 1, 2):
            mum1 = node_order[i]
            ims = node_order[i+1]
            mum2 = node_order[i+2]
            for im_node in ims:
                if im_node in all_empties:
                    links.append((mum2, '-', mum1, '-'))
                else:
                    links.append((mum2, '-', im_node, '-'))
                    links.append((im_node, '-', mum1, '-')) ## invert the IM regions in case they collapse with another allele, so they must be traversed forward
        
    ### write the haplotype walks
    walks = []
    for idx, hap in enumerate(hap_walks):
        if strands[idx] == '+':
            walks.append([idx, mums[oblock[0]][1][idx], mums[oblock[-1]][1][idx] + mums[oblock[-1]][0], [(h, '>') for h in hap]])
        else:
            walks.append([idx, mums[oblock[-1]][1][idx], mums[oblock[0]][1][idx]+ mums[oblock[0]][0], [(h, '<') for h in hap][::-1]])
    return nodes, links, walks

def write_gfa(graphs, hap_ids, out):
    ### graphs is a list of oblocks as graphs
    ### format: dict of node_name -> node sequence
    ### format: list of links: node1, strand, node2, strand
    ### format: list of walks, idx, start, end, walk nodes (with strand)
    nodes, links, walks = graphs
    for n, s in nodes.items():
        out.write(f'S\t{n}\t{s}\n')
    for l in links:
        out.write(f'L\t{'\t'.join(l)}\n')
    for idx, s, e, walk in walks:
        out.write(f'W\t{'\t'.join(hap_ids[idx])}\t{s}\t{e}\t{''.join([b+a for a,b in walk])}\n')

     
def write_simple_graph(graphs, hap_ids):
    def write_gfa_links(links, out):
        ### format: list of links: node1, strand, node2, strand
        for l in links:
            out.write(f'L\t{'\t'.join(l)}\n')
    def write_gfa_nodes(nodes, out):
        ### format: dict of node_name -> node sequence
        for n, s in nodes.items():
            out.write(f'S\t{n}\t{s}\n')
    def write_gfa_walks(walks, out):
        ### format: list of walks, idx, start, end, walk nodes (with strand)
        for idx, s, e, walk in walks:
            out.write(f'W\t{'\t'.join(hap_ids[idx])}\t{s}\t{e}\t{''.join([b+a for a,b in walk])}\n')
    all_simple_links = []
    for i in tqdm(range(len(graphs))):
        nodes, _, walks = graphs[i]
        for w in walks:
            w = w[3]
            last = w[0]
            for idx in range(1, len(w)):
                l, r = w[idx]
                all_simple_links.append((last[0], '+' if last[1]=='>' else '-', l, '+' if r=='>' else '-'))
                last = (l,r)
    all_simple_links = set(all_simple_links)

    out = open('simple.gfa', 'w')
    for nodes, _, walks in tqdm(graphs):
        write_gfa_nodes(nodes, out)
        write_gfa_walks(walks, out)
    write_gfa_links(all_simple_links, out)
    out.close()


def get_gap_seqs(sequences, starting_pos, ending_pos, gap_id, rearranged_oblocks = set(), max_gap_size=100000):
    seqs = {}
    hap_links = []
    for idx in range(len(sequences)):
        if (gap_id, idx) in rearranged_oblocks or (gap_id-1, idx) in rearranged_oblocks:
            hap_links.append(None)
            continue
        seq = get_node_seq(sequences, idx, starting_pos[idx], ending_pos[idx], True)
        if not seq or ((idx != 0) and (len(seq) > max_gap_size)):
            hap_links.append(None)
            continue
        if seq not in seqs.keys():
            seqs[seq] = "interblock%d_a%d"%(gap_id, len(seqs))
        hap_links.append((seqs[seq], '>'))
    return seqs, hap_links
    
def build_graphs(sequences, graphs, oblocks, seq_lengths, rearranged_oblocks = set(), max_gap_size=100000):
    ### get large gaps too
    # assume oblocks is sorted by coords, along with graphs
    NUM_SEQS = len(seq_lengths)
    all_nodes = {}
    all_walks = []
    all_links = []
    cur_walks = [[i, 0,0,[]] for i in range(NUM_SEQS)]
    for i in range(len(oblocks)):
        temp_nodes, temp_links, temp_walks = graphs[i]
        all_nodes.update(temp_nodes)
        all_links += temp_links
        start_pos = [g[2] for g in cur_walks]
        end_pos = [g[1] for g in temp_walks]
        next_walks = [g[3] for g in temp_walks]
        seqs, hap_links = get_gap_seqs(sequences, start_pos, end_pos, i, rearranged_oblocks=rearranged_oblocks, max_gap_size=max_gap_size)
        for k, v in seqs.items():
            all_nodes[v] = k
            
        for idx, h in enumerate(hap_links):
            if h == None:
                if cur_walks[idx][3]:
                    all_walks.append(cur_walks[idx])
                cur_walks[idx] = temp_walks[idx]
                continue
            ## add links
            all_links.append((h[0], '+', next_walks[idx][0][0], '+' if next_walks[idx][0][1] == '>' else '-'))
            if i != 0:
                all_links.append((cur_walks[idx][3][-1][0], '+' if cur_walks[idx][3][-1][1] == '>' else '-', h[0], '+'))
            cur_walks[idx][3] += [h] + next_walks[idx]
            cur_walks[idx][2] = temp_walks[idx][2]
    ###epilogue, add ending sequence
    start_pos = [g[2] for g in cur_walks]
    end_pos = seq_lengths
    seqs, hap_links = get_gap_seqs(sequences, start_pos, end_pos, len(graphs), rearranged_oblocks=rearranged_oblocks, max_gap_size=max_gap_size)
    for k, v in seqs.items():
        all_nodes[v] = k
    for idx, h in enumerate(hap_links):
        if h == None:
            if cur_walks[idx][3]:
                all_walks.append(cur_walks[idx])
            continue
        ## add links
        all_links.append((h[0], '+', next_walks[idx][0][0], '+' if next_walks[idx][0][1] == '>' else '-'))
        if i != 0:
            all_links.append((cur_walks[idx][3][-1][0], '+' if cur_walks[idx][3][-1][1] == '>' else '-', h[0], '+'))      
        cur_walks[idx][3] += [h]
        cur_walks[idx][2] = end_pos[idx]
        all_walks.append(cur_walks[idx])
    return all_nodes, all_links, all_walks

def find_rearranged_blocks(oblocks, mum_orders):
    rearranged_oblocks = []
    for i in tqdm(range(1, len(oblocks))):
        order = (np.where(mum_orders == oblocks[i][0])[1] - np.where(mum_orders == oblocks[i-1][0])[1])
        if np.any(order < 0):
            for x in np.where(order < 0)[0]:
                rearranged_oblocks.append((i-1, x))
                rearranged_oblocks.append((i, x))
    return rearranged_oblocks

def main(args):
    seq_lengths = np.array([int(l.split()[1]) for l in open(args.lens, 'r').read().splitlines()])
    NUM_SEQS = len(seq_lengths)
    hap_ids = get_hap_ids(args)
    mums = get_mums(args, seq_lengths)
    print('read in mums')
    common_gaps, mum_orders = find_common_mum_gaps(mums)
    common_gap_order = list(sorted(common_gaps))
    gap_lengths = get_gap_lengths(common_gap_order, mums)
    print('found collinear mum gaps')
    ## build naive graph
    seqs = [list(SeqIO.parse(f.split()[0], 'fasta'))[0].seq for f in open(args.filelist, 'r')]
    print('loaded sequences')
    trim_mums(gap_lengths, common_gap_order, mums)
    oblocks = get_oblocks(common_gap_order, mums)
    trim_oblocks(oblocks, mums)
    print('trimmed mums and found oblocks')
    graphs = [link_oriented_block(seqs, i, o, mums) for i, o in tqdm(enumerate(oblocks))]
    print('built oblock graphs')
    # write_simple_graph(graphs, hap_ids)
    # print('wrote simple graph')
    rearranged_oblocks = find_rearranged_blocks(oblocks, mum_orders)
    all_nodes, _, all_walks = build_graphs(seqs, graphs, oblocks, seq_lengths, rearranged_oblocks=set(rearranged_oblocks))
    print('built final connected graph')

    new_all_links = []
    for w in tqdm(all_walks):
        w = w[3]
        last = w[0]
        for idx in range(1, len(w)):
            i, j = w[idx]
            new_all_links.append((last[0], '+' if last[1]=='>' else '-', i, '+' if j=='>' else '-'))
            last = (i,j)

    out = open(os.path.join(args.output, args.chr + '_full.gfa'), 'w')
    out.write(f'H\tVN:Z:1.1\tRS:Z:{hap_ids[0][0]}\n')
    write_gfa((all_nodes, list(set(new_all_links)), all_walks), hap_ids, out)
    out.close()

    print('wrote graph')
    
if __name__ == '__main__':
    args = parse_arguments()
    main(args)