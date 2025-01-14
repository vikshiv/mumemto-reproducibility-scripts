import numpy as np
from tqdm.auto import tqdm
from numba import njit
import argparse
import os
import pandas as pd

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Aggregates MEM density for downstream plotting")
    parser.add_argument('--mems', '-m', dest='mems', help='path to mems from mumemto', required=True)
    parser.add_argument('--num-seqs', '-n', dest='num', help='num of sequences', required=True, type=int)
    parser.add_argument('--max-length', '-l', dest='size', help='size of sequence window', required=True, type=int)
    args = parser.parse_args()
    return args

@njit
def update_coverage(coverage, idx_list, start_list, l):
    for i in range(len(idx_list)):
        for start, idx in zip(start_list[i], idx_list[i]):
            coverage[idx, start:start+l] += 1

def convert_to_numpy(x):
    return np.fromstring(x, sep=',')

chunksize = 10 ** 6  

def main(args):
    # file = open(args.mems, 'r')
    chunks = pd.read_csv(args.mems, sep='\t', chunksize=chunksize, usecols=[0,1,2], converters={1: convert_to_numpy, 2: convert_to_numpy})
    coverage = np.zeros((args.num, args.size))
    for chunk in chunks:
        update_coverage(coverage, chunk[1], chunk[2], chunk[0])
    # for m in tqdm(file):
    #     m = m.strip().split()
    #     l = int(m[0])
    #     update_coverage(coverage, m[2], m[2], l)    
    filename = os.path.splitext(args.mems)[0]
    coverage_file = filename + '_coverage.npy'
    np.save(coverage_file, coverage)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
