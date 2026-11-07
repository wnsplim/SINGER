import tskit
import numpy as np
import pandas as pd
import argparse

def compute_pairwise_coalescence_time(ts, leaf_index, query_index, s):
    windows = np.arange(0, ts.sequence_length, s)
    windows = np.append(windows, ts.sequence_length)
    return ts.diversity(sample_sets=[leaf_index, query_index], windows=windows, mode='branch')/2 + (ts.node(leaf_index).time + ts.node(query_index).time)/2

def compute_all_pairwise_coalescence_times(ts, leaf_index, s):
    windows = np.arange(0, ts.sequence_length, s)
    n = len(windows)
    m = ts.num_samples - 1
    all_coalescence_times = np.zeros((n, m))
    index = 0
    for i in range(0, ts.num_samples):
        if i != leaf_index:
            print(i)
            all_coalescence_times[:, index] = compute_pairwise_coalescence_time(ts, leaf_index, i, s)
            index += 1
    return all_coalescence_times

def main():
    parser = argparse.ArgumentParser(description="Calculate pairwise coalescence times for a given leaf node.")
    parser.add_argument("--trees_file", type=str, required=True, help="Path to the tree sequence file.")
    parser.add_argument("--leaf_index", type=int, required=True, help="Index of the leaf node.")
    parser.add_argument("--interval_size", type=int, required=True, help="Size of the interval.")
    parser.add_argument("--output_file", type=str, required=True, help="Output filename.")
    
    args = parser.parse_args()
    
    ts = tskit.load(args.trees_file)
    all_coalescence_times = compute_all_pairwise_coalescence_times(ts, args.leaf_index, args.interval_size)
    np.savetxt(args.output_file, all_coalescence_times, delimiter=",")
    
 
if __name__ == "__main__":
    main()
