#!/usr/bin/env python3

import sys
import os
import argparse
import vcf_io
import arg_io
import numpy as np
import tskit

def read_ts(node_file, edge_file, counts_file):
    node_time = np.atleast_1d(np.loadtxt(node_file))
    edge_all = np.loadtxt(edge_file, ndmin=2)
    keep = edge_all[:, 2] >= 0
    edge_span = edge_all[keep, :]
    length = max(edge_span[:, 1])
    tables = tskit.TableCollection(sequence_length=length)
    node_table = tables.nodes
    parents = set(edge_span[:, 2].astype(int))
    for i, t in enumerate(node_time):
        node_table.add_row(flags=0 if i in parents else tskit.NODE_IS_SAMPLE, time=t)
    parent_indices = np.array(edge_span[:, 2], dtype = np.int32)
    child_indices = np.array(edge_span[:, 3], dtype = np.int32)
    arg_io.add_edges(tables, edge_span, parent_indices, child_indices, 0.0, counts_file, keep)
    for i in range(len(parent_indices)):
        j = int(parent_indices[i])
        k = int(child_indices[i])
        if (node_time[j] < node_time[k]):
            print(j, k, node_time[j], node_time[k])
    for i in range(edge_span.shape[0]):
        if (edge_span[i, 0] >= edge_span[i, 1]):
            print(edge_span[i, :])
    return tables

def read_ARG(node_file, branch_file, mutation_file, counts_file=None, bases=False):
    tables = read_ts(node_file, branch_file, counts_file)
    arg_io.add_mutations(tables, mutation_file, tables.sequence_length, 0.0, lambda n: n, bases)
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    return tables


def write_trees(input_prefix, output_prefix, start, end, step, vcf_file=None):
    bases = arg_io.uses_bases(input_prefix)
    for i in range(start, end, step):
        node_file, rescaled = arg_io.nodes_file(input_prefix, i)
        branch_file = f"{input_prefix}_branches_{i}.txt"
        counts_file = f"{input_prefix}_counts_{i}.txt"
        if not os.path.exists(node_file):
            sys.exit(f"No MCMC sample {i}: {node_file} not found. -start and -end are MCMC sample indices, -end exclusive.")
        draws = arg_io.muts_files(input_prefix, i) or [(None, None)]
        for draw, mutation_file in draws:
            tables = read_ARG(node_file, branch_file, mutation_file, counts_file, bases)
            arg_io.add_provenance(tables, 'convert_to_tskit.py', [input_prefix], i, draw, rescaled)
            ts = tables.tree_sequence()
            if vcf_file is not None:
                ts = vcf_io.add_individuals_from_vcf(ts, vcf_file)
            ts.dump(arg_io.trees_name(output_prefix, i, draw))

def write_fast_trees(input_prefix, output_prefix, start, end, step):
    for i in range(start, end, step):
        trees_file = f"{output_prefix}_{i}.trees"
        node_file = f"{input_prefix}_fast_nodes_{i}.txt"
        branch_file = f"{input_prefix}_fast_branches_{i}.txt"
        mutation_file = f"{input_prefix}_fast_muts_{i}.txt"
        ts = read_ARG(node_file, branch_file, mutation_file).tree_sequence()
        ts.dump(trees_file)

def main():
    parser = argparse.ArgumentParser(description='Convert to tskit format')

    parser.add_argument('-input', type=str, required=True, help='Prefix of ARG files.')
    parser.add_argument('-output', type=str, required=True, help='Prefix of output files: <output>_<i>.trees, or <output>_<r>_<i>.trees for draw r of a run with -mapping_rep above 1. The node times are the rescaled ones (<input>_rescaled_nodes_<i>.txt) when the run wrote them, else the sampled ones.')
    parser.add_argument('-start', type=int, required=True, help='Index of the first MCMC sample to convert.')
    parser.add_argument('-end', type=int, required=True, help='Index after the last MCMC sample to convert (exclusive).')
    parser.add_argument('-step', type=int, default=1, help='Step size of subsampling. Default: 1.')
    parser.add_argument('-vcf', type=str, default=None,
                        help='Input file of the SINGER run (.vcf, .vcf.gz or .bcf). Sample names from the header are '
                             'attached as individuals, and tips are named <sample>_0 / <sample>_1.')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    write_trees(args.input, args.output, args.start, args.end, args.step, args.vcf)

if __name__ == '__main__':
    main()
