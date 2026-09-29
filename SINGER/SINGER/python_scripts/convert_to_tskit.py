#!/usr/bin/env python3

import sys
import os
import gzip
import argparse
import numpy as np
import tskit

def read_ts(node_file, edge_file):
    node_time = np.atleast_1d(np.loadtxt(node_file))
    edge_span = np.loadtxt(edge_file, ndmin=2)
    edge_span = edge_span[edge_span[:, 2] >= 0, :]
    length = max(edge_span[:, 1])
    tables = tskit.TableCollection(sequence_length=length)
    node_table = tables.nodes
    edge_table = tables.edges
    parents = set(edge_span[:, 2].astype(int))
    prev_time = -1
    for i, t in enumerate(node_time):
        if i not in parents:
            node_table.add_row(flags=tskit.NODE_IS_SAMPLE, time=t)
        else:
            t = max(prev_time + 1e-4, t)
            node_table.add_row(time = t)
            prev_time = t
    parent_indices = np.array(edge_span[:, 2], dtype = np.int32)
    child_indices = np.array(edge_span[:, 3], dtype = np.int32)
    edge_table.set_columns(left = edge_span[:, 0], right = edge_span[:, 1], parent = parent_indices, child = child_indices)
    for i in range(len(parent_indices)):
        j = int(parent_indices[i])
        k = int(child_indices[i])
        if (node_time[j] < node_time[k]):
            print(j, k, node_time[j], node_time[k])
    for i in range(edge_span.shape[0]):
        if (edge_span[i, 0] >= edge_span[i, 1]):
            print(edge_span[i, :])
    tables.sort()
    return tables

def read_mutation(tables, mutation_file):
    mutations = np.loadtxt(mutation_file, ndmin=2)
    site_ids = {}
    for i in range(mutations.shape[0]):
        mut_pos = mutations[i, 0]
        if mut_pos >= tables.sequence_length:
            continue
        if mut_pos not in site_ids:
            tables.sites.add_row(position=mut_pos, ancestral_state='0')
            site_ids[mut_pos] = tables.sites.num_rows - 1
        tables.mutations.add_row(site=site_ids[mut_pos], node=int(mutations[i, 1]), derived_state=str(int(mutations[i, 3])))
    return

def read_ARG(node_file, branch_file, mutation_file):
    tables = read_ts(node_file, branch_file)
    read_mutation(tables, mutation_file)
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    ts = tables.tree_sequence()    
    return ts


def read_vcf_sample_names(vcf_file):
    for path in (vcf_file, vcf_file + ".vcf.gz", vcf_file + ".vcf"):
        if not os.path.exists(path):
            continue
        with (gzip.open(path, 'rt') if path.endswith(".gz") else open(path)) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    return line.strip().split("\t")[9:]
    raise FileNotFoundError(f"Could not find VCF file: {vcf_file}")

def add_individuals_from_vcf(ts, vcf_file):
    sample_names = read_vcf_sample_names(vcf_file)
    ploidy = ts.num_samples // len(sample_names)
    if ploidy not in (1, 2) or ploidy * len(sample_names) != ts.num_samples:
        raise ValueError(
            f"{ts.num_samples} sample nodes do not match {len(sample_names)} VCF samples at ploidy 1 or 2."
        )
    tables = ts.dump_tables()
    node_individual = tables.nodes.individual.copy()
    node_metadata = [b""] * tables.nodes.num_rows
    for ind_id, name in enumerate(sample_names):
        tables.individuals.add_row(metadata=name.encode())
        for hap in range(ploidy):
            nid = ind_id * ploidy + hap
            node_individual[nid] = ind_id
            node_metadata[nid] = f"{name}_{hap}".encode()
    tables.nodes.individual = node_individual
    tables.nodes.packset_metadata(node_metadata)
    return tables.tree_sequence()

def write_trees(input_prefix, output_prefix, start, end, step, vcf_file=None):
    for i in range(start, end, step):
        trees_file = f"{output_prefix}_{i}.trees"
        node_file = f"{input_prefix}_nodes_{i}.txt"
        branch_file = f"{input_prefix}_branches_{i}.txt"
        mutation_file = f"{input_prefix}_muts_{i}.txt"
        if not os.path.exists(node_file):
            sys.exit(f"No MCMC sample {i}: {node_file} not found. -start and -end are MCMC sample indices, -end exclusive.")
        ts = read_ARG(node_file, branch_file, mutation_file)
        if vcf_file is not None:
            ts = add_individuals_from_vcf(ts, vcf_file)
        ts.dump(trees_file)

def write_fast_trees(input_prefix, output_prefix, start, end, step):
    for i in range(start, end, step):
        trees_file = f"{output_prefix}_{i}.trees"
        node_file = f"{input_prefix}_fast_nodes_{i}.txt"
        branch_file = f"{input_prefix}_fast_branches_{i}.txt"
        mutation_file = f"{input_prefix}_fast_muts_{i}.txt"
        ts = read_ARG(node_file, branch_file, mutation_file)
        ts.dump(trees_file)

def main():
    parser = argparse.ArgumentParser(description='Convert to tskit format')

    parser.add_argument('-input', type=str, required=True, help='Prefix of ARG files.')
    parser.add_argument('-output', type=str, required=True, help='Prefix of output files.')
    parser.add_argument('-start', type=int, required=True, help='Index of the first MCMC sample to convert.')
    parser.add_argument('-end', type=int, required=True, help='Index after the last MCMC sample to convert (exclusive).')
    parser.add_argument('-step', type=int, default=1, help='Step size of subsampling. Default: 1.')
    parser.add_argument('-vcf', type=str, default=None,
                        help='VCF (or prefix without .vcf or .vcf.gz) used as SINGER input. Sample names from the header are '
                             'attached as individuals, and tips are named <sample>_0 / <sample>_1.')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    write_trees(args.input, args.output, args.start, args.end, args.step, args.vcf)

if __name__ == '__main__':
    main()
