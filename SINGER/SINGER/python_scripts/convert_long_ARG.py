import argparse
import os
import sys
import numpy as np
import tskit
import vcf_io
import arg_io

def read_long_ARG(node_files, branch_files, mutation_files, block_coordinates, counts_files=None, bases=False):
    if len(node_files) != len(branch_files):
        raise ValueError("Lengths of node_files and branch_files must be the same.")

    if len(node_files) != len(block_coordinates):
        raise ValueError("Lengths of node_files and coordinates must be the same.")
    if counts_files is None:
        counts_files = [None] * len(node_files)

    tables = tskit.TableCollection(sequence_length=0)
    node_table = tables.nodes

    node_num = 0
    sample_num = 0

    for node_file_index, (node_file, branch_file, mutation_file, counts_file) in enumerate(zip(node_files, branch_files, mutation_files, counts_files)):
        print(f"Processing segment {node_file_index}")
        node_time = np.atleast_1d(np.loadtxt(node_file))
        edge_all = np.loadtxt(branch_file, ndmin=2)
        keep = edge_all[:, 2] >= 0
        edge_span = edge_all[keep, :]
        if node_file_index == 0:
            sample_num = len(node_time) - len(set(edge_span[:, 2].astype(int)))
            for t in node_time[:sample_num]:
                node_table.add_row(flags=tskit.NODE_IS_SAMPLE, time=t)
        node_num = node_table.num_rows - sample_num
        for t in node_time[sample_num:]:
            node_table.add_row(time=t)

        length = max(edge_span[:, 1])
        tables.sequence_length = max(tables.sequence_length, length + block_coordinates[node_file_index])

        parent_indices = np.array(edge_span[:, 2], dtype=np.int32)
        child_indices = np.array(edge_span[:, 3], dtype=np.int32)

        parent_indices[parent_indices >= sample_num] += node_num
        child_indices[child_indices >= sample_num] += node_num

        arg_io.add_edges(tables, edge_span, parent_indices, child_indices, block_coordinates[node_file_index], counts_file, keep)
        arg_io.add_mutations(tables, mutation_file, block_coordinates[node_file_index],
                             lambda n, s=sample_num, k=node_num: n + k if n >= s else n, bases)

    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    return tables

def block_prefixes(vcf_prefix, output_prefix):
    # Read block_coordinates from the index file
    with open(f"{vcf_prefix}.index", 'r') as f:
        rows = [line.split() for line in f if line.strip()]
    block_coordinates = [int(float(w[-2])) for w in rows if w[0] == rows[0][0]]
    prefixes = [f"{output_prefix}_{i}_{i+1}" for i in range(len(block_coordinates))]
    return prefixes, block_coordinates

def generate_file_lists(prefixes, MCMC_iteration, draw, with_muts=True):
    nodes = [arg_io.nodes_file(p, MCMC_iteration) for p in prefixes]
    node_files = [n[0] for n in nodes]
    rescaled = all(n[1] for n in nodes)
    if any(n[1] for n in nodes) and not rescaled:
        raise SystemExit("Some blocks have rescaled node times and some do not: " + ", ".join(f for f, _ in nodes))
    branch_files = [f"{p}_branches_{MCMC_iteration}.txt" for p in prefixes]
    if not with_muts:
        mutation_files = [None for p in prefixes]
    else:
        mutation_files = [f"{p}_muts_{MCMC_iteration}.txt" if draw is None else f"{p}_muts_{draw}_{MCMC_iteration}.txt" for p in prefixes]
    counts_files = [f"{p}_counts_{MCMC_iteration}.txt" for p in prefixes]
    return node_files, branch_files, mutation_files, counts_files, rescaled

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Generate tskit format for a long ARG.")

    # Add arguments with prefixes
    parser.add_argument("-vcf", required=True, help="Input file of the SINGER run (.vcf, .vcf.gz or .bcf)")
    parser.add_argument("-output", required=True, help="Prefix of the block runs (<output>_<b>_<b+1>) and of the output: <output>_<iteration>.trees, or <output>_<r>_<iteration>.trees for draw r of a run with -mapping_rep above 1. The node times are the rescaled ones when the runs wrote them, else the sampled ones.")
    parser.add_argument("-iteration", type=int, required=True, help="MCMC iteration for generating filenames")

    args = parser.parse_args()

    prefixes, block_coordinates = block_prefixes(args.vcf, args.output)
    draws = arg_io.muts_files(prefixes[0], args.iteration)
    with_muts = bool(draws)
    if not with_muts:
        print(f"Warning: no muts file for iteration {args.iteration}; the tree sequence is written without sites.", file=sys.stderr)
    bases = arg_io.uses_bases(prefixes[0])
    for draw, _ in draws or [(None, None)]:
        node_files, branch_files, mutation_files, counts_files, rescaled = generate_file_lists(prefixes, args.iteration, draw, with_muts)
        tables = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates, counts_files, bases)
        arg_io.add_provenance(tables, 'convert_long_ARG.py', prefixes, args.iteration, draw, rescaled)
        ts = tables.tree_sequence()
        ts = vcf_io.add_individuals_from_vcf(ts, args.vcf)
        name = arg_io.trees_name(args.output, args.iteration, draw)
        print(f"Save to {name}")
        ts.dump(name)

if __name__  == "__main__":
    main()
