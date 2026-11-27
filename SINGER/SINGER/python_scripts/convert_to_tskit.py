#!/usr/bin/env python3

import sys
import os
import argparse
import vcf_io
import arg_io


def read_blocks(input_prefix):
    """The block prefixes, their offsets and the region end from <input_prefix>.blocks (first row 'start end', then
    one row 'index left right' per block); a run without that file is one block at offset 0 with no region end."""
    blocks_file = f"{input_prefix}.blocks"
    if not os.path.exists(blocks_file):
        return [input_prefix], [0.0], 0
    with open(blocks_file) as f:
        rows = [line.split() for line in f if line.strip()]
    return [f"{input_prefix}_block{w[0]}" for w in rows[1:]], [float(w[1]) for w in rows[1:]], float(rows[0][1])


def write_trees(input_prefix, output_prefix, start, end, step, vcf_file=None):
    prefixes, offsets, length = read_blocks(input_prefix)
    for i in range(start, end, step):
        nodes = [arg_io.nodes_file(p, i) for p in prefixes]
        for node_file, _ in nodes:
            if not os.path.exists(node_file):
                sys.exit(f"No MCMC sample {i}: {node_file} not found. -start and -end are MCMC sample indices, -end exclusive.")
        rescaled = all(r for _, r in nodes)
        if any(r for _, r in nodes) and not rescaled:
            sys.exit("Some blocks have rescaled node times and some do not: " + ", ".join(f for f, _ in nodes))
        draws = arg_io.muts_files(prefixes[0], i)
        if not draws:
            print(f"Warning: no muts file for sample {i}; the tree sequence is written without sites.", file=sys.stderr)
        for draw, _ in draws or [(None, None)]:
            if not draws:
                mutation_files = [None for p in prefixes]
            else:
                mutation_files = [f"{p}_muts_{i}.txt" if draw is None else f"{p}_muts_{draw}_{i}.txt" for p in prefixes]
            tables = arg_io.read_long_ARG([f for f, _ in nodes], [f"{p}_branches_{i}.txt" for p in prefixes], mutation_files,
                                          offsets, [f"{p}_counts_{i}.txt" for p in prefixes], length)
            arg_io.add_provenance(tables, 'convert_to_tskit.py', prefixes, i, draw, rescaled)
            ts = tables.tree_sequence()
            if vcf_file is not None:
                ts = vcf_io.add_individuals_from_vcf(ts, vcf_file)
            name = arg_io.trees_name(output_prefix, i, draw)
            ts.dump(name)
            print(f"wrote {name}")

def main():
    parser = argparse.ArgumentParser(description='Convert to tskit format')

    parser.add_argument('-input', type=str, required=True, help='Prefix of ARG files. When <input>.blocks exists (a run of singer_master with -L or -blocks), the blocks it lists are stitched into one tree sequence at their chromosome positions, with no edge between them and the region end of its first row as the sequence length; delete a row to leave that block out.')
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
