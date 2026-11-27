#!/usr/bin/env python3

import sys
import os
import argparse
import vcf_io
import arg_io


def read_blocks(input_prefix):
    blocks_file = f"{input_prefix}.blocks"
    if not os.path.exists(blocks_file):
        return [input_prefix], [0.0], 0
    with open(blocks_file) as f:
        rows = [line.split() for line in f if line.strip()]
    return [f"{input_prefix}_block{w[0]}" for w in rows[1:]], [float(w[1]) for w in rows[1:]], float(rows[0][1])


def write_trees(input_prefix, output_prefix, start, end, step, vcf_file=None):
    prefixes, offsets, length = read_blocks(input_prefix)
    for i in range(start, end, step):
        files = [arg_io.sample_file(p, i) for p in prefixes]
        for f in files:
            if not os.path.exists(f):
                sys.exit(f"No MCMC sample {i}: {f} not found.")
        samples = [arg_io.read_sample(f) for f in files]
        rescaled = all(s['rescaled'] is not None for s in samples)
        if any(s['rescaled'] is not None for s in samples) and not rescaled:
            sys.exit("Some blocks have rescaled node times and some do not: " + ", ".join(files))
        draws = len(samples[0]['draws'])
        if draws == 0:
            print(f"Warning: sample {i} has no mutation draw.", file=sys.stderr)
        for d in range(draws) if draws else [None]:
            draw = None if draws <= 1 else d + 1
            tables = arg_io.read_long_ARG(files, offsets, length, d)
            arg_io.add_provenance(tables, 'convert_to_tskit.py', prefixes, i, draw, rescaled)
            ts = tables.tree_sequence()
            if vcf_file is not None:
                ts = vcf_io.add_individuals_from_vcf(ts, vcf_file)
            name = arg_io.trees_name(output_prefix, i, draw)
            ts.dump(name)
            print(f"wrote {name}")

def main():
    parser = argparse.ArgumentParser(description='Convert to tskit format')

    parser.add_argument('-input', type=str, required=True, help='The .ckp file of the run, or its .blocks file for a run in blocks.')
    parser.add_argument('-output', type=str, required=True, help='Prefix of output files.')
    parser.add_argument('-start', type=int, required=True, help='Start index of the sample.')
    parser.add_argument('-end', type=int, required=True, help='End index of the sample.')
    parser.add_argument('-step', type=int, default=1, help='Step size of subsampling. Default: 1.')
    parser.add_argument('-vcf', type=str, default=None, help='Input VCF of the run, for the sample names.')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    prefix, extension = os.path.splitext(args.input)
    if extension not in ('.ckp', '.blocks'):
        sys.exit(f"-input {args.input}: not a .ckp or .blocks file.")

    write_trees(prefix, args.output, args.start, args.end, args.step, args.vcf)

if __name__ == '__main__':
    main()
