import tskit
import numpy as np
import argparse

def mse(x, y):
    x = np.array(x)
    y = np.array(y)
    return np.mean((x - y)**2)

def count_incompatibility(ts):
    unmapped_sites = 0
    for tree in ts.trees():
        for site in tree.sites():
            num_mutations = len(site.mutations)
            if num_mutations > 2:
                unmapped_sites += 1
            elif num_mutations == 2:
                if site.mutations[0].node != tree.root and site.mutations[1].node != tree.root:
                    unmapped_sites += 1
    return unmapped_sites

def incompatibility_trace(prefix, indices):
    counts = []
    for index in indices:
        file_name = f"{prefix}_{index}.trees"
        try:
            ts = tskit.load(file_name)
            count = count_incompatibility(ts)
            counts.append(count)
        except FileNotFoundError:
            print(f"File not found: {file_name}")
            counts.append(None)
    return counts

def diversity_fit_mse(ts, m):
    windows = np.append(np.arange(0, ts.sequence_length, 1e6), ts.sequence_length)
    site_diversity = ts.diversity(windows=windows, mode='site')
    branch_diversity = ts.diversity(windows=windows, mode='branch')*m
    fit_mse = mse(site_diversity, branch_diversity)
    return fit_mse

def diversity_fit_trace(prefix, m, indices):
    fit_mses = []
    for index in indices:
        file_name = f"{prefix}_{index}.trees"
        try:
            ts = tskit.load(file_name)
            fit_mse = diversity_fit_mse(ts, m)
            fit_mses.append(fit_mse)
        except FileNotFoundError:
            print(f"File not found: {file_name}")
            fit_mses.append(None)
    return fit_mses


def main():
    parser = argparse.ArgumentParser(description='Compute traces for MCMC samples from SINGER.')

    parser.add_argument('-prefix', type=str, required=True, help='Prefix of the .trees files written by convert_to_tskit (PREFIX_INDEX.trees).')
    parser.add_argument('-m', type=float, required=True, help='Mutation rate.')
    parser.add_argument('-start_index', type=int, default=0, help='Index of the first ARG sample. Default: 0.')
    parser.add_argument('-end_index', type=int, required=True, help='Index after the last ARG sample (exclusive, as in convert_to_tskit -end).')
    parser.add_argument('-step', type=int, default=1, help='Step between ARG sample indices. Default: 1.')
    parser.add_argument('-output_filename', type=str, required=True, help='Output filename of the MCMC traces (tab-separated).')

    args = parser.parse_args()

    indices = range(args.start_index, args.end_index, args.step)
    with open(args.output_filename, 'w') as f:
        f.write('index\tincompatible_sites\tdiversity_fit_mse\n')
        for index in indices:
            file_name = f"{args.prefix}_{index}.trees"
            try:
                ts = tskit.load(file_name)
            except FileNotFoundError:
                print(f"File not found: {file_name}")
                f.write(f"{index}\tNA\tNA\n")
                continue
            f.write(f"{index}\t{count_incompatibility(ts)}\t{diversity_fit_mse(ts, args.m)}\n")

if __name__ == "__main__":
    main()
