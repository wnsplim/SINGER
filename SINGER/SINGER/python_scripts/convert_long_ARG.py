import argparse
import numpy as np
import tskit
import vcf_io

def read_long_ARG(node_files, branch_files, mutation_files, block_coordinates):
    if len(node_files) != len(branch_files):
        raise ValueError("Lengths of node_files and branch_files must be the same.")
    
    if len(node_files) != len(block_coordinates):
        raise ValueError("Lengths of node_files and coordinates must be the same.")    
    
    tables = tskit.TableCollection(sequence_length=0)
    node_table = tables.nodes
    branch_table = tables.edges
    
    node_num = 0
    sample_num = 0

    for node_file_index, (node_file, branch_file, mutation_file) in enumerate(zip(node_files, branch_files, mutation_files)):
        print(f"Processing segment {node_file_index}")
        node_time = np.atleast_1d(np.loadtxt(node_file))
        edge_span = np.loadtxt(branch_file, ndmin=2)
        edge_span = edge_span[edge_span[:, 2] >= 0, :]
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
        
        branch_table.append_columns(
            left=edge_span[:, 0] + block_coordinates[node_file_index],
            right=edge_span[:, 1] + block_coordinates[node_file_index],
            parent=parent_indices,
            child=child_indices
        )
        mutations = np.loadtxt(mutation_file, ndmin=2)
        site_ids = {}
        for i in range(mutations.shape[0]):
            mut_pos = mutations[i, 0]
            if mut_pos >= length:
                continue
            if mut_pos not in site_ids:
                tables.sites.add_row(position=mut_pos + block_coordinates[node_file_index], ancestral_state='0')
                site_ids[mut_pos] = tables.sites.num_rows - 1
            mut_node = int(mutations[i, 1])
            if mut_node >= sample_num:
                mut_node += node_num
            tables.mutations.add_row(site=site_ids[mut_pos], node=mut_node, derived_state=str(int(mutations[i, 3])))
    
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    ts = tables.tree_sequence()

    return ts

def generate_file_lists(vcf_prefix, output_prefix, MCMC_iteration):
    # Read block_coordinates from the index file
    with open(f"{vcf_prefix}.index", 'r') as f:
        rows = [line.split() for line in f if line.strip()]
    block_coordinates = [int(float(w[-2])) for w in rows if w[0] == rows[0][0]]

    # Generate node_files and branch_files
    node_files = []
    branch_files = []
    mutation_files = []

    for i in range(len(block_coordinates)):
        node_files.append(f"{output_prefix}_{i}_{i+1}_nodes_{MCMC_iteration}.txt")
        branch_files.append(f"{output_prefix}_{i}_{i+1}_branches_{MCMC_iteration}.txt")
        mutation_files.append(f"{output_prefix}_{i}_{i+1}_muts_{MCMC_iteration}.txt")
    return node_files, branch_files, mutation_files, block_coordinates

def write_output_ts(ts, output_prefix, MCMC_iteration):
    output_ts_filename = f"{output_prefix}_{MCMC_iteration}.trees"
    print(f"Save to {output_ts_filename}")
    ts.dump(output_ts_filename)

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Generate tskit format for a long ARG.")
    
    # Add arguments with prefixes
    parser.add_argument("-vcf", required=True, help="Input file of the SINGER run (.vcf, .vcf.gz or .bcf)")
    parser.add_argument("-output", required=True, help="Output files prefix")
    parser.add_argument("-iteration", type=int, required=True, help="MCMC iteration for generating filenames")
        
    args = parser.parse_args()

    # Generate file lists
    node_files, branch_files, mutation_files, block_coordinates = generate_file_lists(args.vcf, args.output, args.iteration)
    ts = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates)
    ts = vcf_io.add_individuals_from_vcf(ts, args.vcf)
    write_output_ts(ts, args.output, args.iteration)

if __name__  == "__main__":
    main()
