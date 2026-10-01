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
    ts = tables.tree_sequence()

    return ts

def load_file_lists(file_list_path):
    node_files = []
    branch_files = []
    mutation_files = []
    block_coordinates = []

    with open(file_list_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                raise ValueError(f"Invalid line: {line}")
            node_files.append(parts[0])
            branch_files.append(parts[1])
            mutation_files.append(parts[2])
            block_coordinates.append(float(parts[3]))  # or int if desired

    return node_files, branch_files, mutation_files, block_coordinates

def sort_nodes_by_time(ts):
    tables = ts.dump_tables()
    times = tables.nodes.time
    is_sample = (tables.nodes.flags & tskit.NODE_IS_SAMPLE) > 0
    sort_order = np.lexsort((np.where(is_sample, np.arange(ts.num_nodes, dtype=float), times), ~is_sample))
    
    # Remap all references
    node_map = np.full(ts.num_nodes, tskit.NULL, dtype=int)
    for new_id, old_id in enumerate(sort_order):
        node_map[old_id] = new_id
    
    # Reorder nodes
    tables.nodes.set_columns(
        flags=tables.nodes.flags[sort_order],
        time=tables.nodes.time[sort_order],
        population=tables.nodes.population[sort_order],
        individual=tables.nodes.individual[sort_order],
    )
    
    # Remap edges
    edges = tables.edges
    edges.set_columns(
        left=edges.left,
        right=edges.right,
        parent=node_map[edges.parent].astype(np.int32),
        child=node_map[edges.child].astype(np.int32),
    )

    muts = tables.mutations
    muts.set_columns(
        site=muts.site,
        node=node_map[muts.node].astype(np.int32),
        derived_state=muts.derived_state,
        derived_state_offset=muts.derived_state_offset,
        parent=muts.parent,
        time=muts.time,
        metadata=muts.metadata,
        metadata_offset=muts.metadata_offset,
    )

    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    return tables.tree_sequence()

def write_output_ts(ts, output):
    print(f"Save to {output}")
    ts.dump(output)

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Generate tskit format for a long ARG.")
    
    # Add arguments with prefixes
    parser = argparse.ArgumentParser(description="Generate tskit format for a long ARG using file list.")
    parser.add_argument("--file_table", required=True, help="Sub file table")
    parser.add_argument("--output", required=True, help="Output file name")
    parser.add_argument("--vcf", required=False, default=None,
                        help="Input file of the SINGER run (.vcf, .vcf.gz or .bcf). Sample "
                             "names from the header are attached as individuals, and tips "
                             "are named <sample>_0 / <sample>_1.")

    args = parser.parse_args()

    node_files, branch_files, mutation_files, block_coordinates = load_file_lists(args.file_table)
    ts = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates)
    sorted_ts = sort_nodes_by_time(ts)
    if args.vcf is not None:
        sorted_ts = vcf_io.add_individuals_from_vcf(sorted_ts, args.vcf)
    write_output_ts(sorted_ts, args.output)

if __name__  == "__main__":
    main()
