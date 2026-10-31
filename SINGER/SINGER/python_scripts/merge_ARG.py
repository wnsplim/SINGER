import argparse
import numpy as np
import tskit
import vcf_io
import arg_io
from convert_long_ARG import read_long_ARG

def load_file_lists(file_list_path):
    """Rows: nodes branches muts offset, or nodes branches muts counts offset."""
    node_files = []
    branch_files = []
    mutation_files = []
    counts_files = []
    block_coordinates = []

    with open(file_list_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                raise ValueError(f"Invalid line: {line}")
            node_files.append(parts[0])
            branch_files.append(parts[1])
            mutation_files.append(parts[2])
            counts_files.append(parts[3] if len(parts) >= 5 else None)
            block_coordinates.append(float(parts[-1]))

    return node_files, branch_files, mutation_files, counts_files, block_coordinates

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
        metadata=edges.metadata,
        metadata_offset=edges.metadata_offset,
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
    parser = argparse.ArgumentParser(description="Generate tskit format for a long ARG using file list.")
    parser.add_argument("--file_table", required=True, help="Sub file table: one row per block, 'nodes branches muts offset' or 'nodes branches muts counts offset'")
    parser.add_argument("--output", required=True, help="Output file name")
    parser.add_argument("--vcf", required=False, default=None,
                        help="Input file of the SINGER run (.vcf, .vcf.gz or .bcf). Sample "
                             "names from the header are attached as individuals, and tips "
                             "are named <sample>_0 / <sample>_1.")

    args = parser.parse_args()

    node_files, branch_files, mutation_files, counts_files, block_coordinates = load_file_lists(args.file_table)
    tables = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates, counts_files)
    arg_io.add_provenance(tables, 'merge_ARG.py', [f.rsplit('_nodes_', 1)[0] for f in node_files], None, None, None)
    sorted_ts = sort_nodes_by_time(tables.tree_sequence())
    if args.vcf is not None:
        sorted_ts = vcf_io.add_individuals_from_vcf(sorted_ts, args.vcf)
    write_output_ts(sorted_ts, args.output)

if __name__  == "__main__":
    main()
