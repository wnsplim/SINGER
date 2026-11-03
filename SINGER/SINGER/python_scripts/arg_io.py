import glob
import json
import os
import re
import time

import numpy as np
import tskit

EDGE_SCHEMA = tskit.MetadataSchema({
    'codec': 'json',
    'type': 'object',
    'properties': {
        'expected_mutations': {'type': 'number', 'description': 'expected number of mutations on the edge over the sites at which the data can see it'},
        'informed_bases': {'type': 'number', 'description': 'number of bases of the edge at which at least one sample below it and one outside it have data'},
    },
})

MUTATION_SCHEMA = tskit.MetadataSchema({
    'codec': 'json',
    'type': 'object',
    'properties': {
        'informed': {'type': 'integer', 'description': '1 when the data can see the branch of the mutation at this site, 0 when the draw placed it on a branch they cannot'},
    },
})


def muts_files(prefix, i):
    """The mapping draws of sample i as (draw index or None, path): one file <prefix>_muts_<i>.txt, or the files
    <prefix>_muts_<r>_<i>.txt of a run with -mapping_rep above 1, sorted by r."""
    single = f'{prefix}_muts_{i}.txt'
    if os.path.exists(single):
        return [(None, single)]
    found = []
    for f in glob.glob(f'{prefix}_muts_*_{i}.txt'):
        m = re.fullmatch(re.escape(prefix) + r'_muts_(\d+)_' + str(i) + r'\.txt', f)
        if m:
            found.append((int(m.group(1)), f))
    return sorted(found)


def nodes_file(prefix, i):
    """The node times of sample i: the rescaled ones when the run wrote them, else the sampled ones; (path, rescaled)."""
    rescaled = f'{prefix}_rescaled_nodes_{i}.txt'
    if os.path.exists(rescaled):
        return rescaled, True
    return f'{prefix}_nodes_{i}.txt', False


def trees_name(output, i, draw):
    """<output>[_<draw>]_<i>.trees"""
    parts = [output]
    if draw is not None:
        parts.append(str(draw))
    parts.append(str(i))
    return '_'.join(parts) + '.trees'


def edge_metadata(counts_file, keep):
    """The packed edge metadata column for the branches rows selected by the boolean mask keep, from the counts file
    (one row per branches row: expected count, informed bases), or None when the file is absent."""
    if counts_file is None or not os.path.exists(counts_file):
        return None
    counts = np.loadtxt(counts_file, ndmin=2)
    if counts.shape[0] != len(keep):
        raise ValueError(f'{counts_file} has {counts.shape[0]} rows, the branches file {len(keep)}')
    rows = [EDGE_SCHEMA.validate_and_encode_row({'expected_mutations': float(e), 'informed_bases': float(b)}) for e, b in counts[keep]]
    return tskit.pack_bytes(rows)


def add_edges(tables, edge_span, parent, child, offset, counts_file, keep):
    """Appends the branches rows (already filtered by keep) as edges with the counts as metadata when available."""
    tables.edges.metadata_schema = EDGE_SCHEMA
    packed = edge_metadata(counts_file, keep)
    if packed is None:
        metadata, metadata_offset = tskit.pack_bytes([b''] * len(parent))
    else:
        metadata, metadata_offset = packed
    tables.edges.append_columns(left=edge_span[:, 0] + offset, right=edge_span[:, 1] + offset, parent=parent, child=child,
                                metadata=metadata, metadata_offset=metadata_offset)


BASES = 'ACGT'


def add_mutations(tables, mutation_file, offset, node_map, bases=False):
    """The muts rows as sites and mutations; the fifth column, when present, as the mutation's informed flag. Without
    a file, no site. With bases (a run under a four-state model), the row whose upper node is the root gives the
    site's ancestral base and the states of the other rows are bases."""
    tables.mutations.metadata_schema = MUTATION_SCHEMA
    if mutation_file is None:
        return
    mutations = np.loadtxt(mutation_file, ndmin=2)
    by_site = {}
    for row in mutations:
        by_site.setdefault(row[0], []).append(row)
    for pos in sorted(by_site):
        rows = by_site[pos]
        ancestral = '0'
        if bases:
            ancestral = next(BASES[int(row[3])] for row in rows if int(row[2]) == -1)
            rows = [row for row in rows if int(row[2]) != -1]
        tables.sites.add_row(position=pos + offset, ancestral_state=ancestral)
        site = tables.sites.num_rows - 1
        for row in rows:
            meta = {'informed': int(row[4])} if len(row) >= 5 else {}
            state = BASES[int(row[3])] if bases else str(int(row[3]))
            tables.mutations.add_row(site=site, node=node_map(int(row[1])), derived_state=state, metadata=meta)


def uses_bases(prefix):
    """Whether the run at prefix used a four-state model (its .ckp names -model other than CFN)."""
    return read_ckp(prefix).get('-model', 'CFN') != 'CFN'


def read_ckp(prefix):
    """The flags of the run from <prefix>.ckp as a dict, or {} when absent."""
    flags = {}
    path = f'{prefix}.ckp'
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                w = line.split(None, 1)
                if w:
                    flags[w[0]] = w[1].strip() if len(w) > 1 else ''
    return flags


def add_provenance(tables, converter, prefixes, sample, draw, rescaled):
    record = {
        'schema_version': '1.0.0',
        'software': {'name': 'SINGER', 'converter': converter},
        'parameters': {'runs': [{'prefix': p, 'flags': read_ckp(p)} for p in prefixes], 'sample': sample, 'draw': draw, 'rescaled': rescaled},
        'environment': {'tskit': tskit.__version__},
    }
    tables.provenances.add_row(timestamp=time.strftime('%Y-%m-%dT%H:%M:%S'), record=json.dumps(record))
