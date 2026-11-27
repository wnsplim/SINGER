import json
import os
import time

import numpy as np
import tskit

EDGE_SCHEMA = tskit.MetadataSchema({
    'codec': 'json',
    'type': 'object',
    'properties': {
        'expected_mutations': {'type': 'number', 'description': 'expected number of mutations on the edge'},
        'informed_bases': {'type': 'number', 'description': 'number of informed bases of the edge'},
    },
})

MUTATION_SCHEMA = tskit.MetadataSchema({
    'codec': 'json',
    'type': 'object',
    'properties': {
        'informed': {'type': 'integer', 'description': '1 for an informed mutation, else 0'},
    },
})

RECORD = np.dtype([('pos', '<f8'), ('source', '<i4'), ('start', '<f8'), ('target', '<i4')])
MUTATION = np.dtype([('pos', '<f8'), ('lower', '<i4'), ('state', 'i1'), ('informed', 'i1')])
CELL = np.dtype([('pos', '<f8'), ('leaf0', '<i4'), ('state0', 'i1'), ('leaf1', '<i4'), ('state1', 'i1')])


def sample_file(prefix, i):
    return f'{prefix}_{i}.bin'


def read_sample(path):
    with open(path, 'rb') as file:
        b = file.read()
    n, f, nodes, records, scaled, draws = np.frombuffer(b, '<i8', 6).tolist()
    o = 48

    def take(dtype, count):
        nonlocal o
        a = np.frombuffer(b, dtype, int(count), o)
        o += a.nbytes
        return a

    s = {'samples': n, 'length': float(take('<f8', 1)[0]), 'parents': take('<i4', f), 'times': take('<f8', nodes)}
    s['rescaled'] = take('<f8', nodes) if scaled else None
    s['records'] = take(RECORD, records)
    s['draws'] = [take(MUTATION, take('<i8', 1)[0]) for _ in range(draws)]
    s['counts'] = take('<f8', 2 * take('<i8', 1)[0]).reshape(-1, 2)
    cells = take('<i8', 1)[0]
    s['cells'] = take(CELL, cells) if cells else None
    return s


def trees(sample, ended=None):
    par = dict(enumerate(sample['parents'].tolist()))
    kids = {}
    for c, p in par.items():
        kids.setdefault(p, []).append(c)
    left = dict.fromkeys(par, 0.0)
    k = len(par)
    x0 = 0.0

    def end(c, x):
        if ended is not None:
            ended.append((left[c], x, par[c], c))

    def link(c, p, x):
        par[c] = p
        kids.setdefault(p, []).append(c)
        left[c] = x

    for x, s, _, g in sample['records'].tolist():
        yield x0, x, par
        d = par[s]
        sib = kids[d][0] if kids[d][1] == s else kids[d][1]
        pd = par[d]
        end(s, x)
        end(sib, x)
        end(d, x)
        kids[pd].remove(d)
        del par[d], kids[d]
        if g == sib:
            link(sib, k, x)
            link(s, k, x)
            link(k, pd, x)
        else:
            u = par[g]
            end(g, x)
            kids[u].remove(g)
            link(sib, pd, x)
            link(s, k, x)
            link(g, k, x)
            link(k, u, x)
        k += 1
        x0 = x
    yield x0, sample['length'], par
    if ended is not None:
        ended.extend((left[c], sample['length'], p, c) for c, p in par.items())


def trees_name(output, i, draw):
    parts = [output]
    if draw is not None:
        parts.append(str(draw))
    parts.append(str(i))
    return '_'.join(parts) + '.trees'


def add_mutations(tables, mutations, offset, node_map):
    last = None
    for pos, node, state, informed in zip(mutations['pos'].tolist(), node_map(mutations['lower'].astype(np.int64)).tolist(),
                                          mutations['state'].tolist(), mutations['informed'].tolist()):
        if pos != last:
            site = tables.sites.add_row(position=pos + offset, ancestral_state='0')
            last = pos
        tables.mutations.add_row(site=site, node=node, derived_state=str(state), metadata={'informed': informed}, time=tskit.UNKNOWN_TIME)


def finish_mutations(tables):
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    derived = np.array(tskit.unpack_strings(tables.mutations.derived_state, tables.mutations.derived_state_offset))
    ancestral = np.full(tables.sites.num_rows, '0')
    top = (tables.mutations.parent == -1) & (derived == '0')
    ancestral[tables.mutations.site[top]] = '1'
    tables.sites.packset_ancestral_state(ancestral.tolist())
    tables.compute_mutation_times()


def read_long_ARG(sample_files, block_coordinates, sequence_length=0, draw=None):
    tables = tskit.TableCollection(sequence_length=sequence_length)
    tables.edges.metadata_schema = EDGE_SCHEMA
    tables.mutations.metadata_schema = MUTATION_SCHEMA
    for path, offset in zip(sample_files, block_coordinates):
        s = read_sample(path)
        n = s['samples']
        times = s['times'] if s['rescaled'] is None else s['rescaled']
        if tables.nodes.num_rows == 0:
            tables.nodes.append_columns(flags=np.full(n, tskit.NODE_IS_SAMPLE, np.uint32), time=times[:n])
        base = tables.nodes.num_rows - n
        tables.nodes.append_columns(flags=np.zeros(len(times) - n, np.uint32), time=times[n:])
        node_map = lambda a, n=n, base=base: np.where(a >= n, a + base, a).astype(np.int32)
        ended = []
        for _ in trees(s, ended):
            pass
        e = np.array(ended)
        e = e[e[:, 2] >= 0]
        e = e[np.lexsort((e[:, 0], e[:, 3]))]
        metadata, metadata_offset = tskit.pack_bytes([EDGE_SCHEMA.validate_and_encode_row({'expected_mutations': x, 'informed_bases': y})
                                                      for x, y in s['counts'].tolist()])
        tables.edges.append_columns(left=e[:, 0] + offset, right=e[:, 1] + offset, parent=node_map(e[:, 2].astype(np.int64)),
                                    child=node_map(e[:, 3].astype(np.int64)), metadata=metadata, metadata_offset=metadata_offset)
        tables.sequence_length = max(tables.sequence_length, s['length'] + offset)
        if draw is not None:
            add_mutations(tables, s['draws'][draw], offset, node_map)
    finish_mutations(tables)
    return tables


def read_ckp(prefix):
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
