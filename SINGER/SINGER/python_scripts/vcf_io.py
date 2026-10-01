import contextlib
import gzip
import os
import subprocess
import sys


def fail(message):
    sys.exit(f"Error: {message}.")


def records(f, path, chrom):
    seen = False
    try:
        for line in f:
            if not line.strip():
                continue
            if chrom is None or line.startswith("#"):
                yield line
            elif line.startswith(chrom + "\t"):
                seen = True
                yield line
            elif seen:
                return
    except (OSError, EOFError):
        fail(f"cannot read {path}")
    if chrom is not None and not seen:
        fail(f"no record of chromosome {chrom} in {path}")


@contextlib.contextmanager
def open_vcf(path, chrom=None):
    if not os.path.isfile(path):
        fail(f"file {path} not found")
    bcftools = subprocess.Popen(["bcftools", "view", "--no-version", path], stdout=subprocess.PIPE, text=True) if path.endswith(".bcf") else None
    f = bcftools.stdout if bcftools else (gzip.open(path, "rt") if path.endswith(".gz") else open(path))
    try:
        yield records(f, path, chrom)
    finally:
        if bcftools is not None and f.read(1):
            bcftools.kill()
        f.close()
        if bcftools is not None and bcftools.wait() not in (0, -9):
            fail(f"cannot read {path}")


def first_chrom(path):
    with open_vcf(path) as lines:
        for line in lines:
            if not line.startswith("#"):
                return line.split("\t", 1)[0]


def sample_names(path):
    with open_vcf(path) as lines:
        for line in lines:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    fail(f"no #CHROM line in {path}")


def add_individuals_from_vcf(ts, vcf_file):
    names = sample_names(vcf_file)
    ploidy = ts.num_samples // len(names)
    if ploidy not in (1, 2) or ploidy * len(names) != ts.num_samples:
        fail(f"{ts.num_samples} sample nodes do not match the {len(names)} samples of {vcf_file} at ploidy 1 or 2")
    tables = ts.dump_tables()
    node_individual = tables.nodes.individual.copy()
    node_metadata = [b""] * tables.nodes.num_rows
    for ind_id, name in enumerate(names):
        tables.individuals.add_row(metadata=name.encode())
        for hap in range(ploidy):
            nid = ind_id * ploidy + hap
            node_individual[nid] = ind_id
            node_metadata[nid] = f"{name}_{hap}".encode()
    tables.nodes.individual = node_individual
    tables.nodes.packset_metadata(node_metadata)
    return tables.tree_sequence()
