import argparse
import gzip
import os
import sys
import vcf_io

def index_vcf(input_file, segment_length):
    if not os.path.isfile(input_file):
        vcf_io.fail(f"file {input_file} not found")
    if input_file.endswith(".bcf"):
        vcf_io.fail(f"{input_file} is BCF, which SINGER reads without an index")
    index_file = f"{input_file}.index"

    current = (None, -1)
    byte_offset = 0

    with open(index_file, 'w') as f:
        f.write("")

    with (gzip.open if input_file.endswith(".gz") else open)(input_file, 'rb') as f:
        for line in f:
            if line.startswith(b"#") or not line.strip():
                byte_offset += len(line)
                continue

            fields = line.split(b"\t")
            chrom = fields[0].decode()
            pos = int(fields[1])

            segment_start = (pos // segment_length) * segment_length

            if (chrom, segment_start) != current:
                with open(index_file, 'a') as f_idx:
                    f_idx.write(f"{chrom}\t{segment_start}\t{byte_offset}\n")
                current = (chrom, segment_start)

            byte_offset += len(line)

def main():
    parser = argparse.ArgumentParser(description="Index a VCF file by block length.")
    parser.add_argument("vcf_file", type=str, help="VCF file (.vcf or .vcf.gz); the index is written to VCF_FILE.index")
    parser.add_argument("segment_length", type=int, help="Length of segments to index")
    args = parser.parse_args()

    print("Index vcf files")
    print("------------------")
    print(f"VCF file: {args.vcf_file}")
    print(f"Block length: {args.segment_length}")
    print("------------------")

    index_vcf(args.vcf_file, args.segment_length)

if __name__ == "__main__":
    main()
