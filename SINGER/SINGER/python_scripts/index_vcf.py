import argparse
import gzip
import os
import sys

def index_vcf(input_prefix, segment_length):
    input_file = f"{input_prefix}.vcf.gz" if os.path.exists(f"{input_prefix}.vcf.gz") else f"{input_prefix}.vcf"
    index_file = f"{input_prefix}.index"

    current_segment_start = -1
    byte_offset = 0
    last_variant_segment = -1

    with open(index_file, 'w') as f:
        f.write("")

    with (gzip.open if input_file.endswith(".gz") else open)(input_file, 'rb') as f:
        for line in f:
            if line.startswith(b"#"):
                byte_offset += len(line)
                continue

            pos = int(line.split(b"\t")[1])

            segment_start = (pos // segment_length) * segment_length

            if segment_start != current_segment_start:
                with open(index_file, 'a') as f_idx:
                    f_idx.write(f"{segment_start}\t{byte_offset}\n")

                last_variant_segment = segment_start
                current_segment_start = segment_start

            byte_offset += len(line)

def main():
    parser = argparse.ArgumentParser(description="Index a VCF file by block length.")
    parser.add_argument("vcf_file_prefix", type=str, help="VCF file prefix without .vcf or .vcf.gz extension")
    parser.add_argument("segment_length", type=int, help="Length of segments to index")
    args = parser.parse_args()

    print("Index vcf files")
    print("------------------")
    print(f"VCF file prefix: {args.vcf_file_prefix}")
    print(f"Block length: {args.segment_length}")
    print("------------------")

    index_vcf(args.vcf_file_prefix, args.segment_length)

if __name__ == "__main__":
    main()
