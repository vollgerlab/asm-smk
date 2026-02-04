#!/usr/bin/env python3
"""
Merge two haplotype-specific VCFs into a single phased diploid VCF.

Takes VCF files from hap1 and hap2 variant calls and combines them into
a single VCF with phased genotypes (e.g., 1|0 for hap1-only, 0|1 for hap2-only,
1|1 for both haplotypes).
"""

import argparse
import gzip
import sys
from collections import defaultdict


def open_vcf(path):
    """Open VCF file, handling gzip compression."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_vcf_records(vcf_path):
    """Parse VCF and return dict of (chrom, pos, ref, alt) -> record."""
    records = {}
    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _, ref, alt = fields[:5]
            key = (chrom, int(pos), ref, alt)
            records[key] = fields
    return records


def get_vcf_header(vcf_path, sample_name):
    """Extract and modify VCF header for phased output."""
    header_lines = []
    with open_vcf(vcf_path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("##"):
                header_lines.append(line.rstrip())
            elif line.startswith("#CHROM"):
                # Replace sample column
                fields = line.strip().split("\t")
                fields = fields[:9] + [sample_name]
                header_lines.append("\t".join(fields))

    # Add phasing format if not present
    has_gt_format = any("##FORMAT=<ID=GT" in h for h in header_lines)
    if not has_gt_format:
        idx = next(
            (i for i, h in enumerate(header_lines) if h.startswith("##FORMAT")),
            len(header_lines) - 1,
        )
        header_lines.insert(
            idx,
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        )

    return header_lines


def merge_haplotype_vcfs(hap1_path, hap2_path, sample_name, ref_path, output):
    """Merge hap1 and hap2 VCFs into phased diploid VCF."""
    # Parse records from both haplotypes
    hap1_records = parse_vcf_records(hap1_path)
    hap2_records = parse_vcf_records(hap2_path)

    # Get all unique variant positions
    all_keys = set(hap1_records.keys()) | set(hap2_records.keys())

    # Sort by chromosome and position
    sorted_keys = sorted(all_keys, key=lambda x: (x[0], x[1]))

    # Write header
    header = get_vcf_header(hap1_path, sample_name)
    for line in header:
        print(line, file=output)

    # Write records with phased genotypes
    for key in sorted_keys:
        chrom, pos, ref, alt = key

        in_hap1 = key in hap1_records
        in_hap2 = key in hap2_records

        # Determine phased genotype
        if in_hap1 and in_hap2:
            gt = "1|1"
            record = hap1_records[key]
        elif in_hap1:
            gt = "1|0"
            record = hap1_records[key]
        else:
            gt = "0|1"
            record = hap2_records[key]

        # Build output record
        # Use first 8 columns from source, add FORMAT and sample
        out_fields = record[:8]

        # Ensure we have INFO field
        if len(out_fields) < 8:
            out_fields.extend(["."] * (8 - len(out_fields)))

        # Add FORMAT and genotype
        out_fields.append("GT")
        out_fields.append(gt)

        print("\t".join(out_fields), file=output)


def main():
    parser = argparse.ArgumentParser(
        description="Merge haplotype-specific VCFs into phased diploid VCF"
    )
    parser.add_argument("--hap1", required=True, help="Haplotype 1 VCF file")
    parser.add_argument("--hap2", required=True, help="Haplotype 2 VCF file")
    parser.add_argument("--sample", required=True, help="Sample name for output")
    parser.add_argument("--ref", required=True, help="Reference FASTA (for header)")
    parser.add_argument(
        "--output", default="/dev/stdout", help="Output VCF (default: stdout)"
    )

    args = parser.parse_args()

    if args.output == "/dev/stdout":
        output = sys.stdout
    else:
        output = open(args.output, "w")

    try:
        merge_haplotype_vcfs(
            args.hap1, args.hap2, args.sample, args.ref, output
        )
    finally:
        if output != sys.stdout:
            output.close()


if __name__ == "__main__":
    main()
