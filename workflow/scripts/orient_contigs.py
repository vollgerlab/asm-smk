#!/usr/bin/env python3
"""
Analyze PAF alignments to determine contig orientation relative to reference.

For each contig, determines:
1. Primary chromosome (chromosome with most aligned bases)
2. Whether majority of bases align to forward or reverse strand
3. Outputs orientation info for downstream processing
"""

import sys
from collections import defaultdict


def parse_paf(paf_file):
    """Parse PAF file and collect alignment statistics per contig."""
    # contig -> chrom -> {'+': bases, '-': bases}
    alignments = defaultdict(lambda: defaultdict(lambda: {"+": 0, "-": 0}))

    with open(paf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 12:
                continue

            query_name = fields[0]
            query_len = int(fields[1])
            strand = fields[4]
            target_name = fields[5]
            # Use alignment block length (column 11) as the aligned bases
            aln_block_len = int(fields[10])

            alignments[query_name][target_name][strand] += aln_block_len

    return alignments


def determine_orientation(alignments):
    """
    For each contig, determine primary chromosome and orientation.

    Returns dict: contig -> (primary_chrom, should_flip, total_bases)
    """
    orientations = {}

    for contig, chrom_data in alignments.items():
        # Find chromosome with most total aligned bases
        chrom_totals = {}
        for chrom, strand_data in chrom_data.items():
            chrom_totals[chrom] = strand_data["+"] + strand_data["-"]

        if not chrom_totals:
            continue

        primary_chrom = max(chrom_totals, key=chrom_totals.get)
        total_bases = chrom_totals[primary_chrom]

        # Determine strand orientation for primary chromosome
        fwd_bases = chrom_data[primary_chrom]["+"]
        rev_bases = chrom_data[primary_chrom]["-"]

        # Flip if majority of bases are on reverse strand
        should_flip = rev_bases > fwd_bases

        orientations[contig] = (primary_chrom, should_flip, total_bases, fwd_bases, rev_bases)

    return orientations


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.paf> <output.tsv>", file=sys.stderr)
        sys.exit(1)

    paf_file = sys.argv[1]
    output_file = sys.argv[2]

    alignments = parse_paf(paf_file)
    orientations = determine_orientation(alignments)

    with open(output_file, "w") as f:
        # Header
        f.write("contig\tprimary_chrom\tflip\ttotal_bases\tfwd_bases\trev_bases\n")

        for contig, (chrom, flip, total, fwd, rev) in sorted(orientations.items()):
            flip_str = "yes" if flip else "no"
            f.write(f"{contig}\t{chrom}\t{flip_str}\t{total}\t{fwd}\t{rev}\n")


if __name__ == "__main__":
    main()
