#!/usr/bin/env python3
"""
Apply orientation to contigs based on alignment analysis.

Reads orientation TSV and FASTA, outputs oriented FASTA with updated headers.
Contigs that need flipping are reverse complemented.
Headers are updated with full metadata including assembler info and orientation.
"""

import argparse
import gzip


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement)[::-1]


def read_orientation_tsv(tsv_file):
    """Read orientation TSV and return dict of contig -> (chrom, flip)."""
    orientations = {}
    with open(tsv_file) as f:
        header = f.readline()  # Skip header
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 3:
                contig = fields[0]
                chrom = fields[1]
                flip = fields[2] == "yes"
                orientations[contig] = (chrom, flip)
    return orientations


def open_fasta(filepath):
    """Open FASTA file, handling gzip compression."""
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, "rt")
    return open(filepath)


def parse_fasta(fasta_file):
    """Parse FASTA file, yielding (header, sequence) tuples."""
    with open_fasta(fasta_file) as f:
        header = None
        seq_parts = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]  # Remove '>'
                seq_parts = []
            else:
                seq_parts.append(line)

        if header is not None:
            yield header, "".join(seq_parts)


def extract_contig_name(header):
    """Extract the contig name from PanSN-spec header.

    Header format: sample#hap#contig metadata...
    Returns the full PanSN name (sample#hap#contig)
    """
    parts = header.split()
    return parts[0] if parts else header


def main():
    parser = argparse.ArgumentParser(
        description="Apply orientation to contigs and add metadata to headers."
    )
    parser.add_argument("input_fa", help="Input FASTA file (can be gzipped)")
    parser.add_argument("orientation_tsv", help="Orientation TSV from analyze_orientation")
    parser.add_argument("output_fa", help="Output FASTA file (use /dev/stdout for pipe)")
    parser.add_argument("--assembler", default=None, help="Assembler name (e.g., hifiasm)")
    parser.add_argument("--version", default=None, help="Assembler version")
    parser.add_argument("--phasing", default=None, help="Phasing method (trio, hic, none)")
    parser.add_argument("--ultralong", default=None, help="Whether ultralong reads used (yes/no)")

    args = parser.parse_args()

    orientations = read_orientation_tsv(args.orientation_tsv)

    with open(args.output_fa, "w") as out:
        for header, seq in parse_fasta(args.input_fa):
            contig_name = extract_contig_name(header)

            # Build metadata string
            metadata_parts = []
            if args.assembler:
                metadata_parts.append(f"assembler={args.assembler}")
            if args.version:
                metadata_parts.append(f"version={args.version}")
            if args.phasing:
                metadata_parts.append(f"phasing={args.phasing}")
            if args.ultralong:
                metadata_parts.append(f"ultralong={args.ultralong}")

            # Look up orientation info
            if contig_name in orientations:
                chrom, flip = orientations[contig_name]

                # Apply flip if needed
                if flip:
                    seq = reverse_complement(seq)

                # Add orientation info
                flip_str = "yes" if flip else "no"
                metadata_parts.append(f"chrom={chrom}")
                metadata_parts.append(f"flipped={flip_str}")
            else:
                # Contig not in alignments (unmapped), keep as is
                metadata_parts.append("chrom=unplaced")
                metadata_parts.append("flipped=no")

            # Build new header
            metadata_str = " ".join(metadata_parts)
            new_header = f"{contig_name} {metadata_str}" if metadata_str else contig_name

            # Write output
            out.write(f">{new_header}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")


if __name__ == "__main__":
    main()
