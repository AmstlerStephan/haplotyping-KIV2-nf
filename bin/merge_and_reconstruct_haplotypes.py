"""
merge_and_reconstruct_haplotypes.py

Fused process: merge haplotypes AND reconstruct full-length sequences in one pass.
Avoids intermediate I/O and redundant file reads.

Outputs:
  - merged_haplotypes.fasta (fingerprints)
  - reconstructed_haplotypes.fasta (full-length sequences)
  - merged_haplotype_stats.tsv
  - merged_haplotype_log.tsv
"""

import argparse
import logging
from operator import add
import os
import re
import sys
from turtle import pos
import uuid

import pysam
import edlib


def parse_args(argv):
    usage = "Merge haplotypes and reconstruct full-length sequences in one pass"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-l", "--log",
        dest="log",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL",
                 "debug", "info", "warning", "error", "critical"],
        default="INFO",
        help="Logging verbosity",
    )
    parser.add_argument(
        "--fastx_file",
        dest="FASTX_FILE",
        type=str,
        required=True,
        help="Input fastx file (filtered haplotypes)",
    )
    parser.add_argument(
        "--positions",
        dest="POSITIONS",
        type=str,
        required=True,
        help="TSV with polymorphic positions",
    )
    parser.add_argument(
        "--reference",
        dest="REFERENCE",
        type=str,
        required=True,
        help="FASTA file containing the reference sequence",
    )
    parser.add_argument(
        "--reference_start",
        dest="REFERENCE_START",
        type=int,
        default=1,
        help="1-based coordinate offset for reference (default: 1)",
    )
    parser.add_argument(
        "--output_format",
        dest="OUTPUT_FORMAT",
        type=str,
        default="fasta",
        help="Output format of the haplotypes",
    )
    parser.add_argument(
        "--variant_cutoff",
        dest="VARIANT_CUTOFF",
        type=float,
        default=0.0085,
        help="Cutoff for variant cluster merging",
    )
    parser.add_argument(
        "--max_edit_distance",
        dest="MAX_EDIT_DISTANCE",
        type=int,
        default=2,
        help="Maximum edit distance for merging",
    )
    parser.add_argument(
        "-o", "--output",
        dest="OUTPUT",
        default="./",
        help="Output folder",
    )
    return parser.parse_args(argv)


# ============================================================================
# Merging logic (from merge_haplotypes.py)
# ============================================================================

def get_unique_sequences(fasta_file):
    """Extract unique sequences and track quality and read membership."""
    unique_sequences = dict()
    with pysam.FastxFile(fasta_file) as reads:
        for read in reads:

            sequence = read.sequence
            high_qual = all(base.isupper() for base in read.sequence)
            if sequence in unique_sequences:
                unique_sequences[sequence]["reads"][read.name] = read.sequence
            else:
                unique_sequences[sequence] = dict()
                unique_sequences[sequence]["reads"] = dict()
                unique_sequences[sequence]["reads"][read.name] = read.sequence
            unique_sequences[sequence]["high_qual"] = high_qual
    return unique_sequences


def get_number_of_sequences(unique_sequences):
    """Count total number of reads across all sequences."""
    return sum(len(info["reads"]) for info in unique_sequences.values())


def compute_distance_matrix(unique_sequences, max_edit_distance):
    """Precompute all pairwise edit distances."""
    distance_matrix = {}
    sequences = list(unique_sequences.keys())
    for i, seq1 in enumerate(sequences):
        distance_matrix[seq1] = {}
        for seq2 in sequences[i + 1:]:
            result = edlib.align(
                seq1, seq2, mode="NW", task="path", k=max_edit_distance, additionalEqualities=[("A", "a"), ("C", "c"), ("G", "g"), ("T", "t"), ("-", "A"), ("-", "C"), ("-", "G"), ("-", "T"), ("-", "a"), ("-", "c"), ("-", "g"), ("-", "t")]
            )
            distance_matrix[seq1][seq2] = result
    return distance_matrix


def update_distance_matrix(merged_sequences, distance_matrix, close_sequences, max_edit_distance):
    """Remove merged sequences from matrix."""
    for sequence, queries in close_sequences.items():
        for query in queries:
            if query in distance_matrix:
                del distance_matrix[query]
            for key in distance_matrix:
                if query in distance_matrix[key]:
                    del distance_matrix[key][query]
    return distance_matrix


def find_merges_from_matrix(merged_sequences, distance_matrix, variant_cutoff, max_dist, stats_file_path):
    """Find merge candidates using precomputed distance matrix."""
    n_total_sequences = get_number_of_sequences(merged_sequences)
    cluster_cutoff = round(variant_cutoff * n_total_sequences)
    close_sequences = dict()

    for sequence, info in merged_sequences.items():
        n_sequences = len(info["reads"])
        is_bigger_than_cluster_cutoff = n_sequences > cluster_cutoff
        is_high_qual = info["high_qual"]

        if not (is_bigger_than_cluster_cutoff):
            continue

        for query_sequence, query_info in merged_sequences.items():
            if query_sequence == sequence:
                continue

            n_queries = len(query_info["reads"])
            is_smaller_than_cluster_cutoff = n_queries <= cluster_cutoff

            if not is_smaller_than_cluster_cutoff:
                continue

            # Fetch precomputed distance
            if sequence in distance_matrix and query_sequence in distance_matrix[sequence]:
                result = distance_matrix[sequence][query_sequence]
            elif query_sequence in distance_matrix and sequence in distance_matrix[query_sequence]:
                result = distance_matrix[query_sequence][sequence]
            else:
                result = edlib.align(
                    sequence, query_sequence, mode="NW", task="path", k=max_dist, additionalEqualities=[("A", "a"), ("C", "c"), ("G", "g"), ("T", "t"), ("-", "A"), ("-", "C"), ("-", "G"), ("-", "T"), ("-", "a"), ("-", "c"), ("-", "g"), ("-", "t")]
                )

            if result.get("editDistance", float('inf')) == max_dist:
                if sequence in close_sequences:
                    close_sequences[sequence].append(query_sequence)
                else:
                    close_sequences[sequence] = [query_sequence]
                write_merge_log(sequence, query_sequence, n_queries, result, max_dist,
                               cluster_cutoff, variant_cutoff, n_total_sequences, stats_file_path)

    return close_sequences


def merge_sequences(unique_sequences, close_sequences):
    """Merge close sequences into representative sequences."""
    for sequence, queries in close_sequences.items():
        for query in queries:
            unique_sequences[sequence]["reads"].update(
                unique_sequences[query]["reads"]
            )

    for sequence, queries in close_sequences.items():
        for query in queries:
            if query in unique_sequences:
                unique_sequences.pop(query)

    return unique_sequences


def write_merge_log(sequence, query_sequence, n_queries, result, max_dist,
                   cluster_cutoff, variant_cutoff, n_unique_sequences, stats_file_path):
    """Log merge decision details."""
    import re
    edist = result.get("editDistance", -1)
    cigar = result.get("cigar", "")
    with open(stats_file_path, "a+") as stats_file:
        for difference in re.findall(r'\d*=..', cigar):
            change = difference.split("=")[1]
            n_bases = int(re.findall(r"\d*", change)[0])
            pos = int(difference.split("=")[0])
            base = sequence[pos:pos+n_bases]
            query_base = query_sequence[pos:pos+n_bases]
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                sequence, query_sequence, n_queries, edist, pos, base,
                query_base, change, max_dist, cluster_cutoff,
                variant_cutoff, n_unique_sequences), file=stats_file)


# ============================================================================
# Reconstruction logic (from reconstruct_haplotypes.py)
# ============================================================================

def read_reference(reference_fasta):
    """Load reference sequence from FASTA."""
    with pysam.FastxFile(reference_fasta) as fh:
        for record in fh:
            logging.info("Using reference sequence '%s' (length %d bp)",
                        record.name, len(record.sequence))
            return record.sequence
    raise ValueError("Reference FASTA '{}' contains no sequences.".format(reference_fasta))


def read_positions(positions_tsv):
    """Load sorted list of polymorphic positions."""
    positions = []
    with open(positions_tsv) as fh:
        header = fh.readline().strip()
        if header != "position":
            raise ValueError("Expected header 'position' in '{}'".format(positions_tsv))
        for line in fh:
            line = line.strip()
            if line:
                positions.append(int(line))
    return sorted(positions)


def apply_haplotype(reference, positions, fingerprint, reference_start):
    """Reconstruct full-length sequence by applying polymorphic positions."""
    offset = reference_start - 1
    ref_list = list(reference)
    shift = 0

    for index in range(len(positions)):
        fingerprint_index = index + shift
        allele = fingerprint[fingerprint_index]
        ref_pos = positions[index]
        ref_idx = ref_pos - 1 - offset

        # find indels by checking if the next base in the fingerprint is lower case, if so, it is part of an indel and all lower bases should be added at this position
        if fingerprint_index + 1 < len(fingerprint) and fingerprint[fingerprint_index + 1].islower():
            while fingerprint_index + 1 < len(fingerprint) and fingerprint[fingerprint_index + 1].islower():
                fingerprint_index += 1
                shift += 1
                allele += fingerprint[fingerprint_index]
            ref_list[ref_idx] = allele
        if allele == "-":
            # deletion_indices.add(ref_idx)
            ref_list[ref_idx] = ""
        else:
            ref_list[ref_idx] = allele

    reconstructed = "".join(ref_list)
    return reconstructed



def write_haplotype_stats(merged_sequences, output):
    """Write stats for merged haplotypes."""
    haplotype_stats_file = os.path.join(output, "merged_haplotype_stats.tsv")
    with open(haplotype_stats_file, "w") as out_f:
        print("haplotype\thaplotype_occurences\thigh_qual\thaplotype_length",
              file=out_f)
        for sequence, info in merged_sequences.items():
            n_sequences = len(info["reads"])
            haplotype_length = len(sequence)
            high_qual = info["high_qual"]
            print("{}\t{}\t{}\t{}".format(sequence, n_sequences, high_qual,
                                        haplotype_length), file=out_f)


# ============================================================================
# Main fused process
# ============================================================================

def merge_and_reconstruct(args):
    # Full-length reconstruction requires a valid reference sequence.
    if os.path.basename(args.REFERENCE) == "NO_FILE.txt":
        raise ValueError(
            "Full-length reconstruction requires a region reference FASTA. "
            "Please set params.region_references for this region."
        )

    # Load data once
    reference = read_reference(args.REFERENCE)
    positions = read_positions(args.POSITIONS)
    logging.info("Loaded %d polymorphic positions", len(positions))

    # Extract and merge haplotypes
    unique_sequences = get_unique_sequences(args.FASTX_FILE)

    stats_file_path = os.path.join(args.OUTPUT, "merged_haplotype_log.tsv")
    with open(stats_file_path, "w") as f:
        print("sequence\tquery_sequence\tquery_size\tedist\tposition\tbase\t"
              "query_base\tchange\tmax_dist\tcluster_cutoff\tvariant_cutoff\t"
              "n_unique_sequences", file=f)

    distance_matrix = compute_distance_matrix(unique_sequences, args.MAX_EDIT_DISTANCE)
    merged_sequences = unique_sequences.copy()

    for max_dist in range(1, args.MAX_EDIT_DISTANCE + 1):
        close_sequences = find_merges_from_matrix(
            merged_sequences, distance_matrix, args.VARIANT_CUTOFF, max_dist, stats_file_path)
        if not close_sequences:
            break
        merged_sequences = merge_sequences(merged_sequences, close_sequences)
        distance_matrix = update_distance_matrix(
            merged_sequences, distance_matrix, close_sequences, args.MAX_EDIT_DISTANCE)

    # Write merged haplotypes fingerprints and stats
    merged_fasta = os.path.join(args.OUTPUT, "merged_haplotypes.fasta")
    with open(merged_fasta, "w") as f:
        for i, sequence in enumerate(merged_sequences):
            n_reads = len(merged_sequences[sequence]["reads"])
            high_qual = merged_sequences[sequence]["high_qual"]
            unique_id = uuid.uuid4()
            header = "{},size={},high_qual={},uuid={}".format(
                i, n_reads, high_qual, unique_id)
            f.write(">{}\n{}\n".format(header, sequence))

    write_haplotype_stats(merged_sequences, args.OUTPUT)

    # Reconstruct full-length sequences in same pass
    reconstructed_fasta = os.path.join(args.OUTPUT, "reconstructed_haplotypes.fasta")
    n_reconstructed = 0

    with open(reconstructed_fasta, "w") as out_f:
        for i, (fingerprint, info) in enumerate(merged_sequences.items()):
            try:
                reconstructed = apply_haplotype(
                    reference, positions, fingerprint, args.REFERENCE_START)
            except (ValueError, IndexError) as exc:
                logging.warning("Skipping haplotype %d: %s", i, exc)
                continue

            header = "cluster_{},size={},n_reads_reconstructed={},high_qual={}".format(
                i, len(info["reads"]), len(info["reads"]), info["high_qual"])
            out_f.write(">{}\n{}\n".format(header, reconstructed))
            n_reconstructed += 1

    logging.info(
        "Wrote %d merged fingerprints and %d reconstructed sequences",
        len(merged_sequences), n_reconstructed)


def main(argv=sys.argv[1:]):
    args = parse_args(argv)
    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log)
    logging.basicConfig(level=numeric_level, format="%(message)s")
    merge_and_reconstruct(args)


if __name__ == "__main__":
    main()
