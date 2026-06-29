"""
reconstruct_haplotypes.py

Reconstruct full-length haplotype sequences by substituting the polymorphic
positions of a reference sequence with the alleles present in each merged
haplotype fingerprint.

Algorithm
---------
1. Load the per-region reference sequence (FASTA, first record used).
2. Read the sorted polymorphic positions (haplotypes_filtered_positions.tsv).
3. For every merged haplotype fingerprint in merged_haplotypes.fasta:
   - Assert that fingerprint length == number of positions.
   - Copy the reference sequence.
   - For each (i, ref_pos) pair, replace reference base at ref_pos with
     haplotype allele[i]:
       * Uppercase allele  -> direct substitution
       * Lowercase allele  -> substitution, softmask preserved
       * '-'  (deletion)   -> base is removed from the output sequence
4. Write reconstructed_haplotypes.fasta.
"""

import argparse
import logging
import os
import sys

import pysam


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv):
    usage = (
        "Reconstruct full-length sequences by applying merged haplotype "
        "fingerprints onto a reference sequence."
    )
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
        "--merged_haplotypes",
        dest="MERGED_HAPLOTYPES",
        type=str,
        required=True,
        help="FASTA of merged haplotype fingerprints (merged_haplotypes.fasta)",
    )
    parser.add_argument(
        "--positions",
        dest="POSITIONS",
        type=str,
        required=True,
        help="TSV with column 'position' listing the 1-based reference "
             "positions of every polymorphic site "
             "(haplotypes_filtered_positions.tsv)",
    )
    parser.add_argument(
        "--reference",
        dest="REFERENCE",
        type=str,
        required=True,
        help="FASTA file containing the reference sequence for this region. "
             "The first record is used. Positions are 1-based coordinates "
             "within this sequence (adjusted by --reference_start if needed).",
    )
    parser.add_argument(
        "--reference_start",
        dest="REFERENCE_START",
        type=int,
        default=1,
        help="1-based coordinate in the alignment reference that corresponds "
             "to position 1 of the supplied reference FASTA. Use this when "
             "the FASTA is a sub-sequence of the alignment reference. "
             "Default: 1.",
    )
    parser.add_argument(
        "-o", "--output",
        dest="OUTPUT",
        default="./",
        help="Output directory",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_reference(reference_fasta):
    """Return the sequence of the first record in *reference_fasta* (str)."""
    with pysam.FastxFile(reference_fasta) as fh:
        for record in fh:
            logging.info(
                "Using reference sequence '%s' (length %d bp)",
                record.name, len(record.sequence))
            return record.sequence
    raise ValueError(
        "Reference FASTA '{}' contains no sequences.".format(reference_fasta))


def read_positions(positions_tsv):
    """Return a sorted list of 1-based polymorphic reference positions."""
    positions = []
    with open(positions_tsv) as fh:
        header = fh.readline().strip()
        if header != "position":
            raise ValueError(
                "Expected header 'position' in '{}', got '{}'".format(
                    positions_tsv, header))
        for line in fh:
            line = line.strip()
            if line:
                positions.append(int(line))
    return sorted(positions)


def read_merged_haplotypes(merged_fasta):
    """Yield (name, fingerprint_sequence) tuples from *merged_fasta*."""
    with pysam.FastxFile(merged_fasta) as fh:
        for record in fh:
            yield record.name, record.sequence


# ---------------------------------------------------------------------------
# Reconstruction logic
# ---------------------------------------------------------------------------

def apply_haplotype(reference, positions, fingerprint, reference_start):
    """
    Substitute polymorphic positions in *reference* with alleles from
    *fingerprint* and return the reconstructed full-length sequence.

    Parameters
    ----------
    reference : str
        Full reference sequence.
    positions : list[int]
        Sorted 1-based reference positions of polymorphic sites.
        Subtract ``reference_start - 1`` to obtain 0-based indices into
        *reference*.
    fingerprint : str
        Allele string from the merged haplotype (one character per position).
        '-' encodes a deletion; lowercase encodes a softmasked base.
    reference_start : int
        1-based coordinate of reference[0] in the alignment reference.

    Returns
    -------
    str
        Reconstructed full-length sequence.
    """
    if len(fingerprint) != len(positions):
        raise ValueError(
            "Fingerprint length ({}) != number of positions ({}). "
            "Cannot map alleles to reference positions.".format(
                len(fingerprint), len(positions)))

    offset = reference_start - 1          # convert to 0-based index shift
    ref_list = list(reference)            # mutable copy

    deletion_indices = set()

    for allele, ref_pos in zip(fingerprint, positions):
        idx = ref_pos - 1 - offset        # 0-based index into ref_list
        if idx < 0 or idx >= len(ref_list):
            raise IndexError(
                "Position {} maps to index {} which is outside the "
                "reference sequence (length {}).".format(
                    ref_pos, idx, len(ref_list)))

        if allele == "-":
            deletion_indices.add(idx)
        else:
            ref_list[idx] = allele        # uppercase or lowercase preserved

    reconstructed = "".join(
        base for i, base in enumerate(ref_list)
        if i not in deletion_indices
    )
    return reconstructed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def reconstruct_haplotypes(args):
    # Guard: if the reference is the NO_FILE sentinel, skip gracefully.
    if os.path.basename(args.REFERENCE) == "NO_FILE.txt":
        logging.warning(
            "No reference sequence configured for this region -- "
            "skipping full-length reconstruction. "
            "Set params.region_references to enable this step.")
        open(os.path.join(args.OUTPUT, "reconstructed_haplotypes.fasta"), "w").close()
        return

    reference = read_reference(args.REFERENCE)
    positions  = read_positions(args.POSITIONS)
    logging.info("Loaded %d polymorphic positions.", len(positions))

    output_fasta = os.path.join(args.OUTPUT, "reconstructed_haplotypes.fasta")
    n_written = 0

    with open(output_fasta, "w") as out_f:
        for name, fingerprint in read_merged_haplotypes(args.MERGED_HAPLOTYPES):
            try:
                reconstructed = apply_haplotype(
                    reference, positions, fingerprint, args.REFERENCE_START)
            except (ValueError, IndexError) as exc:
                logging.warning("Skipping haplotype '%s': %s", name, exc)
                continue

            out_f.write(">{}\n{}\n".format(name, reconstructed))
            n_written += 1

    logging.info(
        "Wrote %d reconstructed full-length haplotypes to %s",
        n_written, output_fasta)


def main(argv=sys.argv[1:]):
    args = parse_args(argv)
    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log)
    logging.basicConfig(level=numeric_level, format="%(message)s")
    reconstruct_haplotypes(args)


if __name__ == "__main__":
    main()
