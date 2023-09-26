import argparse
import logging
import os
import math
import re
import sys

import pysam
import edlib

def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to telemap"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
            "CRITICAL",
            "debug",
            "info",
            "warning",
            "error",
            "critical",
        ],
        default="INFO",
        help="Print debug information",
    )
    parser.add_argument(
        "--fastx_file", 
        dest="FASTX_FILE", 
        type=str, 
        required=True, 
        help="input fastx file to extract haplotypes from"
    )
    parser.add_argument(
        "--output_format", 
        dest="OUTPUT_FORMAT", 
        type=str, 
        default="fasta",
        help="Output format of the haplotypes"
    )
    parser.add_argument(
        "--variant_cutoff",
        dest="VARIANT_CUTOFF",
        type=float,
        default=0.0085,
        help="Cutoff to merge clusters",
    )
    parser.add_argument(
        "-o", 
        "--output", 
        dest="OUTPUT",
        default="./",
        help="Output folder"
    )
    
    args = parser.parse_args(argv)

    return args


def get_merged_haplotypes(args):
    fasta_file = args.FASTX_FILE
    output_format = args.OUTPUT_FORMAT
    output = args.OUTPUT
    variant_cutoff = args.VARIANT_CUTOFF
    queries_left = True
    max_dist = 1
    stats_file_name = "merged_haplotype_log"
    stats_file_path = os.path.join(output, "{}.tsv".format(stats_file_name))
    
    # write first line of log
    with open(stats_file_path, "w") as stats_file:
        print("sequence\tquery_sequence\tquery_size\tedist\tposition\tbase\tquery_base\tchange\tmax_dist\tcluster_cutoff\tvariant_cutoff\tn_unique_sequences", file = stats_file)
    
    unique_sequences = get_unique_sequences(fasta_file)
    write_haplotypes(unique_sequences, output_format, output, "unique_haplotypes")
    write_subreads(unique_sequences, output_format, output, "unique_haplotypes_subreads")
    while queries_left | max_dist < 100:
        merged_sequences, queries_left = get_merged_sequences(unique_sequences, variant_cutoff, max_dist, stats_file_path)
        max_dist += 1
    write_haplotypes(merged_sequences, output_format, output, "merged_haplotypes")
    write_haplotype_stats(merged_sequences, output, "merged_haplotype_stats")

def write_haplotype_stats(merged_sequences, output, file_name):
    haplotype_stats_file = os.path.join(output, "{}.tsv".format(file_name))
    
    with open(haplotype_stats_file, "w") as out_f:
        print("haplotype\thaplotype_occurences\thigh_qual\thaplotype_length", file=out_f)
        for sequence in merged_sequences:
            n_sequences = len(merged_sequences[sequence]) - 1
            haplotype_length = len(sequence)
            high_qual = merged_sequences[sequence]["high_qual"]
            print("{}\t{}\t{}\t{}".format(sequence, n_sequences, high_qual, haplotype_length), file = out_f)

def get_merged_sequences(unique_sequences, variant_cutoff, max_dist, stats_file_path):
    n_unique_sequences = len(unique_sequences)
    
    # value of "high_qual" is also included in len
    cluster_cutoff = math.ceil(variant_cutoff * n_unique_sequences) + 1
    sequences_to_remove = list()
    unique_sequences_copy = unique_sequences.copy()
    queries_left = False
    for sequence in unique_sequences_copy:
        for query_sequence in unique_sequences:
            n_queries = len(unique_sequences[query_sequence])
            if not unique_sequences[query_sequence]["high_qual"] and n_queries <= cluster_cutoff:
                result = edlib.align(
                    sequence, 
                    query_sequence, 
                    mode="NW", 
                    task="path",
                    k=max_dist
                )
                if result["editDistance"] > 0:
                    queries_left = True
                    sequences_to_remove.append(query_sequence)
                    write_merge_log(sequence, query_sequence, n_queries, result, max_dist, cluster_cutoff, variant_cutoff, n_unique_sequences, stats_file_path)
        
        for sequence_to_remove in sequences_to_remove:
            unique_sequences[sequence_to_remove] = Merge(unique_sequences[sequence_to_remove], unique_sequences[sequence_to_remove])
            unique_sequences.pop(sequence_to_remove)
        sequences_to_remove.clear()
    return unique_sequences, queries_left

def write_merge_log(sequence, query_sequence, n_queries, result, max_dist, cluster_cutoff, variant_cutoff, n_unique_sequences, stats_file_path):
    edist = result["editDistance"]
    cigar = result["cigar"]
    with open(stats_file_path, "a+") as stats_file:
        for difference in re.findall('\d*=..', cigar):
            change = difference.split("=")[1]
            n_bases = int(re.findall("\d*", change)[0])
            pos = int(difference.split("=")[0])
            base = sequence[pos:pos+n_bases]
            query_base = query_sequence[pos:pos+n_bases]
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                sequence, 
                query_sequence,
                n_queries,
                edist,
                pos,
                base,
                query_base,
                change,
                max_dist,
                cluster_cutoff,
                variant_cutoff,
                n_unique_sequences),
                  file = stats_file)

def get_unique_sequences(fasta_file):
    unique_sequences = dict()
    
    with pysam.FastxFile(fasta_file) as reads:
        for read in reads:
            unmasked_sequence = read.sequence.upper()
            contains_lower = not any(base.islower() for base in read.sequence)
            if unmasked_sequence in unique_sequences:
                unique_sequences[unmasked_sequence].update({ read.name : read.sequence})
            else:
                unique_sequences[unmasked_sequence] = {read.name : read.sequence}
                
            unique_sequences[unmasked_sequence].update({ "high_qual" : contains_lower})
    return unique_sequences

def write_subreads(unique_reads, output_format, output, file_name):
    offset = 0
    haplotype_file = os.path.join(output, "{}.{}".format(file_name, output_format))
    with open(haplotype_file, "w") as out_f:
        for i, sequences in enumerate(unique_reads):
            for j, sequence in enumerate(unique_reads[sequences]):
                if sequence != "high_qual":
                    name = "{}_{}".format(i, j - offset)
                    write_fasta_read(name, unique_reads[sequences][sequence], out_f)
                else:
                    offset = 1 
            offset = 0
                    
def write_haplotypes(haplotypes, output_format, output, file_name):
    haplotype_file = os.path.join(output, "{}.{}".format(file_name, output_format))
    with open(haplotype_file, "w") as out_f:
        for i, sequence in enumerate(haplotypes):
            name = "{},size={},high_qual={}".format(i, len(haplotypes[sequence]), haplotypes[sequence]["high_qual"])
            write_fasta_read(name, sequence, out_f) 

def write_fastq_read(read_name, read_seq, read_qual, out_f):
    # print("@{},positions={}".format(read_name, positions), file=out_f)
    print("@{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)
    print("+", file=out_f)
    print("{}".format(read_qual), file=out_f)

def write_fasta_read(read_name, read_seq, out_f):
    # print(">{},positions={}".format(read_name, positions,), file=out_f)
    print(">{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)

# Python code to merge dict using a single
# expression
def Merge(dict1, dict2):
    res = {**dict1, **dict2}
    return res

def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log.upper())
    logging.basicConfig(level=numeric_level, format="%(message)s")

    get_merged_haplotypes(args)


if __name__ == "__main__":
    main()