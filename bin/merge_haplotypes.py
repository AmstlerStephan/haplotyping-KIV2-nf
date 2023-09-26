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
    while queries_left:
        merged_sequences, queries_left = get_merged_sequences(unique_sequences, variant_cutoff, max_dist, stats_file_path)
        max_dist += 1
    write_haplotypes(merged_sequences, output_format, output, "merged_haplotypes")
    write_haplotype_stats(merged_sequences, output, "merged_haplotype_stats")

def get_merged_sequences(unique_sequences, variant_cutoff, max_dist, stats_file_path):
    close_sequences = find_closest_sequences(unique_sequences, variant_cutoff, max_dist, stats_file_path)
    unique_sequences = merge_sequences(unique_sequences, close_sequences)
    return unique_sequences, len(close_sequences) > 1

def merge_sequences(unique_sequences, close_sequences):
    for sequence, queries in close_sequences.items():
        for query in queries:
            unique_sequences[sequence]["reads"] = Merge(unique_sequences[sequence]["reads"], unique_sequences[query]["reads"])
    
    for sequence, queries in close_sequences.items():
        for query in queries:
            if query in unique_sequences:
                unique_sequences.pop(query)
    
    return unique_sequences

def find_closest_sequences(unique_sequences, variant_cutoff, max_dist, stats_file_path):
    n_unique_sequences = len(unique_sequences)
    cluster_cutoff = math.ceil(variant_cutoff * n_unique_sequences)
    close_sequences = dict()
    for sequence, info in unique_sequences.items():
        
        n_sequences = len(info["reads"])
        is_bigger_than_cluster_cutoff = n_sequences > cluster_cutoff
        is_high_qual = info["high_qual"]
        
        if is_bigger_than_cluster_cutoff | is_high_qual:
            for query_sequence, query_info in unique_sequences.items():
                if query_sequence == sequence:
                    continue
                
                n_queries = len(query_info["reads"])
                is_smaller_than_cluster_cutoff = n_queries <= cluster_cutoff
                is_low_qual = not query_info["high_qual"]
                
                if is_low_qual and is_smaller_than_cluster_cutoff:
                    
                    result = edlib.align(
                        sequence, 
                        query_sequence, 
                        mode="NW", 
                        task="path",
                        k=max_dist
                    )
                    if result["editDistance"] > 0:
                        if sequence in close_sequences:
                            close_sequences[sequence].append(query_sequence)
                        else:
                            close_sequences[sequence] = [query_sequence]
                        write_merge_log(sequence, query_sequence, n_queries, result, max_dist, cluster_cutoff, variant_cutoff, n_unique_sequences, stats_file_path)
    return close_sequences

   
def get_unique_sequences(fasta_file):
    unique_sequences = dict()
    
    with pysam.FastxFile(fasta_file) as reads:
        for read in reads:
            unmasked_sequence = read.sequence.upper()
            high_qual = all(base.isupper() for base in read.sequence)
            if unmasked_sequence in unique_sequences:
                unique_sequences[unmasked_sequence]["reads"][read.name] = read.sequence
            else:
                unique_sequences[unmasked_sequence] = dict()
                unique_sequences[unmasked_sequence]["reads"] = dict()
                unique_sequences[unmasked_sequence]["reads"][read.name] = read.sequence
                
            unique_sequences[unmasked_sequence]["high_qual"] = high_qual
    return unique_sequences


def write_subreads(unique_reads, output_format, output, file_name):
    haplotype_file = os.path.join(output, "{}.{}".format(file_name, output_format))
    with open(haplotype_file, "w") as out_f:
        for i, sequence in enumerate(unique_reads.keys()):
            for sub_name, sub_sequence in unique_reads[sequence]["reads"].items():
                name = "{}_{}".format(i, sub_name)
                write_fasta_read(name, sub_sequence, out_f)

                    
def write_haplotypes(haplotypes, output_format, output, file_name):
    haplotype_file = os.path.join(output, "{}.fasta".format(file_name))
    with open(haplotype_file, "w") as out_f:
        for i, sequence in enumerate(haplotypes):
            n_reads = len(haplotypes[sequence]["reads"])
            high_qual = haplotypes[sequence]["high_qual"]
            name = "{},size={},high_qual={}".format(i, n_reads, high_qual)
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


def Merge(dict1, dict2):
    res = {**dict1, **dict2}
    return res


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
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
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



def write_haplotype_stats(merged_sequences, output, file_name):
    haplotype_stats_file = os.path.join(output, "{}.tsv".format(file_name))
    
    with open(haplotype_stats_file, "w") as out_f:
        print("haplotype\thaplotype_occurences\thigh_qual\thaplotype_length", file=out_f)
        for sequence, info in merged_sequences.items():
            n_sequences = len(info["reads"])
            haplotype_length = len(sequence)
            high_qual = info["high_qual"]
            print("{}\t{}\t{}\t{}".format(sequence, n_sequences, high_qual, haplotype_length), file = out_f)


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