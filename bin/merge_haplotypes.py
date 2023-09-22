import argparse
import logging
import os

import pysam
import sys
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
    
    unique_sequences = get_unique_sequences(fasta_file)
    write_haplotypes(unique_sequences, output_format, output, "unique_haplotypes")
    write_subreads(unique_sequences, output_format, output, "unique_haplotypes_subreads")
    merged_sequences = get_merged_sequences(unique_sequences)
    write_haplotypes(merged_sequences, output_format, output, "merged_haplotypes")

def get_merged_sequences(unique_sequences):
    sequences_to_remove = list()
    unique_sequences_copy = unique_sequences.copy()
    for sequence in unique_sequences_copy:
        for query_sequence in unique_sequences:
            if not unique_sequences[query_sequence]["high_qual"] and len(unique_sequences[query_sequence]) < 3:
                result = edlib.align(
                    sequence, 
                    query_sequence, 
                    mode="HW", 
                    task="location",
                    k=1
                )
                if result["editDistance"] > 0:
                    sequences_to_remove.append(query_sequence)
        
        for sequence in sequences_to_remove:
            unique_sequences[sequence] = Merge(unique_sequences[sequence], unique_sequences[query_sequence])
            unique_sequences.pop(sequence)
        sequences_to_remove.clear()
    return unique_sequences



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

def write_fastq_read(read_name, positions, read_seq, read_qual, out_f):
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