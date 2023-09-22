import argparse
import logging
import os

import pysam
import sys

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
                
    print(unique_sequences)
    print(len(unique_sequences))

def get_unique_sequences(fasta_file):
    unique_sequences = dict()
    with pysam.FastxFile(fasta_file) as reads:
        for read in reads:
            unmasked_sequence = read.sequence.upper()
            if unmasked_sequence in unique_sequences:
                unique_sequences[unmasked_sequence].update({ read.name : read.sequence})
            else:
                unique_sequences[unmasked_sequence] = {read.name : read.sequence}
    return unique_sequences        
        
def write_fastq_read(read_name, positions, read_seq, read_qual, out_f):
    # print("@{},positions={}".format(read_name, positions), file=out_f)
    print("@{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)
    print("+", file=out_f)
    print("{}".format(read_qual), file=out_f)

def write_fasta_read(read_name, positions, read_seq, out_f):
    # print(">{},positions={}".format(read_name, positions,), file=out_f)
    print(">{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)


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