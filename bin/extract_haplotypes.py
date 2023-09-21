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
        "--bam_file", 
        dest="BAM_FILE", 
        type=str, 
        required=True, 
        help="input bam file to extract haplotypes from"
    )
    
    args = parser.parse_args(argv)

    return args


def get_haplotypes(args):
    bam_file = args.BAM_FILE
    
    query_names = get_query_names(bam_file)
    haplotypes = extract_haplotypes(bam_file)
    parsed_haplotypes = parse_haplotypes(haplotypes, query_names)
    print(parsed_haplotypes)

def get_query_names(bam_file):
    query_names = dict()
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch():
            query_names[read.query_name] = dict(haplotype = list(), qual = list())
    print(query_names)
    return query_names

def extract_haplotypes(bam_file):
    position_dict = dict()
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        # truncate = True truncates overlapping reads to positions which are aligned to the ref
        # loop over columns of bam_file
        for pileup_column in samfile.pileup(truncate = True):
            # check for polymorphic position
            if len(set(pileup_column.get_query_sequences(add_indels = True))) > 1:
                pos = pileup_column.reference_pos
                position_dict[pos] = list()
            # loop over reads per column
                for pileup_read in pileup_column.pileups:
                    read = pileup_read.alignment
                    read_pos = pileup_read.query_position
                    if pileup_read.is_del:
                        base = "D"
                        qual = 70
                    elif( pileup_read.indel > 0):
                        base = read.query_sequence[read_pos:read_pos + pileup_read.indel]
                        qual = read.query_qualities[read_pos:read_pos + pileup_read.indel]
                    else:
                        base = read.query_sequence[read_pos]
                        qual = read.query_qualities[read_pos]
                        
                    
                    read_info = dict(
                        name = read.query_name, 
                        base = base, 
                        qual = qual)
                    position_dict[pos].append(read_info)

    return position_dict


def parse_haplotypes(haplotypes, query_names):
    for pos, reads in sorted(haplotypes.items()):
        parsed_haplotypes = dict()

        for read in reads:
            name = read.get("name")
            base = read.get("base")
            qual = read.get("qual")
            query_names[name]["haplotype"].append(base) 
            query_names[name]["qual"].append(qual) 
            
    return query_names
            

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

    get_haplotypes(args)


if __name__ == "__main__":
    main()