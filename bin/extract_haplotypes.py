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
        required=True, 
        help="Output folder"
    )
    parser.add_argument(
        "--filter_haplotypes", dest="FILTER", action="store_true", help="Filter Haplotypes"
    )
    parser.add_argument(
        "--min_qscore",
        dest="MIN_CLUSTER_READS",
        type=int,
        default=40,
        help="Reads per cluster. Clusters with less reads will be discarded, clusters with more will be downsampled. 50% must be forward and 50% reverse reads",
    )
    
    args = parser.parse_args(argv)

    return args


def get_haplotypes(args):
    bam_file = args.BAM_FILE
    output_format = args.OUTPUT_FORMAT
    output = args.OUTPUT
    
    query_names = get_query_names(bam_file)
    haplotypes = extract_haplotypes(bam_file, query_names)
    filtered_haplotypes = filter_haplotypes(haplotypes, )
    write_haplotypes(haplotypes, output_format, output)

def get_query_names(bam_file):
    query_names = dict()
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch():
            query_names[read.query_name] = dict(
                position = list(),
                haplotype = list(), 
                quality = list())
    return query_names

def extract_haplotypes(bam_file, query_names):
    
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        # truncate = True truncates overlapping reads to positions which are aligned to the ref
        # loop over columns of bam_file
        for pileup_column in samfile.pileup(truncate = True):
            # check for polymorphic position
            if len(set(pileup_column.get_query_sequences(add_indels = True))) > 1:
                pos = pileup_column.reference_pos
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
                    
                    name = read.query_name  
                    query_names[name]["haplotype"].append(base) 
                    query_names[name]["quality"].append(qual) 
                    query_names[name]["position"].append(pos)
    return query_names


def filter_haplotypes(haplotypes, min_qscore, regions_to_exlude):
 for haplotype_name in haplotypes:
        positions = get_string(haplotypes[haplotype_name].get("haplotype"))
        for i, pos in enumerate(positions):
            if 
            
    return filtered_haplotypes
            
def write_haplotypes(haplotypes, output_format, output):
    haplotype_file = os.path.join(output, "haplotypes.{}".format(output_format))
    with open(haplotype_file, "w") as out_f:
        for haplotype_name in haplotypes:
            sequence = get_string(haplotypes[haplotype_name].get("haplotype"))
            
            if output_format == "fastq":
                qualities = get_string(haplotypes[haplotype_name].get("quality"))
                write_fastq_read(haplotype_name, sequence, qualities, out_f)
            else:
                write_fasta_read(haplotype_name, sequence, out_f) 
                
def write_fastq_read(read_name, read_seq, read_qual, out_f):
    print("@{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)
    print("+", file=out_f)
    print("{}".format(read_qual), file=out_f)

def write_fasta_read(read_name, read_seq, out_f):
    print(">{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)

def get_string(list_to_convert):
    return "".join(list_to_convert)

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