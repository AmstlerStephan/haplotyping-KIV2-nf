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
        dest="MIN_QSCORE",
        type=int,
        default=40,
        help="Positions with less quality will be discarded",
    )
    parser.add_argument(
        "--ranges_to_exclude",
        dest="RANGES_TO_EXCLUDE",
        type=ranges,
        nargs="+",
        help="Positions within the listed ranges will be excluded from the haplotypes",
    )
    
    args = parser.parse_args(argv)

    return args

def ranges(query_range):
    range = dict()
    try:
        start, end = map(int, query_range.split(','))
        range["start"]=start
        range["end"]=end
        return range
    except:
        raise argparse.ArgumentTypeError("range must be start,end")

def get_haplotypes(args):
    bam_file = args.BAM_FILE
    output_format = args.OUTPUT_FORMAT
    output = args.OUTPUT
    min_qscore = args.MIN_QSCORE
    ranges_to_exclude = args.RANGES_TO_EXCLUDE
    
    query_names = get_query_names(bam_file)
    haplotypes = extract_haplotypes(bam_file, query_names)
    write_haplotypes(haplotypes, output_format, output, "haplotypes")
    filtered_haplotypes = filter_haplotypes(haplotypes, min_qscore, ranges_to_exclude)
    write_haplotypes(filtered_haplotypes, output_format, output, "haplotypes_filtered")

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
                    name = read.query_name  
                    
                    if pileup_read.indel > 0:
                        indel_length = pileup_read.indel
                        base = read.query_sequence[read_pos:read_pos + indel_length]
                        qual = read.query_qualities[read_pos:read_pos + indel_length]
                        query_names[name]["haplotype"].extend(base) 
                        query_names[name]["quality"].extend(qual)
                        query_names[name]["position"].extend([pos] * indel_length)
                        continue
                    
                    if pileup_read.is_del:
                        base = "D"
                        qual = 70
                    else:
                        base = read.query_sequence[read_pos]
                        qual = read.query_qualities[read_pos]
                    
                    query_names[name]["haplotype"].append(base) 
                    query_names[name]["quality"].append(qual) 
                    query_names[name]["position"].append(pos)
    return query_names


def filter_haplotypes(haplotypes, min_qscore, regions_to_exlude):
    for haplotype_name in haplotypes:
        positions = haplotypes[haplotype_name].get("position").copy()
        for pos in positions:
            i = haplotypes[haplotype_name].get("position").index(pos)
            qual = haplotypes[haplotype_name].get("quality")[i]
            if exclude_pos(pos, regions_to_exlude):
                haplotypes[haplotype_name].get("position").pop(i)
                haplotypes[haplotype_name].get("haplotype").pop(i)
                haplotypes[haplotype_name].get("quality").pop(i)
            elif qual < min_qscore:
                base = haplotypes[haplotype_name].get("haplotype")[i]
                masked_base = base.lower()
                haplotypes[haplotype_name].get("haplotype")[i] = masked_base
                
    return haplotypes
            
def write_haplotypes(haplotypes, output_format, output, file_name):
    haplotype_file = os.path.join(output, "{}.{}".format(file_name, output_format))
    with open(haplotype_file, "w") as out_f:
        for haplotype_name in haplotypes:
            sequence = get_string(haplotypes[haplotype_name].get("haplotype"), "")
            positions = get_string(haplotypes[haplotype_name].get("position"), " ")
            if output_format == "fastq":
                qualities = get_string(haplotypes[haplotype_name].get("quality"), "")
                write_fastq_read(haplotype_name, positions, sequence, qualities, out_f)
            else:
                write_fasta_read(haplotype_name, positions, sequence, out_f) 
                
def write_fastq_read(read_name, positions, read_seq, read_qual, out_f):
    print("@{},positions={}".format(read_name, positions), file=out_f)
    print("{}".format(read_seq), file=out_f)
    print("+", file=out_f)
    print("{}".format(read_qual), file=out_f)

def write_fasta_read(read_name, positions, read_seq, out_f):
    print(">{},positions={}".format(read_name, positions,), file=out_f)
    print("{}".format(read_seq), file=out_f)

def get_string(list_to_convert, sep):
    return sep.join(map(str, list_to_convert))

def exclude_pos(pos, regions_to_exclude):
    exclude = False
    for range in regions_to_exclude:
        exclude = pos > range["start"] and pos < range["end"] or exclude
    return exclude

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