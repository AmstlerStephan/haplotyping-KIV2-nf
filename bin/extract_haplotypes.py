import argparse
import logging
import os
import re

import pysam
import pandas as pd
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
        "--variant_calling_positions", 
        dest="VARIANT_CALLING_POSITIONS", 
        type=str, 
        help="File with all positions from variant calling (column name must be 'position')"
    )
    parser.add_argument(
        "--use_variant_calling_positions", 
        dest="USE_VARIANT_CALLING_POSITIONS", 
        action="store_true", 
        help="Use positions from variant calling"
    )
    parser.add_argument(
        "-o", 
        "--output", 
        dest="OUTPUT", 
        required=True, 
        help="Output folder"
    )
    parser.add_argument(
        "--filter_haplotypes", 
        dest="FILTER_HAPLOTYPES", 
        default=True,
        action="store_true", 
        help="Filter Haplotypes"
    )
    parser.add_argument(
        "--hardmask", 
        dest="HARDMASK", 
        action="store_true", 
        help="Hardmask low quality bases"
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
        default=[{"start":0, "end":0}],
        nargs="+",
        help="Positions within the listed ranges will be excluded from the haplotypes",
    )
    parser.add_argument(
        "--variant_cutoff",
        dest="VARIANT_CUTOFF",
        type=float,
        default=0.0085,
        help="Cutoff to merge clusters",
    )
    
    
    args = parser.parse_args(argv)

    return args

def get_del():
    return "-"

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
    hardmask = args.HARDMASK
    filter_haplotypes = args.FILTER_HAPLOTYPES
    variant_cutoff = args.VARIANT_CUTOFF
    use_variant_calling_positions = args.USE_VARIANT_CALLING_POSITIONS
    variant_calling_positions = args.VARIANT_CALLING_POSITIONS
    
    query_names = get_query_names(bam_file)
    haplotypes = get_extracted_haplotypes(bam_file, query_names, variant_cutoff, use_variant_calling_positions, variant_calling_positions)
    write_haplotypes(haplotypes, output_format, output, "haplotypes")
    write_stat_file(haplotypes, output, "haplotype_stats")
    
    if filter_haplotypes:
        filtered_haplotypes = get_filtered_haplotypes(haplotypes, min_qscore, ranges_to_exclude, hardmask)
        write_haplotypes(filtered_haplotypes, output_format, output, "haplotypes_filtered")
        write_stat_file(haplotypes, output, "haplotypes_filtered_stats")
    

def get_query_names(bam_file):
    query_names = dict()
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch():
            query_names[read.query_name] = dict(
                position = list(),
                haplotype = list(), 
                quality = list())
    return query_names

def is_variant_position(pos, positions):
    return pos in positions


# Takes all bases from the current line position and counts number of bases that occur in that column 
# If there are no more than one type of base occuring, the position is not polymorphic and the bases are returned 
# Else if there are more than one type of bases it will check whether the variant base occurs more often than the set threshold 
# if so then the position is polymorphic
def is_polymorphic_position(pileup_column, variant_cutoff):
    variants = dict()
    is_polymorphic = False
    variant = None
    ## will not include the information about the next position
    bases = pileup_column.get_query_sequences(add_indels = True)
    n_bases = len(bases)
    
    for base in bases:
        if base in variants:
            variants[base] += 1
        else: 
            variants[base] = 1
    
    if len(variants) > 1:
        is_polymorphic = sum(1 for base in variants if variants[base] / n_bases >= variant_cutoff) >= 2
        # is_polymorphic = all(variants[base] / n_bases >= variant_cutoff for base in variants)

    # get most abundant base
    variant = max(variants, key=variants.get)

    # remove insertions
    if len(variant) > 1:
        variant = variant[0]
    if variant == "*":
        variant = get_del()

    return is_polymorphic, variant

def get_extracted_haplotypes(bam_file, query_names, variant_cutoff, use_variant_calling_positions, variant_calling_positions):
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        
        # Put here to load positions only once
        if use_variant_calling_positions:
            positions = pd.read_csv(variant_calling_positions, sep = "\t")["position"].to_list()
        
        # loop over columns of bam_file
        for pileup_column in samfile.pileup(min_base_quality = 0):
            # python is zero-based
            pos = pileup_column.reference_pos + 1
            polymorphic, variant = is_polymorphic_position(pileup_column, variant_cutoff)
            
            if use_variant_calling_positions and is_variant_position(pos, positions):
                    
                for pileup_read in pileup_column.pileups:
                    
                    read = pileup_read.alignment
                    name = read.query_name
                    read_pos = pileup_read.query_position
                                        
                    if not(polymorphic):
                        base = variant
                        if pileup_read.is_del:
                            base = base.lower()
                            qual = 70
                        else:
                            qual = read.query_qualities[read_pos]
                    # deletion is annotated as - and qual set to 70 (might adjust)
                    elif pileup_read.is_del:
                        base = get_del()
                        qual = 70
                    else:
                        base = read.query_sequence[read_pos]
                        qual = read.query_qualities[read_pos]

                    query_names[name]["haplotype"].append(base) 
                    query_names[name]["quality"].append(qual) 
                    query_names[name]["position"].append(pos)

            elif polymorphic and not(use_variant_calling_positions):
                # loop over reads per column
                for pileup_read in pileup_column.pileups:
                    read = pileup_read.alignment
                    name = read.query_name
                    read_pos = pileup_read.query_position

                    # in case of indel a list is returned and joined with the existing list
                    if pileup_read.indel >= 1:
                        indel_start = read_pos
                        indel_end = read_pos + pileup_read.indel + 1
                        
                        bases = read.query_sequence[indel_start:indel_end]
                        quals = read.query_qualities[indel_start:indel_end]
                        
                        query_names[name]["haplotype"].append(bases)
                        query_names[name]["quality"].append(quals)
                        query_names[name]["position"].append(pos)
                        continue
                    
                    # deletion is annotated as - and qual set to 70 (might adjust)
                    if pileup_read.is_del:
                        base = get_del()
                        qual = 70
                    else:
                        base = read.query_sequence[read_pos]
                        qual = read.query_qualities[read_pos]
                        
                    query_names[name]["haplotype"].append(base) 
                    query_names[name]["quality"].append(qual) 
                    query_names[name]["position"].append(pos)

    return query_names


def get_filtered_haplotypes(haplotypes, min_qscore, regions_to_exclude, hardmask):
    for haplotype_name in haplotypes:
        positions = haplotypes[haplotype_name].get("position").copy()
        for position in positions:
            i = haplotypes[haplotype_name].get("position").index(position)
            qual = get_quality(haplotypes[haplotype_name].get("quality")[i])
            
            if exclude_pos(position, regions_to_exclude):
                haplotypes[haplotype_name]["position"].pop(i)
                haplotypes[haplotype_name]["haplotype"].pop(i)
                haplotypes[haplotype_name]["quality"].pop(i)
            # elif qual < min_qscore:
            #    if hardmask:
            #        haplotypes[haplotype_name]["haplotype"][i] = "N"
            #    else:
            #        base = haplotypes[haplotype_name]["haplotype"][i]
            #        masked_base = base.lower()
            #        haplotypes[haplotype_name]["haplotype"][i] = masked_base
                
    return haplotypes
            
def write_haplotypes(haplotypes, output_format, output, file_name):
    haplotype_file = os.path.join(output, "{}.{}".format(file_name, output_format))
    with open(haplotype_file, "w") as out_f:
        for haplotype_name in haplotypes:
            sequence = get_string(haplotypes[haplotype_name].get("haplotype"), "")
            if output_format == "fastq":
                qualities = get_quality_string(haplotypes[haplotype_name].get("quality"), "")
                write_fastq_read(haplotype_name, sequence, qualities, out_f)
            else:
                write_fasta_read(haplotype_name, sequence, out_f)

def write_stat_file(haplotypes, output, file_name):
    haplotype_file = os.path.join(output, "{}.tsv".format(file_name))
    stats = get_stats(haplotypes)
    with open(haplotype_file, "w") as out_f:
        print("pos\tbase\tcount", file=out_f)
        for pos in stats:
            for base in stats[pos]:
                print("{}\t{}\t{}".format(pos, base, stats[pos][base]), file=out_f)
            
def get_stats(haplotypes):
    position_stats = dict()
    for haplotype_name in haplotypes:
        # loop over every position of each haplotype
        for i, pos in enumerate(haplotypes[haplotype_name]["position"]):
            base = str(haplotypes[haplotype_name]["haplotype"][i])
            if pos in position_stats:
                if base in position_stats[pos]:
                    position_stats[pos][base] += 1
                else:
                    position_stats[pos][base] = 1
            else:
                position_stats[pos] = dict()
                position_stats[pos][base] = 1
                
    return position_stats
    
def write_fastq_read(read_name, read_seq, read_qual, out_f):
    print("@{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)
    print("+", file=out_f)
    print("{}".format(read_qual), file=out_f)

def write_fasta_read(read_name, read_seq, out_f):
    print(">{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)

def get_string(list_to_convert, sep):
    return sep.join(map(str, list_to_convert))

def get_quality_string(qualities, sep):
    qualities_parsed = list()
    for qual in qualities:
        if not isinstance(qual, int):
            for indel_qual in qual:
                qualities_parsed.append(chr(indel_qual + 33))
            continue
        else:
            qualities_parsed.append(chr(qual + 33))
                    
    return get_string(qualities_parsed, sep) 

def get_quality(qual):
    if not isinstance(qual, int):
        qual = min(qual)
    return qual

def exclude_pos(pos, regions_to_exclude):
    exclude = False
    for range in regions_to_exclude:
        exclude = pos >= range["start"] and pos <= range["end"] or exclude
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