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
        "--cluster_stats",
        dest="CLUSTER_STATS",
        type=str,
        help="cluster stats file"
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="OUTPUT",
        required=True,
        help="Output folder"
    )
    parser.add_argument(
        "--min_reads_per_cluster",
        dest="MIN_CLUSTER_SIZE",
        type=int,
        help="Minimal cluster size"
    )
    parser.add_argument(
        "--max_reads_per_cluster",
        dest="MAX_CLUSTER_SIZE",
        type=int,
        help="Maximal cluster size"
    )

    args = parser.parse_args(argv)

    return args

def filter_bam(args):
    bam_file = args.BAM_FILE
    cluster_stats_file = args.CLUSTER_STATS
    min_cluster_size = args.MIN_CLUSTER_SIZE
    max_cluster_size = args.MAX_CLUSTER_SIZE
    output = args.OUTPUT
    outfile_name = outfile_name = os.path.join(output, "filtered_bam.bam")
    clusters_above_threshold = get_clusters(cluster_stats_file, min_cluster_size, max_cluster_size)

    # Load cluster stats DataFrame for read tagging
    cluster_stats = pd.read_csv(cluster_stats_file, sep="\t").copy()

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        outfile = pysam.AlignmentFile(outfile_name, "w", template = bam)
        for read in bam.fetch(until_eof=True):
            name = read.query_name[:-2]
            if read.query_name in clusters_above_threshold:
                # Find matching cluster stats
                cluster_match = cluster_stats[cluster_stats["cluster_id"] == read.query_name]
                if not cluster_match.empty:
                    # add information about cluster size, found reads, written reads fwd and rev to read tag
                    reads_found = int(cluster_match["reads_found"].values[0])
                    reads_written_fwd = int(cluster_match["reads_written_fwd"].values[0])
                    reads_written_rev = int(cluster_match["reads_written_rev"].values[0])

                    read.set_tag("RF", reads_found)
                    read.set_tag("CS", reads_written_fwd + reads_written_rev)
                    read.set_tag("WF", reads_written_fwd)  # Changed from "RF" to "WF" to avoid duplicate
                    read.set_tag("RR", reads_written_rev)
                outfile.write(read)

def get_clusters(cluster_stats_file, min_cluster_size, max_cluster_size):
    cluster_stats = pd.read_csv(cluster_stats_file, sep = "\t")

    cluster_stats_filtered = cluster_stats[(cluster_stats["cluster_written"] == 1) &
                                           (cluster_stats["reads_written_fwd"] + cluster_stats["reads_written_rev"] >= min_cluster_size) &
                                           (cluster_stats["reads_found"] <= max_cluster_size)].copy()

    cluster_stats_filtered_parsed = cluster_stats_filtered["cluster_id"].to_list()

    return cluster_stats_filtered_parsed

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

    filter_bam(args)


if __name__ == "__main__":
    main()
