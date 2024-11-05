#!/usr/bin/env python3

import sys
import gzip
from collections import defaultdict


def process_pair(line_1, line_2, read_repeat_thresh, motif_repeat_thresh, motif_filter, motif_size_filter):
    col_1 = defaultdict(int)
    for i in range(1, len(line_1)):
        motif_res = line_1[i].split(":")
        col_1[motif_res[0]] += int(motif_res[1])
    
    if sum(col_1.values()) < read_repeat_thresh:
        return None

    col_2 = defaultdict(int)
    for i in range(1, len(line_2)):
        motif_res = line_2[i].split(":")
        col_2[motif_res[0]] += int(motif_res[1])
    
    if sum(col_2.values()) < read_repeat_thresh:
        return None

    # combine repeat counts for both reads
    for k, v in col_2.items():
        col_1[k] += v

    filter_motif = [_ for _ in col_1.keys() if _ not in motif_filter] if motif_filter is not None else []
    filter_motif_size = [_ for _ in col_1.keys() if len(_) not in motif_size_filter] if motif_size_filter is not None else []

    # Do not filter motifs based on size if specifically specified
    if motif_filter is not None and motif_size_filter is not None:
        filter_motif_size = [_ for _ in filter_motif_size if _ not in motif_filter]

    filter_motif_thresh = [k for k in col_1.keys() if col_1[k] < motif_repeat_thresh]
    for m in set(filter_motif + filter_motif_size + filter_motif_thresh):
        col_1.pop(m)

    if len(col_1) == 0:
        return None

    # check whether there is still sufficient repetitive reads to report pair after filtering
    if sum(col_1.values()) < read_repeat_thresh:
        return None

    return "-".join(sorted(col_1.keys()))


def parse_pair_reads(input_path, read_repeat_thresh, motif_repeat_thresh, motif_filter, motif_size_filter):
    motif_counts = defaultdict(int)
    with (gzip.open if input_path.endswith(".gz") else open)(input_path, 'rt') as in_file:
        total_line = ""
        last_line = in_file.readline().strip().split()
        for line in in_file:
            if "Total" in line:
                total_line = line.strip()
                continue
            current_line = line.strip().split()
            if current_line[0] == last_line[0]:
                last_pair = last_line
                current_pair = current_line
                pair_motifs = process_pair(last_line, current_line, read_repeat_thresh, motif_repeat_thresh, motif_filter, motif_size_filter)
                if pair_motifs is not None:
                    motif_counts[pair_motifs] += 1
            last_line = current_line
    return motif_counts


def write_output(output_path, motif_counts):
    # Sort entries in decreasing order
    sorted_motif_counts = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)
    with open(output_path, 'w') as out_file:
        out_file.write("Expanded motifs\tNumber of read pairs\n")
        for k, v in sorted_motif_counts:
            out_file.write(f"{k}\t{v}\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    required_group = parser.add_argument_group(title="required arguments")
    required_group.add_argument("-i", "--input", dest="input_path",
                                default=None,
                                action="store",
                                help="Path to sorted superSTR per_read output", 
                                required=True)
    required_group.add_argument("-o", "--output", dest="output_path",
                                default=None,
                                action="store",
                                help="Output path for files.", 
                                required=True)
    optional_parameters = parser.add_argument_group(title="optional parameters")
    optional_parameters.add_argument("--read-repeat-thresh", action="store",
                                    dest="read_repeat_thresh", default=120, type=int,
                                    help="Minimum number of base pairs of repetive sequence required in each read to include read pair in analysis.")
    optional_parameters.add_argument("--motif-repeat-thresh", action="store",
                                    dest="motif_repeat_thresh", default=50, type=int,
                                    help="Minimum number of base pairs of repetive sequence in a read pair to include motif")
    optional_parameters.add_argument("--motif", action="store",
                                    dest="motif_filter", default=None,
                                    help="Comma separate set of motifs to count (if specified then any other motifs will not be reported)")
    optional_parameters.add_argument("--motif-size", action="store",
                                    dest="motif_size_filter", default=None,
                                    help="Comma separated set of motif lengths in base pairs to keep (if then any other motifs will not be reported, unless specified with --motif argument)")
    # Default to help message if no arguments provided.
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    if args.motif_filter is not None:
        args.motif_filter = [_.strip() for _ in args.motif_filter.split(",")]
    if args.motif_size_filter is not None:
        args.motif_size_filter = [int(_) for _ in args.motif_size_filter.split(",")]

    motif_counts = parse_pair_reads(args.input_path, 
                                    args.read_repeat_thresh, 
                                    args.motif_repeat_thresh, 
                                    args.motif_filter, 
                                    args.motif_size_filter)
    write_output(args.output_path, motif_counts)

