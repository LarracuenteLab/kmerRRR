#!/usr/bin/env python3

import os
import sys
import argparse
import time
import datetime
from collections import defaultdict
import numpy as np

def pbr_manipulate(input_file, new_threshold, output_file):
    print(f"This function is only useful to change the per base summary statistic file if the previous one was made using a repeat ratio (rr) less than that of the new threshold: {new_threshold}")
    with open(input_file, "r") as in_f:
        with open(output_file, "w") as out_f:
            count = 0
            out_f.write(f"contig\tposition\tratio\tmean\tmedian\tmode\tmax\tmin\trepetitive\tlunique\n")
            for lines in in_f:
                if lines.startswith("contig"):
                    continue
                parts = lines.strip().split("\t")
                name = parts[0]
                pos = int(parts[1])
                ratio = float(parts[2])
                mean = float(parts[3])
                median = float(parts[4])
                mode = str(parts[5])
                maximum = float(parts[6])
                minimum = float(parts[7])
                repetitive = int(parts[8])
                lunique = int(parts[9])
                if repetitive == 1 and lunique == 1 and mean >= new_threshold:
                    out_f.write(f"{name}\t{pos}\t{ratio}\t{mean}\t{median}\t{mode}\t{maximum}\t{minimum}\t{repetitive}\t{lunique}\n")
                    count += 1
                elif repetitive == 1 and lunique == 1 and mean < new_threshold:
                    lunique = 0
                    out_f.write(f"{name}\t{pos}\t{ratio}\t{mean}\t{median}\t{mode}\t{maximum}\t{minimum}\t{repetitive}\t{lunique}\n")
                else:
                    out_f.write(f"{name}\t{pos}\t{ratio}\t{mean}\t{median}\t{mode}\t{maximum}\t{minimum}\t{repetitive}\t{lunique}\n")
    print(f"{output_file} generated with {count} new manipulated regions\n")

def making_lunique_bed(input_file, cut_off_val, cut_type, bed_output_file):
    with open(bed_output_file, "w") as out_f:
        out_f.write(f'track description="BED file created from pbr file using kmerRRR" visibility=1 useScore=1\n')
        with open(input_file, "r") as in_f:
            count = 0
            contig_pos = defaultdict(list)
            for lines in in_f:
                if lines.startswith("contig"):
                    continue
                parts = lines.strip().split("\t")
                name = parts[0]
                pos = int(parts[1])
                if cut_type == "mean":
                    value = float(parts[3]) # mean
                elif cut_type == "median":
                    value = float(parts[4]) # median
                repetitive = int(parts[8])
                lunique = int(parts[9])
                if value >= cut_off_val and repetitive == 1 and lunique == 1:
                    contig_pos[name].append(pos)
            for contig_name in contig_pos:
                new_list = sorted(contig_pos[contig_name])
                data = np.array(new_list)
                breaks = np.where(np.diff(data) != 1)[0] + 1
                groups = np.split(data, breaks)
                for positions in groups:
                    if len(positions) == 1:
                        out_f.write(f"{contig_name}\t{positions[0]}\t{positions[0] + 1}\n")
                        count += 1
                    else:
                        out_f.write(f"{contig_name}\t{positions[0]}\t{positions[-1]}\n")
                        count += 1
    print(f"{bed_output_file} generated with {count} regions\n")

def making_unique_bed(input_file, cut_off_val, cut_type, bed_output_file):
    with open(bed_output_file, "w") as out_f:
        out_f.write(f'track description="BED file created from pbr file using kmerRRR" visibility=1 useScore=1\n')
        with open(input_file, "r") as in_f:
            count = 0
            contig_pos = defaultdict(list)
            for lines in in_f:
                if lines.startswith("contig"):
                    continue
                parts = lines.strip().split("\t")
                name = parts[0]
                pos = int(parts[1])
                if cut_type == "mean":
                    value = float(parts[3]) # mean
                elif cut_type == "median":
                    value = float(parts[4]) # median
                repetitive = int(parts[8])
                lunique = int(parts[9])
                if value >= cut_off_val and repetitive == 0 and lunique == 1:
                    contig_pos[name].append(pos)
            for contig_name in contig_pos:
                new_list = sorted(contig_pos[contig_name])
                data = np.array(new_list)
                breaks = np.where(np.diff(data) != 1)[0] + 1
                groups = np.split(data, breaks)
                for positions in groups:
                    if len(positions) == 1:
                        out_f.write(f"{contig_name}\t{positions[0]}\t{positions[0] + 1}\n")
                        count += 1
                    else:
                        out_f.write(f"{contig_name}\t{positions[0]}\t{positions[-1]}\n")
                        count += 1
    print(f"{bed_output_file} generated with {count} regions\n")


if __name__ == "__main__":
    start_time = time.time()
    print(f"INFO: {datetime.datetime.now()} pbr_suite started\n")
    parser = argparse.ArgumentParser(description="pbr_suite")
    parser.add_argument("-i", "--input_file", required=True, help="pbr_suite input file")
    parser.add_argument("-nt", "--new_threshold", type=float, help="new threshold value")
    parser.add_argument("-n", "--name", required=True, help="output file name")
    parser.add_argument("-pbr", "--per_base_file", action='store_true', help="Making new pbr file")
    parser.add_argument("-lu", "--lunique_file", action='store_true', help="Extracting local unique regions from pbr file in BED format")
    parser.add_argument("-u", "--unique_file", action='store_true', help="Extracting unique regions from pbr file in BED format")
    parser.add_argument("-cut", "--cut_off_val", default=1, type=float, help="cutoff value")
    parser.add_argument("-ctype", "--cut_type", default='mean', type=str, help="cutoff type")
    args = parser.parse_args()

    input_file = args.input_file
    if os.path.isfile(input_file):
        print(f"pbr_suite: {input_file} exists\n")
    else:
        print(f"pbr_suite: {input_file} does not exist\n")
        sys.exit(1)
    name = args.name

    if args.per_base_file:
        new_threshold = args.new_threshold
        if not new_threshold:
            print(f"pbr_suite: {name} has no new threshold value\nSystem exiting\n")
            sys.exit(1)
        output_file = f"{name}_new_pbr.txt"
        pbr_manipulate(input_file, new_threshold, output_file)
    else:
        print(f"pbr_suite will not create new per base file\n")

    if args.lunique_file:
        bed_file = f"{name}_new_lu_pbr.bed"
        if not args.cut_off_val:
            cut_off = args.cut_off_val
            print(f"pbr_suite will use the default value of 0.95\n")
        else:
            cut_off = args.cut_off_val
            print(f"pbr_suite will use a cut off value of {cut_off}\n")
        if not args.cut_type:
            ctype = args.cut_type
            print(f"pbr_suite will use the default cut type: mean\n")
        else:
            ctype = args.cut_type
            print(f"pbr_suite will use a cut type: {ctype}\n")
        making_lunique_bed(input_file, cut_off, ctype, bed_file)
    else:
        print(f"pbr_suite will not create local unique bed file\n")
    
    if args.unique_file:
        bed_file = f"{name}_new_unique_pbr.bed"
        if not args.cut_off_val:
            cut_off = args.cut_off_val
            print(f"pbr_suite will use the default value of 0.95\n")
        else:
            cut_off = args.cut_off_val
            print(f"pbr_suite will use a cut off value of {cut_off}\n")
        if not args.cut_type:
            ctype = args.cut_type
            print(f"pbr_suite will use the default cut type: mean\n")
        else:
            ctype = args.cut_type
            print(f"pbr_suite will use a cut type: {ctype}\n")
        making_unique_bed(input_file, cut_off, ctype, bed_file)
        
    else:
        print(f"pbr_suite will not create unique bed file\n")
    
    
    total_time = time.time() - start_time
    print(f"pbr_suite ran successfully\nTotal time required: {total_time:.2f} seconds")