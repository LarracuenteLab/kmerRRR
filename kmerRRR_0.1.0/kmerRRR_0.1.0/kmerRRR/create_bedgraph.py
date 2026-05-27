#!/usr/bin/env python3

import sys
import datetime
import time
import os

def create_bed(pbr_file, slide_window, out_file):
    with open(pbr_file, "r") as in_file:
        mean_dict = {}
        median_dict = {}
        for lines in in_file:
            pos_mean = {}
            pos_median = {}
            if lines.startswith("contig"):
                continue
            parts = lines.strip().split()
            contig_name = parts[0]
            position = int(parts[1])
            mean = float(parts[3])
            median = float(parts[4])
            pos_mean[position] = mean
            pos_median[position] = median
            if contig_name in mean_dict:
                mean_dict[contig_name].update(pos_mean)
            else:
                mean_dict[contig_name] = pos_mean
            if contig_name in median_dict:
                median_dict[contig_name].update(pos_median)
            else:
                median_dict[contig_name] = pos_median
    mean_out = f"{out_file}_mean.bedgraph"
    median_out = f"{out_file}_median.bedgraph"
    
    with open(mean_out, "w") as mean_out:
        mean_out.write(
            f'track type=bedGraph name="kmerRRR" description="kmerRRR Mean Ratio" visibility=1 color=0,255,0\n')
        for contig in mean_dict:
            pos_val = mean_dict[contig]
            for i in  range(min(pos_val.keys()), max(pos_val.keys()), slide_window):
                mean_temp = 0
                for j in range(i, i + slide_window):
                    if j in pos_val:
                        mean_temp += pos_val[j]
                mean_to_write = mean_temp / slide_window
                mean_out.write(f"{contig}\t{i}\t{i + slide_window}\t{mean_to_write}\n")

    with open(median_out, "w") as median_out:
        median_out.write(f'track type=bedGraph name="kmerRRR" description="kmerRRR Mean Ratio" visibility=1 color=0,0,255\n')
        for contig in median_dict:
            pos_val = median_dict[contig]
            for i in range(min(pos_val.keys()), max(pos_val.keys()), slide_window):
                median_temp = 0
                for j in range(i, i + slide_window):
                    if j in pos_val:
                        median_temp += pos_val[j]
                median_to_write = median_temp / slide_window
                median_out.write(f"{contig}\t{i}\t{i + slide_window}\t{median_to_write}\n")

def create_bed_mapq(mapq_file, name):
    mapq = []
    luc = []
    with open(mapq_file, 'r') as pbr_f:
        for lines in pbr_f:
            if lines.startswith("contig"):
                continue
            parts = lines.strip().split()
            luc_bool = int(parts[4])
            mapq.append(lines)
            if luc_bool == 1:
                luc.append(lines)
    mapq_out = f"{name}.mapq.bedgraph"
    luc_out = f"{name}.luc.bedgraph"
    with open(mapq_out, 'w') as map_out:
        map_out.write('track type=bedGraph name="MAPQ" description="MAPQ Scores" visibility=1 color=255,165,0\n')
        for line in mapq:
            map_out.write(line)
    if len(luc) > 0:
        with open(luc_out, 'w') as luc_of:
            luc_of.write('track type=bedGraph name="LUC" description="MAPQ Scores for LUC" visibility=1 color=126,209,243\n')
            for line in luc:
                luc_of.write(line)
    else:
        print(f"NO LUC found.\nNO LUC bedgrah will be generated.\n")

def main(args):
    print(f"{datetime.datetime.now()} INFO: create_bed has started\n")
    start_time = time.time()
    if args.per_base_ratio:
        pbr_file = args.per_base_ratio
    slide_window = args.slide_window
    out_file = args.name
    if args.per_base_ratio:
        pbr_file = args.per_base_ratio
        if os.path.exists(pbr_file):
            print(f"Per base ratio file: {pbr_file} found.\n")
        else:
            print(f"Per base ration file: {pbr_file} not found.\nSystem exiting...\n")
            sys.exit(1)
        print(f"Output will be written in {out_file} for mean and median.\n")
        create_bed(pbr_file, slide_window, out_file)
    if args.mapq_file:
        mapq_file = args.mapq_file
        print(f"Creating bedgraph file formats for MAPQ and LUC\n")
        create_bed_mapq(mapq_file, out_file)
    print(f"create_bed ran successfully")
    print(f"Processing time: {time.time() - start_time:.2f} seconds")