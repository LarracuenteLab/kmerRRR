#!/usr/bin/env python3

"""This script is to manipulate MAPQ scores in a given 
bamfile based on a per-base ratio file and user-defined cut-off
values. It will read the per-base ratio file
and identify positions that meet the cut-off criteria,
then replace the MAPQ scores of reads in the bamfile meeting
the criteria with MAPQL (for MAPQ Local adjusted scores). 
The default MAPQL score is 40, but users can specify different 
value. Default cut-off values for the per-base ratios are 
set to 0.75 for mean and 0.85 for median. The output is a new 
bamfile that includes the manipulated MAPQL scores. 

User needs to provide: 
1. a per-base ratio file
2. a bam file
3. a name for the output bam file
"""
import os
import sys
import pysam 
from collections import defaultdict
import time
import datetime
import random

def parsing_ratio_file(per_base_ratio_file, cut_type):
    with open(per_base_ratio_file, 'r') as fratio:
        contig_positions = defaultdict(set)
        position_ctype = {}
        for lines in fratio:
            if lines.startswith("contig"):
                continue
            parts = lines.strip().split("\t")
            contig_name = parts[0]
            position = int(parts[1]) - 1 # pysam's and python's 0-based system
            if cut_type == "mean":
                ctype = float(parts[3])
            elif cut_type == "median":
                ctype = float(parts[4])
            contig_positions[contig_name].add(position)
            position_ctype[contig_name, position] = (ctype)
    return contig_positions, position_ctype

def parsing_bamfiles(bamfile, per_base_ratio_file, cut_off, cut_type, old_mapq, new_mapq, output_bamfile):
    seed_value = int.from_bytes(os.urandom(8), byteorder='big')
    random.seed(seed_value)
    print(f"Seed value used: {seed_value}")
    total_reads = skipped_unmap = skipped_H = manipulated = high_score = not_manipulated = 0
    contig_positions, position_ctype = parsing_ratio_file(per_base_ratio_file, cut_type)
    if not contig_positions:
        print("No contig positions found in the ratio file.")
        sys.exit(1)

    with pysam.AlignmentFile(bamfile, "rb") as bam, pysam.AlignmentFile(output_bamfile, "wb", template=bam) as outbam:
        read_name = defaultdict(list)
        high_score_reads = {}
        for read in bam.fetch(until_eof=True):
            new_read_name = f"{read.query_name}_{read.flag}"
            total_reads += 1
            if read.is_unmapped:
                outbam.write(read)
                skipped_unmap += 1
                continue
            if read.reference_name in contig_positions:
                if read.mapping_quality > old_mapq:
                    high_score += 1
                    outbam.write(read)
                    high_score_reads[new_read_name] = 1
                    if new_read_name in read_name:
                        del read_name[new_read_name]
                elif 'H' in read.cigarstring:
                    outbam.write(read)
                    skipped_H += 1
                    continue
                else:
                    start_pos = read.reference_start
                    end_pos = read.reference_end
                    manipulate = manipulate_mapq(read.reference_name, start_pos, end_pos, contig_positions, position_ctype, cut_off) if any(pos in contig_positions[read.reference_name] for pos in range(start_pos, end_pos)) else False
                    if manipulate:
                        if new_read_name in high_score_reads:
                            continue
                        read.mapping_quality = new_mapq
                        manipulated +=1
                        read.set_tag("KR", "LUC", value_type='Z')
                        read_name[new_read_name].append(read)
                    else:
                        not_manipulated +=1
                        outbam.write(read)
        for qnames in read_name:
            if len(read_name[qnames]) == 1:
                outbam.write(read_name[qnames][0])
                continue
            random_read = random.choice(read_name[qnames])
            outbam.write(random_read)
        
    print(f"Total reads: {total_reads}, Skipped unmapped reads: {skipped_unmap}, Skipped hardclipped reads: {skipped_H}, Manipulated reads: {manipulated}, Not manipulated reads: {not_manipulated}, Reads with high MAPQ: {high_score}")

def manipulate_mapq(contig, start_pos, end_pos, contig_positions, position_ctype, cut_off):
    if not any((contig, pos) in position_ctype for pos in range(start_pos, end_pos)):
        return False
    mean_positions = 0.0
    count = 0
    for pos in range(start_pos, end_pos):
        if pos in contig_positions[contig]:
            mean = position_ctype[contig, pos]
            mean_positions += mean
            count += 1
    if count == 0:
        return False
    if mean_positions/count >= cut_off:
        return True
    else:
        return False 

def alternate_bam(bamfile, per_base_ratio, outfile):
    contig_id = set()
    with open(per_base_ratio, 'r') as contig_file:
        for lines in contig_file:
            if lines.startswith("contig"):
                continue
            parts = lines.strip().split()
            contig = parts[0]
            contig_id.add(contig)
    with pysam.AlignmentFile(bamfile, "rb") as bam, pysam.AlignmentFile(outfile, 'wb', template=bam) as alt_bam:
        for read in bam.fetch(until_eof=True):
            if read.reference_name not in contig_id:
                alt_bam.write(read)

def main(args):
    print(f"{datetime.datetime.now()} INFO: bam_file_manipulate has started\n")
    start_time = time.time()

    cut_type = args.cut_type
    cut_off = args.cut_off
    if args.cut_off is None:
        if args.cut_type == "mean":
            args.cut_off = 1
        elif args.cut_type == "median":
            args.cut_off = 1

    per_base_ratio = args.per_base_ratio
    if not os.path.exists(per_base_ratio):
        print(f"{per_base_ratio} doesn't exist. Program exiting...")
        sys.exit(1)

    bamfile = args.bam_file
    if not os.path.exists(bamfile):
        print(f"{bamfile} doesn't exist. Program exiting...")
        sys.exit(1) 
    
    index_present = f"{bamfile}" + ".bai"
    index_file_created = False

    if os.path.exists(index_present):
        print(f'Index file:{index_present} is present\n')
    else:
        try:
            index_file_created = True
            print(f"Index file is not present for {bamfile}\n") 
            print(f"Creating index file for {bamfile}\n") 
            pysam.index(bamfile)
            print(f"Index file is created for {bamfile} in {index_present}\n")
        except Exception as e:
            print(f"Error while creating index file: {e}")
            sys.exit(1)

    name = args.name

    old_mapq = args.old_mapq
    if old_mapq != 30:
        print(f"User input MAPQ of {old_mapq} cut-off will be used\n")
    else:
        print(f"Default MAPQ of {old_mapq} will be used")

    new_mapq = args.new_mapq
    if new_mapq != 40:
        print(f"User input MAPQ of {new_mapq} will be assigned\n")
    else:
        print(f"Default MAPQ of {new_mapq} will be used")

    alternate_bamfile = True if args.alt_bam else False
    output_bamfile = f"{name}_manipulated.bam"

    print(f"Processing BAM file: {bamfile} with per base ratio file: {per_base_ratio} with cut off value: {cut_off} and cut off method: {cut_type} and with MAPQ cut-off: {old_mapq} and manipulated MAPQ: {new_mapq}\n")
    parsing_bamfiles(bamfile, per_base_ratio, cut_off, cut_type, old_mapq, new_mapq, output_bamfile)
    print(f"Output BAM file created: {output_bamfile}")

    if alternate_bamfile:
        alt_bam= f"{name}_alternate.bam"
        print(f"Alternate bam file without the targeted contigs will be output as {alt_bam}")
        print("Writing alternate bam file....")
        alternate_bam(bamfile, per_base_ratio, alt_bam)
        print(f"{alt_bam} writing finished!")
    else:
        print("No alternate bamfile will be produced")
    
    if index_file_created:
        print(f"Removing {index_present} file\n")
        os.remove(index_present)
        print(f"{index_present} file removed\n")

    print(f"bam_file_manipulate ran successfully!\n")

    print(f"Processing time: {time.time() - start_time:.2f} seconds")