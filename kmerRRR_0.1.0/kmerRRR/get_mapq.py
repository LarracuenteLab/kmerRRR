#!/usr/bin/env python3

import pysam
import sys
import os
from collections import defaultdict
import time
import datetime
import subprocess
import importlib.resources as resources

def get_mapq(bamfile, contig_bed_file, output_file):
    contig = []
    with open(contig_bed_file, 'r') as bedfile:
        for line in bedfile:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) != 4:  # Checking for the locus file format
                print(f"\nLocus file: {bedfile} is not correctly formated\nIt should have contig/chromosome/scaffold\tstart\tend\tlocus_name\n")
                sys.exit(1)
            header_check = True  # Checking for headers
            try:
                int(parts[1])
                header_check = False
            except ValueError:
                header_check = True
            if header_check:
                continue
            parts = line.strip().split()
            contig_id = parts[0]
            contig.append(contig_id)
    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        with open(output_file, 'w') as out:
                out.write("contig\tstart\tend\tMapQ\tLUC\n")
                for contig_id in contig:
                     count = 0
                     skip = 0
                     print(f"Fetching MapQ for contig: {contig_id}")
                     try: #To stop the program from exiting before running completely
                        for read in bam.fetch(contig_id):
                            if read.is_unmapped:
                                skip += 1
                                continue
                            count += 1
                            if read.has_tag("KR"):
                                tag_bool = 1
                            else:
                                tag_bool = 0
                            out.write(f"{contig_id}\t{int(read.reference_start) + 1}\t{int(read.reference_end) + 1}\t{read.mapping_quality}\t{tag_bool}\n")
                        print(f"\nFor {contig_id}, MAPQ extracted for {count} reads, MAPQ skipped for {skip} reads\n")
                     except ValueError as e:
                         print(f"Contig: {contig_id} was not found. Error: {e}\n")
                         pass

def get_mapq_locus(bamfile, contig_bed_file, output_file):
    contig = defaultdict(list)
    with open(contig_bed_file, 'r') as bedfile:
        for line in bedfile:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) != 4:  # Checking for the locus file format
                print(f"\nLocus file: {bedfile} is not correctly formated\nIt should have contig/chromosome/scaffold\tstart\tend\tlocus_name\n")
                sys.exit(1)
            header_check = True  # Checking for headers
            try:
                int(parts[1])
                header_check = False
            except ValueError:
                header_check = True
            if header_check:
                continue
            parts = line.strip().split()
            contig_id = parts[0]
            start_pos = int(parts[1]) - 1 # pysam's and python's 0-based system
            end_pos = int(parts[2]) - 1 # pysam's and python's 0-based system
            contig_unique = contig_id + "#" + str(start_pos)
            if contig_unique not in contig:
                contig[contig_unique].append((start_pos, end_pos))
            else:
                print(f"\nContig: {contig_unique.split("#")[0]} has duplicte entry with the same start position\nSystem exiting...")
                sys.exit(1)
    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        with open(output_file, 'w') as out_f:
            out_f.write("contig\tstart\tend\tMapQ\tLUC\n")
            for contigs, positions in contig.items():
                count = 0
                skip = 0
                contig_name = contigs.split("#")[0]
                start, end = positions[0][0], positions[0][1]
                print(f"\nFetching MapQ for contig: {contig_name}, start position: {start}, end position: {end}\n")
                try:
                    for read in bam.fetch(contig_name):
                        if read.is_unmapped:
                            skip += 1
                            continue
                        if read.reference_start >= start and read.reference_start <= end:
                            count += 1
                            if read.has_tag("KR"):
                                tag_bool = 1
                            else:
                                tag_bool = 0
                            out_f.write(f"{contig_name}\t{int(read.reference_start) + 1}\t{int(read.reference_end) + 1}\t{read.mapping_quality}\t{tag_bool}\n")
                except ValueError as e:
                    print(f"Contig: {contig_name} was not found. Error: {e}\n")
                    pass
                print(f"\nFor {contig_name}, MAPQ extracted for {count} reads and skipped for {skip} reads\n")


def plotting_data(mapq_file, slide_window, name_plot):
    with resources.path("kmerRRR.R", "mapq_plot.r") as r_script:
        subprocess.run(["Rscript", str(r_script), str(mapq_file), str(slide_window), str(name_plot)], check=True)  

def main(args):
    print(f"{datetime.datetime.now()} INFO: get_mapq has started\n")
    start_time = time.time()

    bamfile = args.bam_file
    if not os.path.exists(bamfile):
        print(f"{bamfile} doesn't exist. Program exiting...")
        sys.exit(1)

    bedfile = args.bed_file
    if not os.path.exists(bedfile):
        print(f"{bedfile} doesn't exist. Program exiting...")
        sys.exit(1)

    name = args.name
    locus = args.locus_based
    output_file = f"{name}_MapQ.txt"

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

    print(f"Bam file: {bamfile}\nBed file: {bedfile}\nName of the program: {name}\nOutput file to be generated: {output_file}\n")
    if locus:
        get_mapq_locus(bamfile, bedfile, output_file)
        print(f"Locus based MAPQL scores are being extracted...\n")
    else:
        get_mapq(bamfile, bedfile, output_file)
        print(f"MAPQ scores are being extracted for the full contig (not locus based)...")
    print(f"MAPQ scores have been extracted from {bamfile} and output into {output_file}\n")

    plot = args.plot
    sliding_window = args.slide_window
    if plot:
        print(f"\nGenerating plots for MAPQ score\n")
        plotting_data(output_file, sliding_window, name)
    else:
        print(f"\nNo plots will be generated\n")
    
    if index_file_created:
        print(f"Removing {index_present} file\n")
        os.remove(index_present)
        print(f"{index_present} file removed\n")

    print(f"Get mapq program ran successfully!\n")

    print(f"Processing time: {time.time() - start_time:.2f} seconds")