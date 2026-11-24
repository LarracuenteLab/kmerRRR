#!/usr/bin/env python3

"""This code is to calculate the mean of LUCs from one or more 
files to test what value of k is suitable for the repeat structure.
This script takes a list of per base ratio files or a single 
per base ratio file and a BED locus file containing the locus of 
interest and name or path of an output file. It outputs a file with
the mean for all LUCs for that locus.
ONE CAVEAT, make sure the file name has only one number that corresponds
to the size k of k-mers."""

import sys
import argparse
from collections import defaultdict
import datetime
import time
import re, os
import matplotlib.pyplot as plt

def kmean(file_list, bed_file, name):
    new_name = f"{name}_LUCmean.txt"
    locus_dict = defaultdict(list)
    with open(bed_file, 'r') as bedf:
        for lines in bedf:
            parts = lines.strip().split()
            if len(parts) != 4:  # Checking for the locus file format
                print(f"\nLocus file: {bed_file} is not correctly formated\nIt should have contig/chromosome/scaffold\tstart\tend\tlocus_name\n")
                sys.exit(1)
            header_check = True  # Checking for headers
            try:
                int(parts[1])
                header_check = False
            except ValueError:
                header_check = True
            if header_check:
                continue
            contig_name = parts[0]
            start = 0 if int(parts[1]) == 0 else int(parts[1]) - 1 
            end = int(parts[2])
            locus = parts[3]
            contig_locus = f"{contig_name}#{locus}"
            locus_dict[contig_locus].extend([start, end])
    print(f"{datetime.datetime.now()} INFO: {len(locus_dict)} locus read from {os.path.basename(bed_file)}\n")
    with open(new_name, 'w') as out_f:
        out_f.write(f"contig\tlocus\tkmer\tstart\tend\tlu_cumulative_score\tlu_mean\trep_mean\tuni_mean\tuncertain_counts\n")
        for files in file_list:
            file_name = os.path.basename(files)
            match_digit = re.search(r"\d+", file_name)
            if match_digit:
                kmer = match_digit.group()
                print(f"{datetime.datetime.now()} INFO: For {files}, size of k is {kmer}\n")
            else:
                print(f"{datetime.datetime.now()} INFO: Error in file name, no numbers of k found\n")
            for key, items in locus_dict.items():
                locus_unique = 0
                repetitive = 0
                unique = 0
                count = 0
                uncertain = 0
                cname = key.split("#")[0]
                lname = key.split("#")[1]
                start_pos = items[0]
                end_pos = items[1]
                locus_size = end_pos - start_pos
                with open(files, 'r') as in_f:
                    for lines in in_f:
                        parts = lines.strip().split()
                        header_check = True  # Checking for headers
                        try:
                            int(parts[1])
                            header_check = False
                        except ValueError:
                            header_check = True
                        if header_check:
                            continue
                        if parts[0] == cname:
                            if int(parts[1]) >= start_pos or int(parts[1]) <= end_pos :
                                if int(parts[8]) == 1 and int(parts[9]) == 1:
                                    locus_unique +=1
                                    count += 1
                                elif int(parts[8]) == 0 and int(parts[9]) == 1:
                                    unique += 1
                                    count += 1
                                elif int(parts[8]) == 1 and int(parts[9]) == 0:
                                    repetitive += 1
                                    count += 1
                                else:
                                    uncertain += 1
                                    count += 1
                if count == locus_size:
                    print(f"{datetime.datetime.now()} INFO: The locus size and the total count is same\n")
                else:
                    if count < locus_size:
                        print(f"{datetime.datetime.now()} INFO: The locus size parsed is less than the size from BED locus file by {locus_size - count}\n")
                    else:
                        print(f"{datetime.datetime.now()} INFO: The locus size parsed is greater than the size from BED locus file by {count - locus_size}\n")
                print(f"{datetime.datetime.now()} INFO: Size of {count} is used to calculate mean\n")
                lunique_mean = locus_unique/count
                repeat_mean = repetitive/count
                unique_mean = unique/count
                out_f.write(f"{cname}\t{lname}\t{kmer}\t{start_pos}\t{end_pos}\t{locus_unique}\t{lunique_mean}\t{repeat_mean}\t{unique_mean}\t{uncertain}\n")

def plotting_data(input_file, output_file):
    dir_name = os.path.dirname(output_file)
    out_name = os.path.basename(output_file)
    if dir_name:
        new_dir = f"{dir_name}/plotting_data"
    else:
        new_dir="plotting_data"
    os.makedirs(new_dir, exist_ok=True)
    print(f"{datetime.datetime.now()} INFO: Made {new_dir}\n")
    conloc_kmers = defaultdict(list)
    conloc_lunique = defaultdict(list)
    with open(input_file, 'r') as in_f:
        for lines in in_f:
            parts = lines.strip().split()
            header_check = True  # Checking for headers
            try:
                int(parts[2])
                header_check = False
            except ValueError:
                header_check = True
            if header_check:
                continue
            conloc = f"{parts[0]}#{parts[1]}"
            conloc_kmers[conloc].append(int(parts[2]))
            conloc_lunique[conloc].append(float(parts[6]))
    for keys, items in conloc_kmers.items():
        kmers = items
        lunique_mean = conloc_lunique[keys]
        plt.figure(figsize=(12,8))
        plt.plot(kmers, lunique_mean, color="darkblue", linestyle="--", marker='o')
        plt.title("LUC mean for different kmers")
        plt.xlabel("k size")
        plt.ylabel("LUC locus unique mean")
        plt.tight_layout()
        plt.savefig(os.path.join(new_dir, f"{out_name}.{keys}.pdf"), format="pdf", dpi = 300)
        plt.close()

def main(args, arg_parser):
    start_time = time.time()
    if args.name:
        output = args.name
    else:
        output = "output"
    if args.bed_file or args.files:
        bed_file = args.bed_file
        files_list = args.files
        print(f"{datetime.datetime.now()} INFO: Starting to calculate LUC means\n")
        kmean(files_list, bed_file, output)
    else:
        print(f"{datetime.datetime.now()} INFO: Not all arguments are found\n")
        arg_parser.print_help()
        sys.exit(1)

    if args.plot:
        print(f"{datetime.datetime.now()} INFO: Plot will be output\n")
        if args.name:
            plot_out = f"{args.name}_plot"
        else:
            plot_out=f"output_plot"
        if not args.plot_file:
            if os.path.exists(f"{output}_LUCmean.txt"):
                print(f"{datetime.datetime.now()} INFO: Input file: {output}_LUCmean.txt will be used\n")
                plotting_data(f"{output}_LUCmean.txt", plot_out)
            else:
                print(f"{datetime.datetime.now()} INFO: {output}_LUCmean.txt file not found\nSystem exiting...\n")
                sys.exit(1)
        else:
            plotting_data(args.plot_file, plot_out)
            print(f"{datetime.datetime.now()} INFO: Input file: {args.plot_file}\n")
    else:
        print(f"{datetime.datetime.now()} INFO: No plot will be output\n")
    print(f"{datetime.datetime.now()} INFO: Total time reqcuired: {time.time() - start_time:.2f} seconds")

if __name__=="__main__":
    arg_parser = argparse.ArgumentParser("LUCmean")
    arg_parser.add_argument("-bed", "--bed_file", help="Name or Path of the BED locus file")
    arg_parser.add_argument("-f", "--files", nargs="+", help="Name or path of the per base ratio files")
    arg_parser.add_argument("-n", "--name", help="Name of the outputfile, default: output_LUCmean.txt")
    arg_parser.add_argument("-p", "--plot", action='store_true', help="PLotting the data")
    arg_parser.add_argument("--plot_file", help="Name or path of the file with which plotting will be done")
    args = arg_parser.parse_args()
    main(args, arg_parser)