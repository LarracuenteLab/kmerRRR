#!/usr/bin/env python3

import os
import gzip
import time
import datetime
import subprocess

# This function is to get the summary stats on the global kmers list
# this is also a stand alone program that the user can use
# to report the summary in std out or in a file.
def gkmer_info_txt(global_kmer_file):
    tmp_dict = {}
    with gzip.open(global_kmer_file, 'rt') as gkmer:
        for lines in gkmer:
            if lines.startswith("kmer"):
                continue
            kmer, count = lines.strip().split()
            tmp_dict[kmer] = int(count)
    singletons =0 
    total_count = 0
    distinct_kmers = len(tmp_dict)
    max_count = 0
    max_kmer = ""
    for key, values in tmp_dict.items(): 
        total_count += values
        if values == 1:
            singletons +=1
        if values > max_count:
            max_count = values
            max_kmer = key
    print(f"Total kmers count: {total_count}\nTotal number of distinct kmers: {distinct_kmers}\nNumber of singletons: {singletons}\nMax kmer:{max_kmer} with max count: {max_count}\n")


# This function is to get the summary stats on the global kmers list
# this is also a stand alone program that the user can use
# to report the summary in std out or in a file
def gkmer_info_jf(global_kmer_file):
    if global_kmer_file.endswith(".jf"):
        print(f"\n {global_kmer_file} is a jellyfish output file")
    else:
        print(f"\nGlobal k-mer ({global_kmer_file}) is not output using jellyfish. System exiting...\n")
    jellyfish_stats_command = ["jellyfish", "stats", f"{global_kmer_file}"]
    print(f"\nPrinting the basic statistics of global k-mers from {global_kmer_file}\n")
    output = subprocess.run(jellyfish_stats_command, capture_output=True, text= True).stdout
    print(f"\n{output}\n")


def main(args):
    start_time = time.time()
    print(f"{datetime.datetime.now()} Info: gkmer_info started")
    global_kmer_file = args.global_kmer_file
    if os.path.exists(global_kmer_file):
        if global_kmer_file.endswith(".jf"):
            print(f"{global_kmer_file} is made with jellyfish\n")
            print(f"Global kmers info is extracting from {global_kmer_file} jellyfish k-mer count file\n")
            gkmer_info_jf(global_kmer_file)
        elif global_kmer_file.endswith(".txt.gz"):
            print(f"Text file for global kmers found {global_kmer_file}.txt.gz\n")
            print(f"Extracting global kmers info from {global_kmer_file}.txt.gz\n")
            gkmer_info_txt(global_kmer_file)
    print(f"gkmer_info ran successfully!\n")
    print(f"Processing time: {time.time() - start_time:.2f} seconds")