#!/usr/bin/env python3

"""This script is the main core of the tool.
It uses a sequence file and a bed file with locus information
to create a user definded local kmer ratio file based on a global 
kmer file. It then calculates a per base kmer ratio
and outputs a file with per base statistics including the
mean, median, mode, max, min and poisson distribution
for each position in the contig-locus context that the user defined.

It requires the following input files:
1. A sequence file in FASTA or FASTQ format (can be gzipped)
2. A locus file in BED format with contig ID, start position, end position, and locus name
3. A global kmer file with kmer counts
4. name for the program that would also be in the output files
"""

import sys
import gzip
import os
from collections import defaultdict, Counter
import numpy as np 
import statistics
from Bio import SeqIO  
import time 
import datetime
import subprocess
import pyjellyfish as pyjf 
import importlib.resources as resources
import tempfile
import shutil

# This function reads a sequence file in FASTA or FASTQ format
# can be gzipped or not and return a dictionary with contig IDs as keys
# using Biopython's SeqIO module
def read_sequence_file(seq_file):
    sequences = {}
    gzipped = True if (seq_file.endswith(".gz")) else False
    if gzipped:
        if seq_file.endswith (".fasta.gz") or seq_file.endswith ("fa.gz"): 
            format = "fasta" 
        elif seq_file.endswith (".fastq.gz") or seq_file.endswith ("fq.gz"):
            format = "fastq"
        else:
            print(f"Sequence file not formatted correctly.\nProgram exiting...")
            sys.exit(1)
    else:
        if seq_file.endswith (".fasta") or seq_file.endswith ("fa"): 
            format = "fasta" 
        elif seq_file.endswith (".fastq") or seq_file.endswith ("fq"):
            format = "fastq" 
        else:
            print(f"Sequence file not formatted correctly.\nProgram exiting...")
            sys.exit(1)   
    open_func = gzip.open if gzipped else open
    with open_func(seq_file, "rt") as sfile:
        for record in SeqIO.parse(sfile, format):
            sequences[record.id] = str(record.seq).upper() #ensuring all sequences are in uppercase
    print(f"Read {len(sequences)} contigs from {seq_file}\n")
    return sequences

# This function is to produce the reverse complement of a sequence
# it can handle N bases as well
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def kmer_ratio_txt(genome_sequence, locus_file, global_kmer_file, canonical, kmer_ratio_outfile):
    valid_bases = {"A", "T", "G", "C"}
    print(f"\nCalculating kmer ratio from {locus_file} and {global_kmer_file}\n")
    sequence = read_sequence_file(genome_sequence)
    kmer_length = 0
    g_kmer = {}
    func_open = gzip.open if global_kmer_file.endswith(".gz") else open
    with func_open(global_kmer_file, "rt") as gkmer:
        for lines in gkmer:
            if lines.startswith('kmer'):
                continue
            kmer, count = lines.strip().split()
            g_kmer[kmer] = int(count)
            if kmer_length == 0:
                kmer_length = len(kmer)

    print(f"\nLoaded {len(g_kmer)} global kmers from {global_kmer_file}\n")
    print(f"k-mer length is: {kmer_length}")
    with open(locus_file, 'r') as locus:
        locus_dict = defaultdict(list)
        contig_dict = defaultdict(list)
        for line in locus:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) != 4:  # Checking for the locus file format
                print(f"\nLocus file: {locus_file} is not correctly formated\nIt should have contig/chromosome/scaffold\tstart\tend\tlocus_name\n")
                sys.exit(1)
            header_check = True  # Checking for headers
            try:
                int(parts[1])
                header_check = False
            except ValueError:
                header_check = True
            if header_check:
                continue
            start_pos = int(parts[1]) - 1 # To ensure python's 0-based system
            end_pos = int(parts[2])
            locus_name = parts[3]
            contig_id = parts[0] + "#" + str(start_pos)
            contig_dict[contig_id].append(start_pos)
            contig_dict[contig_id].append(end_pos)
            locus_dict[locus_name].append(contig_id)
    tmp_count = 0
    tmp_file_list = []
    for loc_name in locus_dict:
        lkmer_dict = defaultdict(Counter)
        locus_contig = defaultdict(set)
        tmp_count += 1
        temp_file_create = tempfile.NamedTemporaryFile(mode='w', delete=False) # f"tmp_lkmer.{loc_name}.txt"
        tmp_locus_kmer_file = temp_file_create.name
        tmp_file_list.append(tmp_locus_kmer_file)
        for contigs in locus_dict[loc_name]:
            contig_name = contigs.split("#")[0]
            start_pos = contig_dict[contigs][0]
            end_pos = contig_dict[contigs][1]
            seq = sequence[contig_name]
            seq_parse = seq[start_pos : end_pos]
            for i in range(len(seq_parse) - kmer_length + 1):
                kmer = seq_parse[i: i + kmer_length]
                if set(kmer) <= valid_bases:
                    if canonical == 1:
                        rev_kmer = min(kmer, reverse_complement(kmer))
                        lkmer_dict[loc_name][rev_kmer] += 1
                    else:
                        lkmer_dict[loc_name][kmer] += 1
            locus_contig[loc_name].add(contigs)
        with open(tmp_locus_kmer_file, "a") as tmp_lkmer:
            for locus_names in lkmer_dict:
                for contig_pos in locus_contig[locus_names]:
                    contig_name = contig_pos.split("#")[0]
                    start_pos = contig_dict[contig_pos][0]
                    end_pos = contig_dict[contig_pos][1]
                    for kmers, counts in lkmer_dict[locus_names].items():
                        g_kmer_count = g_kmer.get(kmers, 0) # global_kmer_query.query(kmers)
                        if g_kmer_count > 0:
                            global_value = g_kmer_count
                            ratio = counts / global_value if global_value > 0 else 0
                            tmp_lkmer.write(f"{contig_name}\t{locus_names}\t{start_pos + 1}\t{end_pos}\t{kmers}\t{counts}\t{ratio}\n")
        print(f"\nKmer ratios calculated and written to {tmp_locus_kmer_file}\n")
    with open(kmer_ratio_outfile, "w") as out_file:
        out_file.write(f"contig_id\tlocus_name\tstart\tend\tkmer\tcount\tratio\n")
        for temp_files in tmp_file_list:
            with open(temp_files, "r") as tmp_f:
                contents = tmp_f.read()
                out_file.write(contents)
    print(f"Removing all the local k-mer temp files\n")
    for temp_files in tmp_file_list:
        os.remove(temp_files)
    print(f"Removed local k-mer temp files\n")

# This function takes the genome sequence, user defined bed file with locus information,
# and a global kmer file output from jellyfish to calculate the kmer ratio for each contig-locus pair
# and outputs a file with the kmer ratio for each contig-locus pair
def kmer_ratio_jf(genome_sequence, locus_file, global_kmer_file, canonical, kmer_size, threads, kmer_ratio_outfile):
    print(f"\nCalculating kmer ratio from {locus_file} and {global_kmer_file}\n")
    sequence = read_sequence_file(genome_sequence)
    kmer_length = kmer_size
    if global_kmer_file.endswith(".jf"):
        print(f"\nFound global k-mers file made using jellyfish\n")
    print(f"k-mer length is: {kmer_length}")
    temp_dir = tempfile.TemporaryDirectory()
    temp_dir_name = temp_dir.name
    with open(locus_file, 'r') as locus:
        locus_dict = defaultdict(list)
        contig_dict = defaultdict(list)
        for line in locus:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) != 4:  # Checking for the locus file format
                print(f"\nLocus file: {locus_file} is not correctly formated\nIt should have contig/chromosome/scaffold\tstart\tend\tlocus_name\n")
                sys.exit(1)
            header_check = True  # Checking for headers
            try:
                int(parts[1])
                header_check = False
            except ValueError:
                header_check = True
            if header_check:
                continue
            start_pos = int(parts[1]) - 1 # To ensure python's 0-based system
            end_pos = int(parts[2])
            locus_name = parts[3]
            contig_id = parts[0] + "#" + str(start_pos)
            contig_dict[contig_id].append(start_pos)
            contig_dict[contig_id].append(end_pos)
            locus_dict[locus_name].append(contig_id)
    tmp_count = 0
    tmp_file_list = []
    global_kmer_query = pyjf.Jellyfish(global_kmer_file)
    for loc_name in locus_dict:
        lkmer_dict = defaultdict(Counter)
        locus_contig = defaultdict(set)
        tmp_count += 1
        temp_file_create = tempfile.NamedTemporaryFile(mode='w', delete=False) # f"tmp_lkmer.{loc_name}.txt"
        tmp_locus_kmer_file = temp_file_create.name
        tmp_file_list.append(tmp_locus_kmer_file)
        for contigs in locus_dict[loc_name]:
            contig_name = contigs.split("#")[0]
            start_pos = contig_dict[contigs][0]
            end_pos = contig_dict[contigs][1]
            seq_size = int(np.ceil((end_pos - start_pos)/1_000_000)) if (int(np.ceil((end_pos - start_pos)/1_000_000)) > 1) else 1
            seq = sequence[contig_name]
            seq_parse = seq[start_pos : end_pos]
            get_locus_kmer = get_parsed_kmer_jf(seq_parse, contig_name, kmer_length, seq_size, canonical, threads, temp_dir_name)
            if os.path.isfile(get_locus_kmer) and os.path.getsize(get_locus_kmer) > 0:
                print(f"Locus based k-mer file found\n")
            else:
                print(f"Locus based k-mer not found\n")
                sys.exit(1)
            with open(get_locus_kmer, 'r') as tmp_kmer:
                for lines in tmp_kmer:
                    parts = lines.strip().split()
                    kmer = parts[0]
                    counts = int(parts[1])
                    lkmer_dict[loc_name][kmer] += counts
            locus_contig[loc_name].add(contigs)
        with open(tmp_locus_kmer_file, "a") as tmp_lkmer:
            for locus_names in lkmer_dict:
                for contig_pos in locus_contig[locus_names]:
                    contig_name = contig_pos.split("#")[0]
                    start_pos = contig_dict[contig_pos][0]
                    end_pos = contig_dict[contig_pos][1]
                    for kmers, counts in lkmer_dict[locus_names].items():
                        g_kmer_count = global_kmer_query.query(kmers)
                        if g_kmer_count > 0:
                            global_value = g_kmer_count
                            ratio = counts / global_value if global_value > 0 else 0
                            tmp_lkmer.write(f"{contig_name}\t{locus_names}\t{start_pos + 1}\t{end_pos}\t{kmers}\t{counts}\t{ratio}\n")
        print(f"\nKmer ratios calculated and written to {tmp_locus_kmer_file}\n")
    with open(kmer_ratio_outfile, "w") as out_file:
        out_file.write(f"contig_id\tlocus_name\tstart\tend\tkmer\tcount\tratio\n")
        for temp_files in tmp_file_list:
            with open(temp_files, "r") as tmp_f:
                contents = tmp_f.read()
                out_file.write(contents)
    print(f"Removing all the local k-mer temp files\n")
    for temp_files in tmp_file_list:
        os.remove(temp_files)
    print(f"Removed local k-mer temp files\n")
    shutil.rmtree(temp_dir_name)

# This is the subfunction for kmer_ratio_jf to count local k-mer counts for the locus-based sequences
def get_parsed_kmer_jf(sequence, contig_id, kmer_size, seq_size, canonical, threads, temp_dir_name):
    temp_f_fasta = tempfile.NamedTemporaryFile(suffix=".fa", mode='w', delete=False)
    temp_fasta_path = os.path.join(temp_dir_name, temp_f_fasta.name)
    temp_seq = temp_fasta_path
    temp_f_jf = tempfile.NamedTemporaryFile(suffix=".jf", mode='w', delete=False)
    output_file = os.path.join(temp_dir_name, temp_f_jf.name) # f"tmp_parsed_kmer.jf"
    kmer_length = str(kmer_size)
    hash_size_command = f"{seq_size}M"
    thread_command = str(threads)
    out_command = f"{output_file}"  
    canonical_command = "False"
    with open(temp_seq, "w") as temp_fa:
        temp_fa.write(f">temp_seq\n")
        temp_fa.write(f"{sequence}\n")
    seq_command = temp_seq
    jellyfish_command = ["jellyfish", "count", "-m", f"{kmer_length}", "-t", f"{thread_command}", "-s", f"{hash_size_command}",  f"{seq_command}", "-o", f"{out_command}"]
    if canonical == 1:
        canonical_command = "True"
        jellyfish_command.append("-C")
    try:
        print(f"Parsed k-mers counting for {contig_id} started with {len(sequence)} parsed sequence\n")
        print(f"Options are:\nCanonical:{canonical_command}\nk-mer size: {kmer_length}\nhash table size: {hash_size_command}, threads: {thread_command}\n")
        print(f"Temporary file generated: {output_file}\n")
        subprocess.run(jellyfish_command, check= True)
    except subprocess.CalledProcessError as e:
        print(f"Error occured while calling jellyfish: {e}\n")
    
    print(f"Finished counting parsed k-mer for {contig_id}\nStarting jellyfish dump command\n")
    temp_parsed_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
    temp_parsed_filename = temp_parsed_file.name
    jellyfish_dump_command = ["jellyfish", "dump", "-c", f"{out_command}", "-o", f"{temp_parsed_filename}"]
    subprocess.run(jellyfish_dump_command, check=True)
    print(f"Temporary parsed k-mers dumped into {out_command} file\n")
    os.remove(temp_fasta_path)
    os.remove(output_file)
    return temp_parsed_filename

# This function calculates the per base kmer ratio
# for each contig-locus pair based on the kmer ratio file
# and the genome sequence file
# It outputs a file with per base statistics
# including mean, median, mode, max, min and poisson distribution
def per_base_ratio(sequence_file, lkmer_ratio_file, jf, global_kmer_file, canonical, repeat_ratio, per_base_outfile):
    print(f"\nCalculating per base kmer ratio from {lkmer_ratio_file} and {sequence_file}\n")
    genome = read_sequence_file(sequence_file)
    contig_kmer_ratio = defaultdict(list)
    contig_positions = defaultdict(list)
    contig_count = defaultdict(list)
    if jf:
        global_kmer_query = pyjf.Jellyfish(global_kmer_file)
    else:
        with gzip.open(global_kmer_file, 'rt') as gkmer_file:
            gkmer = {}
            for lines in gkmer_file:
                if lines.startswith("kmer"):
                    continue
                parts = lines.strip().split()
                kmer = parts[0]
                counts = parts[1]
                if kmer in gkmer:
                    gkmer[kmer] += int(counts)
                else:
                    gkmer[kmer] = int(counts)
    kmer_length = None
    with open(lkmer_ratio_file, "r") as lratio:
        for line in lratio:
            if line.startswith("contig_id"):
                continue
            contig, locus, start_pos, end_pos, kmer, count, ratio = line.strip().split()
            start_pos = int(start_pos) + 1
            contig_pos = contig + "#" + str(start_pos)
            if kmer_length is None:
                kmer_length = len(kmer)
            contig_kmer_ratio[(contig_pos, locus)].append((kmer, float(ratio)))
            contig_positions[(contig_pos, locus)].append((int(start_pos), int(end_pos)))
            contig_count[(contig_pos, locus)].append((kmer))
    print(f"\nLoaded kmer ratios for {len(contig_kmer_ratio)} contigs from {lkmer_ratio_file}\n")
    with open(per_base_outfile, "w") as ofile:
        ofile.write("contig\tposition\tratio\tmean\tmedian\tmode\tmax\tmin\trepetitive\tlunique\n")
        for (contig_pos, locus), kmer_ratios in contig_kmer_ratio.items():
            for start_pos, end_pos in contig_positions[(contig_pos, locus)]:
                start_pos = int(start_pos)
                end_pos = int(end_pos)
            contig_name = contig_pos.split("#")[0]
            print(f"\nStart coordinate: {start_pos + 1}, End coordinate: {end_pos}\n")
            if contig_name in genome:
                sequence = genome[contig_name][start_pos:end_pos]
                kmer_pos = defaultdict(list)
                for i in range(len(sequence) - kmer_length + 1):
                    kmer = sequence[i:i + kmer_length]
                    if canonical == 1:
                        kmer_rc = reverse_complement(kmer)
                        kmer = min(kmer, kmer_rc)
                    kmer_pos[kmer].append(i + start_pos)
                pos_ratio = defaultdict(float)
                pos_total_count = defaultdict(int)
                pos_ratio_counter = defaultdict(list)
                pos_counts = defaultdict(list)
                for values in kmer_ratios:
                    kmer, ratio = values
                    kmer = min(reverse_complement(kmer), kmer) if canonical == 1 else kmer
                    positions = kmer_pos.get(kmer, [])
                    for p in positions:
                        for i in range(kmer_length):
                            pos = p + i
                            pos_ratio[pos] = pos_ratio.get(pos, 0) + ratio
                            pos_total_count[pos] = pos_total_count.get(pos, 0) + 1
                            pos_ratio_counter[pos].append(ratio)
                print(f"\nCalculated per base kmer ratios for contig {contig_name} with {len(pos_ratio)} positions\n")
                for kmer in contig_count[(contig_pos, locus)]:
                    kmer = min(reverse_complement(kmer), kmer) if canonical == 1 else kmer
                    positions_count = kmer_pos.get(kmer, [])
                    if jf:
                        gkmer_count = global_kmer_query.query(kmer)
                    else:
                        gkmer_count = gkmer[kmer]
                    for p in positions_count:
                        for i in range(kmer_length):
                            pos = p + i
                            pos_counts[pos].append(gkmer_count)
                # Calculate statistics for each position
                for nt_pos, ratio in pos_ratio.items():
                    total_events = pos_total_count[nt_pos] if pos_total_count[nt_pos] > 0 else 1
                    total_count = np.sum(pos_counts[nt_pos])
                    count_qual = 1 if total_count > total_events else 0
                    mean_val = ratio / total_events
                    ratio_qual = 1 if mean_val >= repeat_ratio else 0
                    median_val = np.median(pos_ratio_counter[nt_pos])
                    max_val = np.max(pos_ratio_counter[nt_pos])
                    min_val = np.min(pos_ratio_counter[nt_pos])
                    mode_val = statistics.multimode(pos_ratio_counter[nt_pos])
                    if count_qual == 1 and ratio_qual == 1:
                        repetitive = 1
                        local_unique = 1
                    elif count_qual == 1 and ratio_qual != 1:
                        repetitive = 1
                        local_unique = 0
                    elif count_qual != 1 and ratio_qual == 1:
                        repetitive = 0
                        local_unique = 1
                    else:
                        repetitive = 0
                        local_unique = 0
                    ofile.write(f"{contig_name}\t{nt_pos + 1}\t{ratio}\t{mean_val}\t{median_val}\t{",".join(str(item) for item in mode_val)}\t{max_val}\t{min_val}\t{repetitive}\t{local_unique}\n")
    print(f"\nPer base kmer ratios written to {per_base_outfile}\n")

def plotting_data(per_base_file, swindow, name_plot):
    with resources.path("kmerRRR.R", "pbr_plot.r") as r_script:
        subprocess.run(["Rscript", str(r_script), str(per_base_file), str(swindow), str(name_plot)], check=True)

def main(args):
    start_time = time.time()
    print(f"{datetime.datetime.now()} INFO: summary_stat has started\n")
    sequence_file = args.sequence_file
    if not os.path.exists(sequence_file):
        print(f"Error: Sequence file {sequence_file} does not exist\n")
        sys.exit(1)

    locus_file = args.bed_file
    if not os.path.exists(locus_file):
        print(f"Error: Locus file {locus_file} does not exist\n")
        sys.exit(1)

    global_kmer_file = args.global_kmer_file
    if not os.path.exists(global_kmer_file):
        print(f"Error: Global kmer file {global_kmer_file} does not exist\n")
        sys.exit(1)
    
    canonical = args.canonical
    repeat_ratio = args.repeat_ratio
    name = args.name
    kmer_size = args.kmer_size

    kmer_ratio_outfile = f"{name}_kmer_ratio.txt"
    per_base_outfile = f"{name}_per_base_ratio.txt"
    print(f"Starting kmer ratio calculation for {sequence_file} with locus {locus_file}, global kmer file {global_kmer_file} and canonical: {canonical}\n")
    if global_kmer_file.endswith(".txt.gz"):
        print(f"Text file for global kmers found {global_kmer_file}.txt.gz\n")
        print(f"Kmer summary statistics will be calculated from {global_kmer_file}.txt.gz\n")
        kmer_ratio_txt(sequence_file, locus_file, global_kmer_file, canonical, kmer_ratio_outfile)
    elif global_kmer_file.endswith(".jf"):
        print(f"{global_kmer_file} is a jellyfish data structure\n")
        threads = args.threads
        print(f"Kmers summary statistics will be calculated from {global_kmer_file} jellyfish file\n")
        kmer_ratio_jf(sequence_file, locus_file, global_kmer_file, canonical,kmer_size, threads, kmer_ratio_outfile)
    else:
        print(f"\n{global_kmer_file} is in wrong format. System exiting...\n")
        sys.exit(1)

    print(f"Kmer ratio calculation completed. Now calculating per base ratio.\n")
    print(f"Repeat parameter to be used: {repeat_ratio}\n")
    jf = True if global_kmer_file.endswith(".jf") else False
    per_base_ratio(sequence_file, kmer_ratio_outfile, jf, global_kmer_file, canonical, repeat_ratio, per_base_outfile)
    print(f"Per base ratio calculation completed. Output written to {per_base_outfile}\n")
    plot = args.plot
    if plot:
        slide_window = args.slide_window
        print(f"\nGenerating plots for mean and median\n")
        plotting_data(per_base_outfile, slide_window, name)
    else:
        print("No plots will be generated\n")
        
    print(f"Summary stat run completed successfully!\n")
    print(f"Processing time: {time.time() - start_time:.2f} seconds")
