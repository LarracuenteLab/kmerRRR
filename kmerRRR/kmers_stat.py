#!/usr/bin/env python3

"""This script is the main core of the tool.
It uses a sequence file to extract all the global kmers
and a bed file with locus information to create a user-defined
local kmer ratio summary statitstics file based on a global kmer 
counts. It calculates a per-base kmer ratio and outputs a file
with per-base statistics including mean, median, mode, max, min, 
and poisson distribution for each position in the contig-locus 
context that the user defined.

It requires the following inputs from the users:
1. A sequence file in FASTA or FASTQ format (can be gzipped)
2. A locus file in BED format with contig ID, start position, end position, and locus name
3. kmer size based on the user input (integer value)
4. canonical value (0 or 1 for no and yes, respectively)
5. name for the program that would also be in the output files
"""

import sys
import gzip
import os
import shutil
import tempfile
from collections import defaultdict, Counter
import numpy as np  
import statistics
from Bio import SeqIO    
import io
import time
import datetime
import subprocess
import pyjellyfish as pyjf 
import importlib.resources as resources

# This function reads a sequence file in FASTA or FASTQ format.
# Can be gzipped or not and returns a dictionary with contig IDs as keys
# using Biopython's SeqIO module.
def read_sequence_file(seq_file):
    sequences = {}
    gzipped = True if (seq_file.endswith(".gz")) else False
    if gzipped:
        if seq_file.endswith(".fasta.gz") or seq_file.endswith("fa.gz"):
            format = "fasta"
        elif seq_file.endswith(".fastq.gz") or seq_file.endswith("fq.gz"):
            format = "fastq"
        else:
            print(f"\nSequence file not formatted correctly.\nProgram exiting...")
            sys.exit(1)
    else:
        if seq_file.endswith(".fasta") or seq_file.endswith("fa"):
            format = "fasta"
        elif seq_file.endswith(".fastq") or seq_file.endswith("fq"):
            format = "fastq"
        else:
            print(f"\nSequence file not formatted correctly.\nProgram exiting...")
            sys.exit(1)
    open_func = gzip.open if gzipped else open
    with open_func(seq_file, "rt") as sfile:
        for record in SeqIO.parse(sfile, format):
            sequences[record.id] = str(record.seq).upper()  # ensuring all sequences are in uppercase
    print(f"Read {len(sequences)} contigs from {seq_file}\n")
    return sequences


# This function is to produce the reverse complement of a sequence.
# It can handle N bases as well.
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# This function is to get all the local kmers
# which will be written in a temporary file.
def contig_kmers(genome_file, kmer_size, canonical, lkmer_out_file):
    valid_bases = {"A", "T", "G", "C"}
    print(f"\nExtracting local kmers from genome file {genome_file} with kmer size {kmer_size} and canonical {canonical}\n")
    genome = read_sequence_file(genome_file)
    # contig_dict = defaultdict(list)
    tmp_file_names = []
    i = 1
    for contig, seq in genome.items():
        kmer_dict = defaultdict(int)
        contig_temp_file_create = tempfile.NamedTemporaryFile(mode='wb', delete=False) # f"tmp_{i}.gz"
        tmp_file = contig_temp_file_create.name
        tmp_file_names.append(tmp_file)
        with gzip.GzipFile(tmp_file, "w") as tmp:
            with io.TextIOWrapper(tmp, encoding='utf-8') as tmp_file:
                for j in range(len(seq) - kmer_size + 1):
                    kmer = seq[j:j + kmer_size]
                    if set(kmer) <= valid_bases:
                        if canonical == 1:
                            kmer = min(kmer, reverse_complement(kmer))
                        kmer_dict[kmer] += 1
                    else:
                        continue
                for kmers, count in kmer_dict.items():
                    tmp_file.write(f"{contig}\t{kmers}\t{count}\n")
                print(f"\nExtracted {len(kmer_dict)} kmers for contig {contig}\n")
        i += 1
    
    Buffer_size = 1024*1024*1024
    print(f"Buffer size:{Buffer_size} bytes or {Buffer_size/1024**3} GB")
    with gzip.open(lkmer_out_file, 'wb') as f:
        for tmp_f in tmp_file_names:
            with gzip.open(tmp_f, 'rb') as in_f:
                shutil.copyfileobj(in_f, f, length= Buffer_size)
        print(f"\nAll local kmers written to {lkmer_out_file}\n")
        print(f"\nTotal temporary file created: {len(tmp_file_names)}\n")
    print(f"\nExtracted contig kmers for all the {len(genome)} contigs with kmer size {kmer_size}\n")
    return tmp_file_names

def global_kmers_txt(lkmers_file, global_kmer_file):
    gkmer_dict = defaultdict(int)
    with gzip.open(lkmers_file, 'rt') as lkmer:
        for lines in lkmer:
            if lines.startswith("contig"):
                continue
            contig, kmer, count = lines.strip().split()
            gkmer_dict[kmer] += int(count)
    buffer_size = 1_000_000
    buffer =[]
    with gzip.open(global_kmer_file, 'wb') as out_f:
        with io.TextIOWrapper(out_f, encoding= 'utf-8') as gkmer:
            gkmer.write(f"kmer\tcount\n")
            for kmers, counts in gkmer_dict.items():
                buffer.append(f"{kmers}\t{counts}\n")
                if len(buffer) >= buffer_size:
                    gkmer.writelines(buffer)
                    buffer.clear()
            if buffer:
                gkmer.writelines(buffer)

# This function is to get the summary stats on the global kmers list.
# This is also a stand-alone program that the user can use
# to get the summary in std out or in a file.
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
    distinc_kmers = len(tmp_dict)
    max_count = 0
    max_kmer = ""
    for key, values in tmp_dict.items(): 
        total_count += values
        if values == 1:
            singletons +=1
        if values > max_count:
            max_count = values
            max_kmer = key
    print(f"Total kmers count: {total_count}\nTotal number of distinct kmers: {distinc_kmers}\nNumber of singletons: {singletons}\nMax kmer:{max_kmer} with max count: {max_count}")

# Jellyfish implementation for global_kmers
def global_kmers_jf(sequence, canonical, kmer_size, genome_size, threads, output_file):
    kmer_length = int(kmer_size)
    hash_size_command = f"{genome_size}M"
    thread_command = int(threads)
    out_command = f"{output_file}"  
    canonical_command = "False"
    seq_command = sequence
    jellyfish_command = ["jellyfish", "count", "-m", f"{kmer_length}", "-t", f"{thread_command}", "-s", f"{hash_size_command}",  "-o", f"{out_command}"]
    if canonical == 1:
        canonical_command = "True"
        jellyfish_command.append("-C")
    if sequence.endswith(".gz"):
        jellyfish_command.append("/dev/stdin")
        zcat_process = subprocess.Popen(["zcat", f"{seq_command}"], stdout=subprocess.PIPE)
        try:
            #print(f"Global k-mers counting for {sequence} started\n")
            print(f"Options are:\nCanonical:{canonical_command}\nk-mer size: {kmer_length}\nhash table size: {hash_size_command}, threads: {thread_command}\n")
            print(f"Output file generated: {output_file}")
            subprocess.run(jellyfish_command, stdin=zcat_process.stdout, check=True)
            zcat_process.stdout.close()
            zcat_process.wait()
        except subprocess.CalledProcessError as e:
            print(f"Error occured while calling jellyfish: {e}")
    else:
        jellyfish_command.append(seq_command)
        try:
            #print(f"Global k-mers counting for {sequence} started\n")
            print(f"Options are:\nCanonical:{canonical_command}\nk-mer size: {kmer_length}\nhash table size: {hash_size_command}, threads: {thread_command}\n")
            print(f"Output file generated: {output_file}")
            subprocess.run(jellyfish_command, check= True)
        except subprocess.CalledProcessError as e:
            print(f"Error occured while calling jellyfish: {e}")

# To get global_kmer stats
def gkmer_info_jf(global_kmer_file):
    if global_kmer_file.endswith(".jf"):
        print(f"\n {global_kmer_file} is a jellyfish output file")
    else:
        print(f"\nGlobal k-mer ({global_kmer_file}) is not output using jellyfish. System exiting...\n")
    jellyfish_stats_command = ["jellyfish", "stats", f"{global_kmer_file}"]
    print(f"\nPrinting the basic statistics of global k-mers from {global_kmer_file}\n")
    output = subprocess.run(jellyfish_stats_command, capture_output=True, text= True).stdout
    print(f"\n{output}\n")

# This function takes the genome sequence, user defined bed file with locus information,
# and a global kmer file output from jellyfish to calculate the kmer ratio for each contig-locus pair
# and outputs a file with the kmer ratio for each contig-locus pair.
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

# This is the subfunction for kmer_ratio_jf to count local k-mer counts for the locus-based sequences.
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

# This function takes the genome sequence, user-defined bed file with locus information,
# and a global kmer file to calculate the kmer ratio for each contig-locus pair
# and outputs a file with the kmer ratio for each contig-locus pair.
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

# This function calculates the per-base kmer ratio
# for each contig-locus pair based on the kmer ratio file
# and the genome sequence file.
# It outputs a file with per base statistics
# including the mean, median, mode, max, min and poisson distribution.
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
            start_pos = int(start_pos) - 1
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
            print(f"\nStart coordinate: {start_pos}, End coordinate: {end_pos}\n")
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

# To plot the mean and median summary stat for per base ratio file.
def plotting_data(per_base_file, swindow, name_plot):
    with resources.path("kmerRRR.R", "pbr_plot.r") as r_script:
        subprocess.run(["Rscript", str(r_script), str(per_base_file), str(swindow), str(name_plot)], check=True)

def main(args):
    start_time = time.time()
    print(f"{datetime.datetime.now()} INFO: kmers_stat has started\n")

# user arguments processing and checking

    # Input files checking and processing
    sequence_file = args.sequence_file
    if not os.path.exists(sequence_file):
        print(f"\nError: Sequence file {sequence_file} does not exist.\n")
        sys.exit(1)
    locus_file = args.bed_file
    if not os.path.exists(locus_file):
        print(f"\nError: Locus file {locus_file} does not exist.\n")
        sys.exit(1)

    kmer_size = args.kmer_size
    if not kmer_size:
        print(f"\nKmer size input is required\n")
        sys.exit(1)

    canonical = args.canonical
    repeat_ratio = args.repeat_ratio
    name = args.name
    if not name:
        print(f"\nName for the output files is required\n")
        sys.exit(1)

    jellyfish = args.jellyfish
    genome_size = args.genome_size

    print(f"\nInput sequence file: {sequence_file}\nKmer size: {kmer_size}\nCanonical: {canonical}\nLocus file: {locus_file}\nOutput name: {name}\n")
    if jellyfish:
        print(f"\nJellyfish will be used to count global k-mers\n")
    else:
        print(f"\nPython's dictionary or hash algorithm will be used to count the global k-mers\n")

    # For extracting the kmers from the sequence
    if not jellyfish:
        lkmer_out_file = f"{name}.contig_kmers.txt.gz"
        print(f"\nTemporary local kmers file will be created as {lkmer_out_file}\n")

    if jellyfish:
        threads = args.threads
        if threads == None:
            print(f"Threads command was not input. System exiting...\n")
            sys.exit(1)
        global_kmer_outfile =f"{name}.jf"
        if os.path.exists(global_kmer_outfile):
            print(f"{global_kmer_outfile} already exists, removing the previous database...\n")
            os.remove(global_kmer_outfile)
            print(f"{global_kmer_outfile} REMOVED!")
        print(f"\nGlobal kmers will be written to {global_kmer_outfile}")
    else:
        global_kmer_outfile = f"{name}_global_kmers.txt.gz"
        if os.path.exists(global_kmer_outfile):
            print(f"{global_kmer_outfile} already exists, removing the previous database...\n")
            os.remove(global_kmer_outfile)
            print(f"{global_kmer_outfile} REMOVED!")
        print(f"\nGlobal kmers will be written to {global_kmer_outfile}\n")

    print("\nStarting kmer extraction...\n")
    
    if jellyfish:
        global_kmers_jf(sequence_file, canonical, kmer_size, genome_size, threads, global_kmer_outfile)
    else:
        temp_file_names = contig_kmers(sequence_file, kmer_size, canonical, lkmer_out_file)
        print(f"\nTemporary local kmers is created as {lkmer_out_file}\n")
        for tmp_files in temp_file_names:
            os.remove(tmp_files)
        print(f"\nExtracting global kmers from local kmers...\n")
        
        global_kmers_txt(lkmer_out_file, global_kmer_outfile)
        print(f"\nRemoving temporary local kmers file {lkmer_out_file}\n")
        os.remove(lkmer_out_file)
        print(f"\nTemporary local kmers file {lkmer_out_file} removed\n")

    print(f"\nKmer extraction completed successfully.\n")
    
    # For the kmer ratio summary statistics
    kmer_ratio_outfile = f"{name}_kmer_ratio.txt"
    per_base_outfile = f"{name}_per_base_ratio.txt"
    print(f"\nStarting kmer ratio calculation for {sequence_file} with locus {locus_file}, global kmer file {global_kmer_outfile} and canonical {canonical}\n")
    
    if jellyfish:
        threads = args.threads
        kmer_ratio_jf(sequence_file, locus_file, global_kmer_outfile, canonical, kmer_size, threads, kmer_ratio_outfile)
    else:
        kmer_ratio_txt(sequence_file, locus_file, global_kmer_outfile, canonical, kmer_ratio_outfile)
        print(f"\nKmer ratio calculation completed. Now calculating per base ratio\n")

    print(f"Repeat parameter to be used: {repeat_ratio}\n")
    jf = True if jellyfish else False
    per_base_ratio(sequence_file, kmer_ratio_outfile, jf, global_kmer_outfile, canonical, repeat_ratio, per_base_outfile)
    print(f"\nPer base ratio calculation completed. Output written to {per_base_outfile}\n")

    print(f"\nPrinting global kmers general statistics\n")
    
    if jellyfish:
        gkmer_info_jf(global_kmer_outfile)
    else:
        gkmer_info_txt(global_kmer_outfile)

    plot = args.plot
    if plot:
        slide_window = args.slide_window
        print(f"\nGenerating plots for mean and median\n")
        plotting_data(per_base_outfile, slide_window, name)
    else:
        print(f"\nNo plots will be generated\n")

    # Removing all the temporary files
    if os.path.isfile("parsed_output_kmer.txt"):
        print(f"Temporary parsed locus k-mer file found\nRemoving it...\n")
        os.remove("parsed_output_kmer.txt")
        print(f"File removed\n")
    if os.path.isfile("temp.file.fa"):
        print(f"Temporary parsed sequence file found\nRemoving it...\n")
        os.remove("temp.file.fa")
        print(f"File removed\n")
    if os.path.exists("tmp_parsed_kmer.jf"):
        print(f"Temporary parsed jf file found\nRemoving it...\n")
        os.remove("tmp_parsed_kmer.jf")
        print(f"File removed\n")
    if os.path.exists("temp_query_output.txt"):
        print(f"Temporary file found\nRemoving it...\n")
        os.remove("temp_query_output.txt")
        print(f"File removed\n")
    if os.path.exists("temp.query.fasta"):
        print(f"Temporary file found\nRemoving it...\n")
        os.remove("temp.query.fasta")
        print(f"File removed\n")
    
    # Finishing the program
    print(f"\nKmer Stat program ran successfully!\n")
    print(f"{datetime.datetime.now()} INFO: Finished running kmers_stat from kmerRRR tool\n")
    print(f"Processing time: {time.time() - start_time:.2f} seconds")