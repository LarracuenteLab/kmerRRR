#!/usr/bin/env python3
"""This code is to find all the kmers in a genome assembly file.
It can handle fastq or fasta file, gzipped or not.
It extracts the local kmers at first for each contig
and then merges all the kmers to calculate the total kmers.
It can output canonical kmers if specified.

It takes the following input files:
1. A genome assembly file
2. kmer size based on the user
3. canonical value (0 or 1)
4. name of the program that would be also the output file name"""

import sys
import gzip
import os
from collections import defaultdict
from Bio import SeqIO 
import io
import shutil
import time
import datetime
import subprocess
import tempfile

# This function reads a sequence file in FASTA or FASTQ format
# can be gzipped or not and return a dictionary with contig IDs as keys
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
# This is also a stand alone program that the user can use
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

#Jellyfish implementation for global_kmers
def global_kmers_jf(sequence, canonical, kmer_size, genome_size, threads, output_file):
    kmer_length = int(kmer_size)
    hash_size_command = f"{genome_size}M"
    thread_command = int(threads)
    out_command = f"{output_file}"  
    canonical_command = "False"
    seq_command = sequence
    if sequence.endswith(".gz"):
        zcat_process = subprocess.Popen(["zcat", f"{seq_command}"], stdout=subprocess.PIPE)
    jellyfish_command = ["jellyfish", "count", "-m", f"{kmer_length}", "-t", f"{thread_command}", "-s", f"{hash_size_command}",  "-o", f"{out_command}"]
    if canonical == 1:
        canonical_command = "True"
        jellyfish_command.append("-C")
    if sequence.endswith(".gz"):
        try:
            print(f"Global k-mers counting for {sequence} started\n")
            print(f"Options are:\nCanonical:{canonical_command}\nk-mer size: {kmer_length}\nhash table size: {hash_size_command}, threads: {thread_command}\n")
            print(f"Output file generated: {output_file}")
            subprocess.run(jellyfish_command, stdin=zcat_process.stdout, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occured while calling jellyfish: {e}")
    else:
        jellyfish_command.append(seq_command)
        try:
            print(f"Global k-mers counting for {sequence} started\n")
            print(f"Options are:\nCanonical:{canonical_command}\nk-mer size: {kmer_length}\nhash table size: {hash_size_command}, threads: {thread_command}\n")
            print(f"Output file generated: {output_file}")
            subprocess.run(jellyfish_command, check= True)
        except subprocess.CalledProcessError as e:
            print(f"Error occured while calling jellyfish: {e}")

#To get global_kmer stats
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
    print(f"{datetime.datetime.now()} INFO: global_kmers has started\n")
    start_time = time.time()

    sequence_file = args.sequence_file
    if not os.path.exists(sequence_file):
        print(f"Sequence file: {sequence_file} does not exist. Program exiting...")
        sys.exit(1)

    kmer_size = args.kmer_size
    canonical = args.canonical
    name = args.name

    genome_size = args.genome_size
    jellyfish = args.jellyfish
    
    print(f"Input sequence file: {sequence_file}, Kmer size: {kmer_size}, Canonical: {canonical}, Output name: {name}\n")

    if jellyfish:
        print(f"\nJellyfish will be used to count the global k-mers\n")
    else:
        print(f"\nPython's dictionary or hash algorithm will be used to count the global k-mers\n")

    #For extracting the kmers from the sequence
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
    if not jellyfish:
        temp_file_names = contig_kmers(sequence_file, kmer_size, canonical, lkmer_out_file)
        print(f"\nTemporary local kmers is created as {lkmer_out_file}\n")
        for tmp_files in temp_file_names:
            os.remove(tmp_files)
        print(f"\nExtracting global kmers from local kmers...\n")
    
    if jellyfish:
        global_kmers_jf(sequence_file, canonical, kmer_size, genome_size, threads, global_kmer_outfile)
        gkmer_info_jf(global_kmer_outfile)
    else:
        global_kmers_txt(lkmer_out_file, global_kmer_outfile)
        print(f"\nRemoving temporary local kmers file {lkmer_out_file}\n")
        os.remove(lkmer_out_file)
        print(f"\nTemporary local kmers file {lkmer_out_file} removed\n")
        gkmer_info_txt(global_kmer_outfile)
    
    print(f"\nKmer extraction completed successfully.\n")
    print(f"global_kmers program ran successfully!\n")

    print(f"Processing time: {time.time() - start_time:.2f} seconds")