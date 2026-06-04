#!/usr/bin/env python3 

import subprocess #type: ignore
import gzip
import io
import os
import time 
import datetime
import sys

def kmer_dump(global_kmer_file, out_kmer_file):
    print(f"\nGLobal k-mers will be extracted from {global_kmer_file} database and written to {out_kmer_file}\n")
    Buffer_size = 1_000_000
    buffer = []
    counting_buffer = 0
    if global_kmer_file.endswith(".jf"):
        print(f"Global k-mers file contains k-mer counts counted using jellyfish")
    jellyfish_dump_command = ["jellyfish", "dump", "-c", f"{global_kmer_file}"]
    kmer_counts = subprocess.run(jellyfish_dump_command, capture_output=True, text=True).stdout
    with gzip.GzipFile(out_kmer_file, "w") as out_f:
        with io.TextIOWrapper(out_f, encoding= 'utf-8') as gkmer_out:
            gkmer_out.write(f"kmer\tcount\n")
            for lines in kmer_counts.strip().splitlines():
                key, count = lines.split()
                buffer.append(f"{key}\t{count}\n")
                if len(buffer) >= Buffer_size:
                    counting_buffer += 1
                    gkmer_out.writelines(buffer)
                    buffer.clear()
            if buffer:
                temp_count = len(buffer)
                gkmer_out.writelines(buffer)
                counting_buffer += temp_count/Buffer_size
    print(f"Total buffer chunks written: {counting_buffer}\n")
    print(f"Total k-mers extracted: {counting_buffer * Buffer_size}\n")

def main(args):
    start_time = time.time()
    print(f"{datetime.datetime.now()} INFO: kmer_dump has started\n")
    global_kmer_file = args.global_kmer_jf
    name = args.name
    outfile = f"{name}.txt.gz"
    if os.path.exists(global_kmer_file):
        print(f"{global_kmer_file} found..\nGlobal kmers extracting from {global_kmer_file} jellyfish data structure\n")
    else:
        print(f"{global_kmer_file} does not exist\nSystem exiting...")
        sys.exit(1)
    print(f"Starting extracting k-mers...\n")
    kmer_dump(global_kmer_file, outfile)
    print(f"Finished extracting k-mers\nGlobal kmers output in {outfile}\n")
    print(f"gkmers_dump ran successfully!\n")
    print(f"Processing time: {time.time() - start_time:.2f} seconds")