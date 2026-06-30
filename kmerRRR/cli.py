#!/usr/bin/env python3

import argparse
import sys
import time
import datetime
from kmerRRR.kmers_stat import main as kmers_stat_main 
from kmerRRR.summary_stat import main as summary_stat_main
from kmerRRR.global_kmers import main as global_kmers_main
from kmerRRR.get_mapq import main as get_mapq_main
from kmerRRR.bam_file_manipulate import main as bam_file_manipulate_main
from kmerRRR.gkmer_info import main as gkmer_info
from kmerRRR.kmer_dump import main as kmer_dump
from kmerRRR.plotting_data import main as plot_data
from kmerRRR.create_bedgraph import main as create_bed
from kmerRRR import __version__

def main():
    parser = argparse.ArgumentParser(prog="kmerRRR")
    #subcommand for version
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}", help= "Print version of kmerRRR")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # subcommand for kmers_stat
    parser_kmers_stat = subparsers.add_parser("kmers_stat", help= "To extract global k-mer counts and to calculate summary statistics from a bed file with locus information")
    parser_kmers_stat.add_argument("-seq", "--sequence_file", required=True, help="Input genome assembly file (fasta or fastq format), can be gzipped")
    parser_kmers_stat.add_argument("-k", "--kmer_size", type=int, required=True, help="k-mer size")
    parser_kmers_stat.add_argument("-c", "--canonical", type=int, choices=[0, 1], default=1, help="Whether to output canonical kmers (0 for no, 1 for yes), default: 0")
    parser_kmers_stat.add_argument("-bed", "--bed_file", required=True, help="Path or name to the locus file in BED format (see example files in test.files folder)")
    parser_kmers_stat.add_argument("-n", "--name", required=True, help="Name for the output files")
    parser_kmers_stat.add_argument("-p", "--plot", action="store_true", help="Generate plot of mean and median from per base ratio file, default: False")
    parser_kmers_stat.add_argument("-jf", "--jellyfish", action="store_true", help="To use jellyfish to count k-mers instead of using the python's dictionary (default = False)")
    parser_kmers_stat.add_argument("-g", "--genome_size", required= True, type= int, help = "Genome size of the organism in Mb")
    parser_kmers_stat.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use, only use when using jellyfish")
    parser_kmers_stat.add_argument("-sw", "--slide_window", type=int, default=10000, help="Sliding window for plotting per base ratio file, default is 10kb")
    parser_kmers_stat.add_argument("-rr", "--repeat_ratio", type=float, default=1, help="Repeat parameter to state repetitiveness, default = 1")
    parser_kmers_stat.set_defaults(func=kmers_stat_main)

    # subcommand for summary_stat
    parser_summary_stat = subparsers.add_parser("summary_stat", help= "To calculate summary statistics from a bed file with locus information, given that you already have global k-mer counts from kmers_stat or global_kmers program")
    parser_summary_stat.add_argument("-seq", "--sequence_file", required=True, help="Path or name to the sequence file in FASTA or FASTQ format, can be gzipped")
    parser_summary_stat.add_argument("-bed", "--bed_file", required=True, help="Path or name to the locus file in BED format (see example files in test.files folder)")
    parser_summary_stat.add_argument("-c", "--canonical", type=int, choices=[0, 1], default=1, help="Whether to consider canonical k-mers, only if global k-mers were counted with canonical option, default: 0")
    parser_summary_stat.add_argument("-gk", "--global_kmer_file", required=True, help="Path or name to the global kmer file with kmer counts")
    parser_summary_stat.add_argument("-n", "--name", required=True, help="Path or name for the output files")
    parser_summary_stat.add_argument("-p", "--plot", action="store_true", help="Generate plot of mean and median from per base ratio file, default: False")
    parser_summary_stat.add_argument("-k", "--kmer_size", type=int, required=True, help="k-mer size")
    parser_summary_stat.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use, only use when using jellyfish")
    parser_summary_stat.add_argument("-sw", "--slide_window", type=int, default=10000, help="Sliding window for plotting per base ratio file, default is 10kb")
    parser_summary_stat.add_argument("-rr", "--repeat_ratio", type=float, default=1, help="Repeat parameter to state repetitiveness, default = 1")
    parser_summary_stat.set_defaults(func=summary_stat_main)

    # subcommand for global_kmer
    parser_global_kmers = subparsers.add_parser("global_kmers", help= "To extract k-mers and their counts from genome assembly, if kmers_stat program was not run before or to run kmer counting again with different options separately")
    parser_global_kmers.add_argument("-seq", "--sequence_file", required=True, help="Path or name to the sequence file in FASTA or FASTQ format, can be gzipped")
    parser_global_kmers.add_argument("-k", "--kmer_size", type=int, required=True, help="K-mer size")
    parser_global_kmers.add_argument("-c", "--canonical", type=int, choices=[0, 1], default=0, help="Whether to output canonical kmers (0 for no, 1 for yes), default: 0")
    parser_global_kmers.add_argument("-n", "--name", required=True, help="Path or name for the output files")
    parser_global_kmers.add_argument("-jf", "--jellyfish", action="store_true", help="To use jellyfish to count k-mers instead of using the python's dictionary (default = False)")
    parser_global_kmers.add_argument("-g", "--genome_size", required= True, type= int, help = "Genome size of the organism in Mb")
    parser_global_kmers.add_argument("-t", "--threads", type=int, help="Number of threads to use, only use when using jellyfish")
    parser_global_kmers.set_defaults(func=global_kmers_main)

    # subcommand for bam_file_manipulate
    parser_bam_file_manipulate = subparsers.add_parser("bam_file_manipulate", help= "Manipulating MAPQ score in bamfile based on mean or median k-mer per base summary statistics from kmers_stat or summary_stat output")
    parser_bam_file_manipulate.add_argument("-pbr", "--per_base_ratio", required= True, help="Path or name of per base kmer ratio file with the summary statistics")
    parser_bam_file_manipulate.add_argument("-bam", "--bam_file", required= True, help= "Path or name of the bam file to be manipulated")
    parser_bam_file_manipulate.add_argument("-n", "--name", required= True, help="Path or name for the output files")
    parser_bam_file_manipulate.add_argument("-cut", "--cut_off", nargs='?', type= float, default = 1, help= "Cut-off value to manipulate the bamfile, default: 1 for mean and 1 for median")
    parser_bam_file_manipulate.add_argument("-ctype", "--cut_type", type=str, nargs= '?', default= "mean", help= "Use either mean or median to manipulate the bamfile, default: mean")
    parser_bam_file_manipulate.add_argument("-alt", "--alt_bam", action='store_true', help="To output a bamfile without contigs present in the locus file for future merge and use in the downstream analysis, default: False")
    parser_bam_file_manipulate.add_argument("--old_mapq", nargs='?', type=int, default=30, help= "Minimum MAPQ score required for manipulating the read. Reads with score below this will be manipulated, default: 30")
    parser_bam_file_manipulate.add_argument("--new_mapq", nargs='?', type=int, default=40, help= "Manipulated new MAPQL score, default: 40")
    parser_bam_file_manipulate.add_argument("-sd", "--seed", type = int, help = "Seed value for sorting multimappers randomly")
    parser_bam_file_manipulate.set_defaults(func=bam_file_manipulate_main)
 
    # subcommand for get_mapq
    parser_get_mapq = subparsers.add_parser("get_mapq", help= "To extract MAPQ scores from a bamfile using with or without locus-based bed file")
    parser_get_mapq.add_argument("-bam", "--bam_file", required= True, help= "Path or name of the bamfile")
    parser_get_mapq.add_argument("-bed", "--bed_file", required= True, help= "Path or name of the bedfile")
    parser_get_mapq.add_argument("-n", "--name", required= True, help= "Path or name for the output files")
    parser_get_mapq.add_argument("-locus", "--locus_based", action="store_true", help= "Whether MAPQ will be extracted based on locus-based BED file, default: False")
    parser_get_mapq.add_argument("-p", "--plot", action="store_true", help="Generate plot of MAPQ scores per contig from mapq file, default: False")
    parser_get_mapq.add_argument("-sw", "--slide_window", type=int, default=10000, help="Sliding window for plotting per base ratio file, default is 10kb")
    parser_get_mapq.set_defaults(func=get_mapq_main)

    # subcommand for kmer info
    parser_gkmer_info = subparsers.add_parser("gkmer_info", help= "To get the stat report from global k-mer databse")
    parser_gkmer_info.add_argument("-gk", "--global_kmer_file", required = True, help = "Path or name of the global k-mer file (or folder if it's a database)")
    parser_gkmer_info.set_defaults(func=gkmer_info)

    # subcommand for kmer dump
    parser_kmer_dump = subparsers.add_parser("kmer_dump", help= "To extract all the global k-mers from the lmdb file")
    parser_kmer_dump.add_argument("-jf", "--global_kmer_jf", required= True, help= "Path or name of the global k-mer database folder")
    parser_kmer_dump.add_argument("-n" "-name", required= True, help = "Path or name for the output file")
    parser_kmer_dump.set_defaults(func=kmer_dump)

    # subcommand for plotting_data
    parser_plot_data = subparsers.add_parser("plotting_data", help= "To plot the data from per base ration file")
    parser_plot_data.add_argument("-pbr", "--per_base_ratio", help= "Path or name of the per base ratio file")
    parser_plot_data.add_argument("-mapq", "--mapq_file", help="Path or name of the mapq file generated by get_mapq script")
    parser_plot_data.add_argument("-n", "--name", required=True, help="Name or path of the output files")
    parser_plot_data.add_argument("-sw", "--slide_window", type=int, default=10000, help="Sliding window for plotting per base ratio file, default is 10kb")
    parser_plot_data.set_defaults(func=plot_data)

    # subcommand for create_bedgraph
    parser_create_bed = subparsers.add_parser("create_bedgraph", help="Create a bedgraph format file from per base ratio file")
    parser_create_bed.add_argument("-pbr", "--per_base_ratio", help= "Path or name of the per base ratio file")
    parser_create_bed.add_argument("-mapq", "--mapq_file", help="Path or name of the MAPQ fil generae from get_mapq")
    parser_create_bed.add_argument("-sw", "--slide_window", type=int, default=10000, help="Sliding window for plotting per base ratio file, default is 10kb")
    parser_create_bed.add_argument("-n", "--name", required=True, help="Name or path of the output files")
    parser_create_bed.set_defaults(func=create_bed)

   
    # Calling args
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    print(f"{datetime.date.today()}")
    start_time = time.time()
    main()
    print(f"Time required: {time.time() - start_time:.2f} seconds")
    sys.exit(0)