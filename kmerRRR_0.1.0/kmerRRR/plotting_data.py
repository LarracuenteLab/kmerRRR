#!/usr/bin/env python3

import time
import datetime
import subprocess
import importlib.resources as resources


# To plot the mean and median summary stat for per base ratio file
def plotting_pbr_data(per_base_file, swindow, name_plot):
    with resources.path("kmerRRR.R", "pbr_plot.r") as r_script:
        subprocess.run(["Rscript", str(r_script), str(per_base_file), str(swindow), str(name_plot)], check=True)

# To plot the MAPQ
def plotting_mapq_data(mapq_file, slide_window, name_plot):
    with resources.path("kmerRRR.R", "mapq_plot.r") as r_script:
        subprocess.run(["Rscript", str(r_script), str(mapq_file), str(slide_window), str(name_plot)], check=True)  

def main(args):
    start_time = time.time()
    print(f"{datetime.datetime.now()} INFO: kmers_stat has started\n")

    if args.per_base_ratio:
        pbr_file = args.per_base_ratio
    else:
        pbr_file = False
    if args.mapq_file:
        mapq_file = args.mapq_file
    else:
        mapq_file = False
    name = args.name

    if pbr_file:
        plotting_file = pbr_file
        slide_window = args.slide_window
        print(f"\nGenerating plots for mean and median\n")
        plotting_pbr_data(pbr_file, slide_window, name)
    if mapq_file:
        plotting_file = mapq_file
        slide_window = args.slide_window
        print(f"\nGenerating plots for MAPQ\n")
        plotting_mapq_data(mapq_file, slide_window, name)

    if not pbr_file and mapq_file:
        print(f"Error: Need either per base ratio or mapq file")
        
    if pbr_file and mapq_file:
        print(f"\nProgram to plot {pbr_file} and {mapq_file} ran successfully!\n")
    else:
        print(f"\nProgram to plot {plotting_file} ran successfully!\n")
    print(f"{datetime.datetime.now()} INFO: Finished running plotting_data from kmerRRR tool\n")
    print(f"Processing time: {time.time() - start_time:.2f} seconds")