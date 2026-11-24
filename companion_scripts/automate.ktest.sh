#!/usr/bin/bash

# This script is to use to automate k test whereby
# users can estimate what value(s) of k might be useful
# for their locus of interest

date

sequence_file=$1
bed_file=$2
genome_size=$3
threads=$4
python_script=$5

if [[ $# -lt 5 ]]; then
    echo -e "\nUsage: bash automate.ktest.sh [sequence_file] [bed_file] [genome_size] [threads] [path to LUCmean.py]\n"
    exit 1
else
    echo -e "\n5 arguments given\n"
fi

if [[ -f  "$sequence_file" ]]; then
    echo -e "\nFound $sequence_file\n"
else
    echo -e "\n$equence_file not found!!!\n"
    exit 1
fi

if [[ -f  "$bed_file" ]]; then
    echo -e "\nFound $bed_file\n"
else
    echo -e "\n$bed_file not found!!!\n"
    exit 1
fi

if [[ -f "$python_script" ]]; then
    echo -e "\n Found $python_script\n"
else
    echo -e "\n$python_script not found\n"
    exit 1
fi

declare -a pbr_list
for k in 21 31 41 51 61 71 81 91 101 111 121 131 141 151 161 171 181
    do 
    echo -e "\nStarting ${k}\n" 
    mkdir -p ${k}mer 
    kmerRRR kmers_stat -seq "$sequence_file" -bed "$bed_file" -jf -k ${k} -n "${k}mer/${k}mer.ktest" -c 1 -t "$threads" -g "$genome_size"
    file_name="${k}mer/${k}mer.ktest_per_base_ratio.txt"
    if [[ -f  "$file_name" ]]; then
        echo -e "\nFound $file_name\n"
        pbr_list+=("$file_name")
    else
        echo -e "\n$file_name not found!!!\n"
    fi
    echo -e "\nFinished with ${k}\n" 
done 

echo -e "\nFinished running kmers_stat from kmerRRR\nStarting to work on LUCmean.py script to calculate and plot kmean per locus\n"

python3 "$python_script" -bed "$bed_file" -f "${pbr_list[@]}" --plot

echo -e "\nFinished Running\n"