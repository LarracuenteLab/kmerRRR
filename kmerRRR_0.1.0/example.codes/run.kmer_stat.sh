#!/bin/bash
#SBATCH -J job.name
#SBATCH -o job.name.log.out%j
#SBATCH -e job.name.log.err%j
#SBATCH -c 1 #No multithreading incorporated in any of the scripts
#SBATCH -t 6:00:00 
#SBATCH -A alarracu_lab #only if you are using charlesworths as partition
#SBATCH -p charlesworths

#Call the program you want to use from kmerRRR
user=$USER

kmerRRR kmers_stat -seq "../test.files/test.fasta" -bed "../test.files/contig.bed" -c 1 -k 61 -n "../results.example/test.${user}" -g 200 --plot &> "../results.example/test.${user}.log"