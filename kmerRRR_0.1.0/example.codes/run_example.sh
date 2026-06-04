#!/bin/bash
#SBATCH -J job.name
#SBATCH -o job.name.log.out%j
#SBATCH -e job.name.log.err%j
#SBATCH -c 1 #No multithreading incorporated in any of the scripts
#SBATCH -t 6:00:00 
#SBATCH -A alarracu_lab #only if you are using charlesworths as partition
#SBATCH -p charlesworths
#SBATCH --mail-user=ALL
#SBATCH --user=jrahmat

source /path/to/your/directory/miniforge3/bin/activate #/path/to/your/directory/ is where you have installed miniforge3

source /path/to/your/directory/venv/bin/activate #/path/to/your/directory/ is where you have opened the venv in the beginning of this documentation

#Call the python script you want to use

kmerRRR <program_name> -seq <sequence.fasta> -bed <locus.bed> -c 1 -k 61 -n “name” –-plot -g 200 &> output.log