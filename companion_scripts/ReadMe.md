# Documentation on how to use the companion scripts.

There are three scripts in this folder. pbr_suite.py is used to manipulate a per base ratio file to generate a more stringent LUCs for a locus when the original per base ratio file generated with a relaxed thresholds. This script is also used to create locus-unique and unique BED files that can be used with SAMtools to retain reads overalapping the regions within the BED file. This script was used to generate BED files to retain small RNA-seq data for the manuscript. 
The options for pbr_suite.py are:
```
pbr_suite.py [-h] -i INPUT_FILE [-nt NEW_THRESHOLD] -n NAME [-pbr] [-lu] [-u] [-cut CUT_OFF_VAL]
                    [-ctype CUT_TYPE]

pbr_suite

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        pbr_suite input file
  -nt NEW_THRESHOLD, --new_threshold NEW_THRESHOLD
                        new threshold value
  -n NAME, --name NAME  output file name
  -pbr, --per_base_file
                        Making new pbr file
  -lu, --lunique_file   Extracting local unique regions from pbr file in BED format
  -u, --unique_file     Extracting unique regions from pbr file in BED format
  -cut CUT_OFF_VAL, --cut_off_val CUT_OFF_VAL
                        cutoff value
  -ctype CUT_TYPE, --cut_type CUT_TYPE
                        cutoff type
```

The second pythonscript files in the LUCmean.py file, which is used along with the automate.ktest.sh bash script to assess which values of k is optimal for a given genome or locus. To run the bash script along with the LUCmean.py use:

```
bash automate.ktest.sh <path/to/sequence_file.fasta> <path/to/bed_file.bed> <genome_size> <threads> <kmerRRR/cmpanion_scripts/LUCmean.py>
```

This will run the script for k values: 17 and 21-181 with an increment of 10. To chnage the values of k, you should chnage it withint the script. The output will generate folders for each k value used. In addition, LUCmean.py will generate a spearate folder with the mean values of locus-unique classifications across each locus. 

Users can also use the LUCmean.py script seperately if they have more than one per base ratio files with different values of k. 
The options for using LUCmean.py is:
```
LUCmean [-h] [-bed BED_FILE] [-f FILES [FILES ...]] [-n NAME] [-p] [--plot_file PLOT_FILE]

options:
  -h, --help            show this help message and exit
  -bed BED_FILE, --bed_file BED_FILE
                        Name or Path of the BED locus file
  -f FILES [FILES ...], --files FILES [FILES ...]
                        Name or path of the per base ratio files
  -n NAME, --name NAME  Name of the outputfile, default: output_LUCmean.txt
  -p, --plot            PLotting the data
  --plot_file PLOT_FILE
                        Name or path of the file with which plotting will be done
```
