# qbio304-student-work/scripts

This folder contains all script files written for the QBio304 (Applied
Bioinformatics) student work assigned in WS 2022/2023.


## Script execution order:

1. fastqc_untrimmed.sh  
Check the HTML quality reports. Consider changing the parameters before executing
the trimmomatic script.

2. trimmomatic  
Check the Trimmomatic output.

3. fastqc_trimmed.sh  
Check the HTML quality reports before continuing.

4. kallisto_index.sh

5. kallisto_quant.sh

6. multiqc.sh  
Check the MultiQC report.

7. dge-analysis-PRJCA004229.R


**Note:** all scripts must be executed with the data folder as the current working
directory. The programs invoked by these scripts must already be installed
within the execution path.



## List of scripts:

### dge-analysis-PRJCA004229.R

R script for performing a differential gene expression (DGE) analysis.

Prerequisite: the package "biomaRt" is installed via
BiocManager::install("biomaRt"). The other packages listed in the section
"Load required packages" are installed as usual.


### fastqc_trimmed.sh

Shell script to perform a FastQC analysis of all trimmed FASTQ files.
For each FASTQ file a corresponding "*fastqc.html" report is written. This
shell script must be executed with the data folder as the current working
directory.

Prerequisite: fastqc is installed within the execution path.


### fastqc_untrimmed.sh

Shell script to perform a FastQC analysis of all untrimmed (raw) FASTQ files.
For each FASTQ file a corresponding "*fastqc.html" report is written. This
shell script must be executed with the data folder as the current working
directory.

Prerequisite: fastqc is installed within the execution path.


### kallisto_index.sh

Shell script to create kallisto index files for the Oryza nivara and Oryza
sativa reference genes. Logs the output to kallisto-index-Oryza_nivara.log.
and kallisto-index-Oryza_sativa.log. This shell script must be executed with
the data folder as the current working directory.

Prerequisite: kallisto is installed within the execution path.


### kallisto_quant.sh

Shell script to execute "kallisto quant..." for all trimmed paired-end FASTQ
files. This shell script must be executed with the data folder as the current
working directory.

The kallisto program creates on folder [run accession] for each pair of
FASTQ files. Addionally a log file kallisto-quant-[run accession].log is
written.

Prerequisite: kallisto is installed within the execution path.


### multiqc.sh

Shell script to execute MultiQC on the current working directory. The MultiQC
report is written to the results/multiqc subdirectory. Any preexisting reports
are overwritten.

Prerequisite: multiqc is installed within the execution path.


### trimmomatic

Shell script. Performs a paired-end trimming of all of the illumina FASTQ
raw data files in the data subdirectory.

This script calls the trimmomatic java program within the

qbio304-student-work/programs/trimmomatic

subdirectory on all [run accession]_[f1|r2]_fastq.gz data files and creates
the trimmed files [run accession]_[f1|r2].trim.[p|u|_fastq.gz.

For details of the trimmomatic parameters see the script file. The output of
trimmomatic was (manually) stored in trimmomatic.log file within the data
subdirectory.

This shell script is written in Python3 and therefore requires Python3 to be
installed and accessible in the current environment (Python3 is called via the
shebang "#! /usr/bin/env python3").