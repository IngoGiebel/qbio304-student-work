# qbio304-student-work/data
This folder contains all raw and processed data files used in the QBio304
(Applied Bioinformatics) student work assigned in WS 2022/2023:


## [run accession]/
Subfolders created by kallisto which contain the abundance.tsv files (plaintext
files of the abundance estimates).


## gene-sets-files/
GMT (Gene Matrix Transposed file format) files for different gene sets:

- Oryza_sativa_all.gmt  
Downloaded from
http://structuralbiology.cau.edu.cn/PlantGSEA/database/Osa.DetailInfo)

- Oryza_sativa_Cyc-orig.gmt  
PlantCyc gene sets downloaded from
http://structuralbiology.cau.edu.cn/PlantGSEA/database/Osa_Cyc

- Oryza_sativa_Cyc.gmt  
This file is a manually corrected version of Oryza_sativa_Cyc-orig.gmt because
the original file is not fully complient to the GMT file format (the original
file uses "," instead of tabs to separate the genes which belong to a gene set.

- Oryza_sativa_GFam.gmt  
Gene Family based gene sets downloaded from
http://structuralbiology.cau.edu.cn/PlantGSEA/database/Osa_GFam

- Oryza_sativa_GO.gmt  
GO (Gene Ontology) gene sets downloaded from
http://structuralbiology.cau.edu.cn/PlantGSEA/database/Osa_GO

- Oryza_sativa_KEGG.gmt  
KEGG gene sets downloaded from
http://structuralbiology.cau.edu.cn/PlantGSEA/database/Osa_KEGG

- Oryza_sativa_MIR.gmt  
MIR gene sets downloaded from
http://structuralbiology.cau.edu.cn/PlantGSEA/database/Osa_MIR

- Oryza_sativa_PO.gmt  
PO gene sets downloaded from
http://structuralbiology.cau.edu.cn/PlantGSEA/database/Osa_PO


## ref-genomes/oryza-nivara
Reference genome files for Oryza nivara (wild rice), cultivar BJ278,
downloaded from
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/fasta/oryza_nivara/cdna/


## ref-genomes/oryza-sativa
Reference genome files for Oryza sativa (cultivated rice), cultivar Nipponbare
(reference for Japonica group), downloaded from
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/fasta/oryza_sativa/cdna/


## [run accession]_[f1|r2].fastq.gz
Raw FASTQ files downloaded from
https://ngdc.cncb.ac.cn/gsa/browse/CRA003736

## [run accession]_[f1|r2].trim.[p|u].fastq.gz
Trimmed FASTQ files created by Trimmomatic.

## *_fastqc.html
Report files created by FastQC.

## *_fastqc.zip
Files created by MultiQC.

## kallisto-index-Oryza_nivara.log
kallisto log file of the index file creation for the Oryza nivara reference
genome.

## kallisto-index-Oryza_sativa.log
kallisto log file of the index file creation for the Oryza sativa reference
genome.

## kallisto-quant-[run accession].log
kallisto log file of the pseudoalignment run.

## md5checksums.txt
MD5 checksums for the downloaded raw FASTQ files. The integrety of the
downloaded files was checked by `md5sum -c md5checksums.txt`.

## Oryza_nivara.Oryza_nivara_v1.0.cdna.all.index
kallisto index file of the Oryza nivara reference genome.

## Oryza_sativa.IRGSP-1.0.cdna.all.index
kallisto index file of the Oryza sativa reference genome.

## studydesign-PRJCA004229.tsv
Studydesign file for the data.

## trimmomatic.log
Log file of the Trimmomatic program execution.