#! /usr/bin/env bash

# This shell script must be executed with the data folder as the current working directory.
# Prerequisite: kallisto is installed within the execution path.

kallisto index -i Oryza_nivara.Oryza_nivara_v1.0.cdna.all.index ref-genomes/oryza-nivara/Oryza_nivara.Oryza_nivara_v1.0.cdna.all.fa.gz &> kallisto-index-Oryza_nivara.log
kallisto index -i Oryza_sativa.IRGSP-1.0.cdna.all.index ref-genomes/oryza-sativa/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz &> kallisto-index-Oryza_sativa.log
