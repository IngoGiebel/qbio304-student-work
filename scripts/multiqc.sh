#! /usr/bin/env bash

# This shell script must be executed with the data folder as the current working directory.
# Prerequisite: multiqc is installed within the execution path.

multiqc --force --dirs --filename MultiQC-Report.html --outdir ../results/multiqc --export --interactive .
