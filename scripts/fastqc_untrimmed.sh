#! /usr/bin/env bash

# This shell script must be executed with the data folder as the current working directory.
# Prerequisite: fastqc is installed within the execution path.

fastqc *_f1.fastq.gz *_r2.fastq.gz -t 8
