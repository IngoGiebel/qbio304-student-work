#! /usr/bin/env bash

# This shell script must be executed with the data folder as the current working directory.
# Prerequisite: fastqc is installed within the execution path.

fastqc *_f1.trim.p.fastq.gz *_r2.trim.p.fastq.gz -t 8
