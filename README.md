# tools for sRNA analysis

This repo is used to deposit some scripts I used in sRNA-related analysis



### (1) v1.fastq2fasta.py

This script is used to convert fastq file (single-ended Illumina reads) to fasta format, but with additional demands.

The sequence id of fasta file will be reformatted like below:

`seq<serial number>_<frequency>_<length>`

