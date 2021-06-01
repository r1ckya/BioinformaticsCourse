#!/bin/bash
python gen_reads.py > reads.fastq
python de_brujin.py 50 < reads.fastq
