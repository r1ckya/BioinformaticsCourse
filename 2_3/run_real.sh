#!/bin/bash

n=10000
k=10

head -n $((4 * n)) data.fastq > sample.fastq

java -jar trimmomatic-0.39.jar \
    SE \
    sample.fastq \
    trimmed.fastq \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:10:3 MINLEN:36

python de_brujin.py $k < trimmed.fastq
