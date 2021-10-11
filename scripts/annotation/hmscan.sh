#!/bin/sh

sample=$1
describer=$(echo ${sample} | sed 's/.fq.gz//')

hmmscan --cpu 1 --domtblout ${sample}.pfam.domtblout ~/bioinformatics/pfam/Pfam-A.hmm ${sample}
