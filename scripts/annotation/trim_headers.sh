#!/bin/sh

sample=$1
describer=$(echo ${sample} | sed 's/.pfam.domtblout//')

tail -n +4 ${describer}.pfam.domtblout > ${describer}.out
