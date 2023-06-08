#!/bin/bash

sample=$1
mat=$2
pat=$3

cat $mat | grep  "^>" | sed 's/>//g' > $sample/$sample.mat.contigs.txt
cat $pat | grep  "^>" | sed 's/>//g' > $sample/$sample.pat.contigs.txt
