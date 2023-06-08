#!/bin/bash

inbam=$1
contig=$2
outdir=$3
base=`basename $inbam .bam`

samtools=/opt/hall-lab/samtools*/bin/samtools


outbam=$outdir/$base.$contig.bam
$samtools view -h $inbam \
| awk -v contig="$contig" '{
  if($1~/@/) print $0;
  else if ($1==contig) print $0;
}' | $samtools view -S -b - > $outbam


