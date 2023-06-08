#!/bin/bash

sample=$1
outdir=$2

mkdir -p $outdir/$sample/diploid/annot

onekg=vcfs_clean_filtered_0811/$sample/diploid/$sample.ill.vcf
svtools=/opt/hall-lab/python-2.7.15/bin/svtools

echo "$outdir/$sample/diploid/$sample.merged.2.100.vcf" > $outdir/$sample/diploid/round3.list
echo "$onekg" >> $outdir/$sample/diploid/round3.list

$svtools lsort -f $outdir/$sample/diploid/round3.list -r > $outdir/$sample/diploid/$sample.sorted.3.vcf
cat $outdir/$sample/diploid/$sample.sorted.3.vcf | $svtools lmerge -i /dev/stdin -f 100 -w carrier_wt > $outdir/$sample/diploid/$sample.merged.3.100.vcf
