#!/bin/bash

sample=$1
outdir=$2

svtools=/opt/hall-lab/python-2.7.15/bin/svtools

 
for caller in paf svimasm sv  pav
do
  echo "$outdir/$sample/pat/$caller.$sample.pat.vcf"
  echo "$outdir/$sample/mat/$caller.$sample.mat.vcf"
done > $outdir/$sample/diploid/var.list

for caller in pbsv sniffles svim
do
  echo "$outdir/$sample/diploid/$caller.$sample.diploid.vcf"
done >> $outdir/$sample/diploid/var.list


pad=20
$svtools lsort -f $outdir/$sample/diploid/var.list -r > $outdir/$sample/diploid/$sample.sorted.1.vcf
cat  $outdir/$sample/diploid/$sample.sorted.1.vcf | $svtools lmerge -i /dev/stdin -f $pad > $outdir/$sample/diploid/$sample.merged.1.$pad.vcf
echo "$outdir/$sample/diploid/$sample.merged.1.$pad.vcf"  > $outdir/$sample/diploid/round2.list
$svtools lsort -f $outdir/$sample/diploid/round2.list -r > $outdir/$sample/diploid/$sample.sorted.2.vcf
cat $outdir/$sample/diploid/$sample.sorted.2.vcf | $svtools lmerge -i /dev/stdin -f 100 -w carrier_wt > $outdir/$sample/diploid/$sample.merged.2.100.vcf
