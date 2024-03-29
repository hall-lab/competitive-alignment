#!/bin/bash

sample=$1
hap=$2
vcf1=$3
outdir=$4

vawk="/opt/hall-lab/python-2.7.15/bin/python2 /opt/hall-lab/vawk/vawk"
bgzip=/opt/hall-lab/htslib-1.9/bin/bgzip
bedtools=/opt/hall-lab/bin/bedtools
zjoin="/opt/hall-lab/python-2.7.15/bin/python2 /opt/hall-lab/io/zjoin"


dd=$outdir/$sample/$hap
mkdir -p $dd/filter/logs

exclude=assembly_cov/$sample/$sample.$hap.exclude.bed 


cat $vcf1 \
  | $vawk  'BEGIN{OFS="\t"; FS="\t"}{
   end=I$END;
   start=$2;
   svlen=I$SVLEN
   if(I$SVTYPE=="BND") {
     split($5, spl, /[\[\]:]/);
     end=spl[3];
     svlen=end-start;
   }
   if(end<$2) {temp=end; end=start; start=temp;}
   print $1, start, end+1, $3, I$SVTYPE, 20, svlen;
}' \
> $dd/filter/$sample.padded.bed


cat $dd/filter/$sample.padded.bed \
| awk 'BEGIN{OFS="\t"; FS="\t"}{
  print $1, $2, $2+2, $4;
  print $1, $3, $3+2, $4;
}' | $bedtools  annotate -counts  -i  stdin -files $exclude \
| sort -k4,4 | $bedtools groupby -g 4 -c 5,5 -o min,max \
| $zjoin -a $dd/filter/$sample.padded.bed -b stdin -1 4 -2 1 \
| cut -f -7,9- \
| $bedtools annotate -i stdin -files $exclude \
> $dd/filter/$sample.annot.bed


cat $dd/filter/$sample.annot.bed | awk '{if($8>0 || $10>0.5) print $0;}' \
> $dd/filter/$sample.annot.filtered.bed

cat $vcf1 | $zjoin -p "#" -v -a stdin -b $dd/filter/$sample.annot.filtered.bed -1 3 -2 4 \
> $dd/sv.$sample.$hap.vcf

