#!/bin/bash

sample=$1
outdir=$2

vawk="/opt/hall-lab/python-2.7.15/bin/python2 /opt/hall-lab/vawk/vawk"
bgzip=/opt/hall-lab/htslib-1.9/bin/bgzip
bedtools=/opt/hall-lab/bin/bedtools
zjoin="/opt/hall-lab/python-2.7.15/bin/python2 /opt/hall-lab/io/zjoin"


vcf3=$outdir/$sample/diploid/$sample.merged.3.100.vcf
dd=`dirname $vcf3`
mkdir -p $dd/annot/logs


cat $vcf3 \
  | $vawk  'BEGIN{OFS="\t"; FS="\t"}{
   if(I$SVTYPE!="BND") {
     end=I$END;
     start=$2;
     svlen=end-start;
     if(end<$2) {temp=end; end=start; start=temp;}
     ct=split(I$SNAME, spl, ",");
     for (ii=1; ii<=ct; ii++) {
       split(spl[ii], spl1, ":"); 
       split(spl1[2], spl2, "_"); 
       print $1, start, end+1, $3, I$SVTYPE, 200, svlen, spl2[1];
     }
   }
}' \
| sort -k1,8 -u  \
> $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.bed
#
#
#
cat $vcf3 \
  | grep -v SECONDARY \
  | $vawk  'BEGIN{OFS="\t"; FS="\t"}{
   if(I$SVTYPE=="BND") {
     start=$2;
     svlen=I$SVLEN;
     split($5, spl, /[\[\]:]/);
     end=spl[3];
     chr2=spl[2];
     if ($1!=chr2) svlen="inter";
     else svlen=end-start;
     ct=split(I$SNAME, spl, ",");
     for (ii=1; ii<=ct; ii++) {
       split(spl[ii], spl1, ":"); 
       split(spl1[2], spl2, "_"); 
       print $1, start-20, start+20, chr2, end-20, end+20, $3, I$SVTYPE, 200, svlen, spl2[1];
     }
  }
}' \
| sort -k1,11 -u  \
> $dd/annot/$sample.post_merge.all.round3.bnd.padded.bed
#
#
#
#
cat $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.bed \
| sort -k1,7 \
| $bedtools groupby -g 1,2,3,4,5,6,7 -c 8,8 -o count_distinct,distinct \
> $dd/annot/$sample.post_merge.all.round3.no_bnd.callers.bed
#
cat $dd/annot/$sample.post_merge.all.round3.bnd.padded.bed \
| sort -k1,10 \
| $bedtools groupby -g 1,2,3,4,5,6,7,8,9,10 -c 11,11 -o count_distinct,distinct \
> $dd/annot/$sample.post_merge.all.round3.bnd.callers.bed
#
#
cat $dd/annot/$sample.post_merge.all.round3.no_bnd.callers.bed \
| awk 'BEGIN{OFS="\t"; FS="\t"}{
  if($3<$2) {temp=$2; $2=$3; $3=temp;}
  if($2<0) $2=0;
  print $0, $5"_"$1"_"$2"_"$3;}' \
> $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.id.bed
#
cat $dd/annot/$sample.post_merge.all.round3.bnd.callers.bed \
| awk 'BEGIN{OFS="\t"; FS="\t"}{
  if( $9==200) {
    if($3<$2) {temp=$2; $2=$3; $3=temp;}
    if($2<0) $2=0;
    if($6<$5) {temp=$5; $5=$6; $6=temp;}
    if($5<0) $5=0;
    print $1, $2, $3, $7, $8, $9, $10, $11, $12,  $8"_"$1"_"$2"_"$3"_"$4"_"$5"_"$6;
    print $4, $5, $6, $7, $8, $9, $10, $11, $12,  $8"_"$1"_"$2"_"$3"_"$4"_"$5"_"$6;
  }
}' > $dd/annot/$sample.post_merge.all.round3.bnd.padded.id.bed
#
#
#
#
cat $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.id.bed \
| awk 'BEGIN{OFS="\t"; FS="\t"}{
  print $1, $2-20, $2+20, $10;
  print $1, $3-10, $3+20, $10;
}' \
| awk 'BEGIN{OFS="\t"}{if($2<0) $2=0; print $0;}' \
| $bedtools  annotate -counts  -i  stdin -files annot_tracks/str.bed annot_tracks/segdup.bed annot_tracks/alpha_sat.bed \
| sort -k4,4 | $bedtools groupby -g 4 -c 5,5,6,6,7,7 -o min,max,min,max,min,max \
| $zjoin -a $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.id.bed -b stdin -1 10 -2 1 \
| cut -f -10,12- \
| $bedtools annotate -i stdin -files annot_tracks/str.bed annot_tracks/segdup.bed annot_tracks/alpha_sat.bed \
> $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.id.annot.bed


cat $dd/annot/$sample.post_merge.all.round3.bnd.padded.id.bed \
| awk 'BEGIN{OFS="\t"; FS="\t"}{
  print $1, $2, $3, $10;
}' | $bedtools  annotate -counts  -i  stdin -files annot_tracks/str.bed annot_tracks/segdup.bed annot_tracks/alpha_sat.bed \
| sort -k4,4 | $bedtools groupby -g 4 -c 5,5,6,6,7,7 -o min,max,min,max,min,max \
| $zjoin -a $dd/annot/$sample.post_merge.all.round3.bnd.padded.id.bed -b stdin -1 10 -2 1 \
| cut -f -10,12- \
| awk 'BEGIN{OFS="\t"}{print $0, 0, 0, 0;}' \
> $dd/annot/$sample.post_merge.all.round3.bnd.padded.id.annot.bed



cat $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.id.bed \
| awk 'BEGIN{OFS="\t"; FS="\t"}{
  print $1, $2-20, $2+20, $10;
  print $1, $3-10, $3+20, $10;
}' \
| awk 'BEGIN{OFS="\t"}{if($2<0) $2=0; print $0;}' \
| $bedtools  annotate -counts  -i  stdin -files ../str.new.bed ../vntr.new.bed ../segdup.bed ../alpha_sat.bed \
| sort -k4,4 | $bedtools groupby -g 4 -c 5,5,6,6,7,7,8,8 -o min,max,min,max,min,max,min,max \
| $zjoin -a $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.id.bed -b stdin -1 10 -2 1 \
| cut -f -10,12- \
| $bedtools annotate -i stdin -files ../str.new.bed ../vntr.new.bed ../segdup.bed ../alpha_sat.bed \
> $dd/annot/$sample.post_merge.all.round3.no_bnd.padded.id.annot.new.bed


cat $dd/annot/$sample.post_merge.all.round3.bnd.padded.id.bed \
| awk 'BEGIN{OFS="\t"; FS="\t"}{
  print $1, $2, $3, $10;
}' | $bedtools  annotate -counts  -i  stdin -files ../str.new.bed  ../vntr.new.bed ../segdup.bed ../alpha_sat.bed \
| sort -k4,4 | $bedtools groupby -g 4 -c 5,5,6,6,7,7,8,8 -o min,max,min,max,min,max,min,max \
| $zjoin -a $dd/annot/$sample.post_merge.all.round3.bnd.padded.id.bed -b stdin -1 10 -2 1 \
| cut -f -10,12- \
| awk 'BEGIN{OFS="\t"}{print $0, 0, 0, 0, 0;}' \
> $dd/annot/$sample.post_merge.all.round3.bnd.padded.id.annot.new.bed
