#!/bin/bash

python=/opt/hall-lab/python-2.7.15/bin/python

cromdir=/storage1/fs1/ccdg/Active/analysis/abelhj/pangenome/cromwell
wwcalls=/storage1/fs1/ccdg/Active/analysis/wen-wei.liao/callsets/Yale_HPP_Year1_Variant_Calls
pavdir=/storage1/fs1/ccdg/Active/analysis/abelhj/pangenome/pav_calls/20210518_PAV_HPRC

cat samps.txt |  awk  -v cromdir="$cromdir" -v wwdir="$wwcalls" -v pavdir="$pavdir" 'BEGIN{OFS="\t"}{
  sample=$1; snum=$2;
  print sample, snum, "mat", 1, "paf", cromdir"/small_var_vcfs/sv."sample".1.vcf.gz", "haploid_vcf_converters/doctor_paf.hap.py";
  print sample, snum, "pat", 2, "paf", cromdir"/small_var_vcfs/sv."sample".2.vcf.gz", "haploid_vcf_converters/doctor_paf.hap.py";
  print sample, snum, "mat", 1, "sv", cromdir"/sv_bedpes/sv."sample".1.bedpe", "haploid_vcf_converters/doctor_sv.hap.py";
  print sample, snum, "pat", 2, "sv", cromdir"/sv_bedpes/sv."sample".2.bedpe", "haploid_vcf_converters/doctor_sv.hap.py";
  print sample, snum, "diploid", 0, "pbsv", wwdir"/"sample"/"sample".GRCh38_no_alt.pbsv.vcf.gz", "haploid_vcf_converters/doctor_sniffles.py";
  print sample, snum, "diploid", 0, "svim", wwdir"/"sample"/"sample".GRCh38_no_alt.svim.q10.vcf.gz", "haploid_vcf_converters/doctor_sniffles.py";
  print sample, snum, "diploid", 0, "sniffles", wwdir"/"sample"/"sample".GRCh38_no_alt.sniffles.vcf.gz", "haploid_vcf_converters/doctor_sniffles.py";
  print sample, snum, "mat", 1, "pav", pavdir"/hprc_pav_"sample".vcf.gz", "haploid_vcf_converters/doctor_pav.hap.py";
  print sample, snum, "pat", 2, "pav", pavdir"/hprc_pav_"sample".vcf.gz", "haploid_vcf_converters/doctor_pav.hap.py";
  print sample, snum, "mat", 2, "svimasm", wwdir"/"sample"/"sample".GRCh38_no_alt.svim-asm.vcf.gz", "haploid_vcf_converters/doctor_sniffles.hap.py";
  print sample, snum, "pat", 1, "svimasm", wwdir"/"sample"/"sample".GRCh38_no_alt.svim-asm.vcf.gz", "haploid_vcf_converters/doctor_sniffles.hap.py";
}' > sample.map.haps



#converts vcfs and/or bedpe to svtools compatible formats
#note that ins are treated like dups for merging but this is reverted later.  ins are only merged based on position and length


outdir=vcfs_clean_filtered_0811

while read sample snum hap hapnum caller input prog
do
  mkdir -p $outdir/$sample/$hap
  mkdir -p $outdir/$sample/diploid
  printf "%s %s -i %s -c %s -p %s -r %s > %s/%s/%s/%s.%s.%s.vcf \n" $python $prog $input $caller $hapnum $hap $outdir $sample $hap  $caller $sample $hap
done < <(cat sample.map.haps | awk '{if($5!="sv") print $0;}') > make_vcfs.1.sh

while read sample snum hap hapnum caller input prog
do
  base=`basename $input .bedpe`
  pafdir=/storage1/fs1/ccdg/Active/analysis/abelhj/pangenome/cromwell/small_var_vcfs
  paf=$pafdir/$base.vcf.gz
  printf "%s %s -i %s -b %s  -p %s -r %s > %s/%s/%s/%s.%s.%s.unfiltered.vcf \n" $python $prog $paf $input $hapnum $hap $outdir $sample $hap  $caller $sample $hap
done < <(cat sample.map.haps | awk '{if($5=="sv") print $0;}') > make_vcfs.2.sh

while read sample
do
  mantavcf=../../1kg_allison/per_sample/renamed/redo/$sample.1kg.ill.vcf.gz
  outfile=$outdir/$sample/diploid/$sample.ill.vcf
  printf "%s haploid_vcf_converters/revert.manta.py -i %s > %s \n"  $python $mantavcf $outfile
done < <(cat sample.map.haps | cut -f 1 | sort -u) >> make_vcfs.2.sh

mkdir -p $outdir/logs


LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q general \
-oo $outdir/logs/tovcf.%J.log \
 -a 'docker(apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
"bash ./make_vcfs.1.sh"


LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q general \
-oo $outdir/logs/tovcf.%J.log \
 -a 'docker(apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
"bash ./make_vcfs.2.sh"



while read sample
do
  mkdir -p $outdir/$sample/$hap/logs
  LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub -G compute-ccdg \
  -oo $outdir/$sample/$hap/logs/make_vcfs.%J.log \
  -q general  -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
  "/bin/bash ./per_sample_merge.sh $sample $outdir"

done < <(cat sample.map.haps | cut -f 1 | sort -u ) 
##


#
while read sample 
do
  cat  $outdir/$sample/diploid/annot/$sample.post_merge.all.round3.bnd.padded.id.annot.bed $outdir/$sample/diploid/annot/$sample.post_merge.all.round3.no_bnd.padded.id.annot.bed \
   | awk -v sample="$sample"  'BEGIN{OFS="\t"}{
     print $0, sample;
   }' 
done < <( cat sample.map.haps | cut -f 1 | sort -u  )  | gzip -c > $outdir/all_samples.annot.2.bed.gz
exit 0

while read sample
do 
  cat  $outdir/$sample/diploid/annot/$sample.post_merge.all.round3.bnd.padded.id.annot.bed $outdir/$sample/diploid/annot/$sample.post_merge.all.round3.no_bnd.padded.id.annot.bed \
   | awk -v sample="$sample"  'BEGIN{OFS="\t"}{
      print $0, sample;
    }' \
   | cut -f -5,7-11,13,15,17-20 \
   | awk 'BEGIN{OFS="\t"}{
  hits_str=0; hits_segdup=0; hits_sat=0; has_ill=0; has_lr=0; has_ass=0; has_dip=0; has_mat=0; has_pat=0;
  ct=split($8, spl, ","); 
  for(ii=1; ii<=ct; ii++) {
    if(spl[ii]=="ill") {has_ill=1; has_dip=has_dip+1;}
    else {
      split(spl[ii], spl1, ".");      
      if(spl1[1]=="svim" || spl1[1]=="pbsv" || spl1[1]=="sniffles") has_lr=has_lr+1;
      if(spl1[1]=="sv" || spl1[1]=="pav" || spl1[1]=="paf" || spl1[1]=="svimasm") has_ass=has_ass+1;
      if(spl1[2]=="mat") has_mat=has_mat+1;
      if(spl1[2]=="pat") has_pat=has_pat+1;
      if(spl1[2]=="diploid") has_dip=has_dip+1;
    }
  }
  if($10>0 || $13>0.5) hits_str=1;
  if($11>0 || $14>0.5) hits_segdup=1;
  if($12>0 || $15>0.5) hits_sat=1;
  print $0, hits_str, hits_segdup, hits_sat, has_ill, has_lr, has_ass, has_mat, has_pat, has_dip;
}' | cut -f -9,16- | awk '{if($1!~/chr[UMXYE]/ && $1!~/random/) print $0}' \
| gzip -c > $outdir/$sample/diploid/annot/$sample.post_merge.4.annot.bed.gz
done < <(cat sample.map.haps | cut -f 1 | sort -u)
#

zjoin="/usr/bin/python2 /home/abelhj/zjoin"

while read sample
do

  zcat vcfs_clean_filtered_0811/$sample/diploid/annot/$sample.post_merge.4.annot.bed.gz | cut -f 4- \
  | sort -u \
  | awk 'BEGIN{OFS="\t"}{
   vcfid=$1; 
   svtype=$2;
   n_callers=$4;
   callers=$5;
   id=$6;
   sample=$7;
   hits_str=$8
   hits_segdup=$9;
   hits_sat=$10;
   has_ill=$11;
   has_lr=$12;
   has_ass=$13;
   has_mat=$14;
   has_pat=$15;
   has_dip=$16;
   if(has_ill>0) has_ill=1;
   if(has_lr>0) has_lr=1;
   if(has_ass>0) has_ass=1;
   if(has_mat>0) has_mat=1;
   if(has_pat>0) has_pat=1;
   if(has_dip>0) has_dip=1;
   print vcfid, "NCALLERS="n_callers";CALLERS="callers";HITS_STR="hits_str";HITS_SEGDUP="hits_segdup";HITS_SAT="hits_sat";HAS_ILL="has_ill";HAS_LR="has_lr";HAS_ASS="has_ass";HAS_MAT="has_mat";HAS_PAT="has_pat";HAS_DIP="has_dip;
   if(svtype=="BND") {
     gsub("_1", "_2", vcfid); 
     print vcfid, "NCALLERS="n_callers";CALLERS="callers";HITS_STR="hits_str";HITS_SEGDUP="hits_segdup";HITS_SAT="hits_sat";HAS_ILL="has_ill";HAS_LR="has_lr";HAS_ASS="has_ass";HAS_MAT="has_mat";HAS_PAT="has_pat";HAS_DIP="has_dip;
   }
  }' > vcfs_clean_filtered_0811/$sample/diploid/annot/$sample.temp.txt


  cat vcfs_clean_filtered_0811/$sample/diploid/$sample.merged.3.100.vcf \
  | $zjoin -p "#" -a stdin -b vcfs_clean_filtered_0811/$sample/diploid/annot/$sample.temp.txt -1 3 -2 1 \
  | awk 'BEGIN{OFS="\t"}{if($1~/#/) print $0; else {$8=$8";"$10; $9="GT"; $10="0/1"; print $0;}}' \
 > vcfs_clean_filtered_0811/$sample/diploid/$sample.merged.3.100.annot.vcf

LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q general \
-oo $outdir/logs/tovcf.%J.log \
-a 'docker(apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
"/opt/hall-lab/python-2.7.15/bin/python2  haploid_vcf_converters/revert.final.py \
-i vcfs_clean_filtered_0811/$sample/diploid/$sample.merged.3.100.annot.vcf | /opt/hall-lab/htslib*/bin/bgzip -c > vcfs_clean_filtered_0811/$sample/diploid/$sample.merged.3.100.annot.final.vcf.gz"

done < <(cat sample.map.haps | cut -f 1 | sort -u)



   

gunzip -c  *0811/all_samples.annot.2.bed.gz | cut -f -5,7-11,13,15,17-20   | awk 'BEGIN{OFS="\t"}{
  hits_str=0; hits_segdup=0; hits_sat=0; has_ill=0; has_lr=0; has_ass=0; has_dip=0; has_mat=0; has_pat=0;
  ct=split($8, spl, ","); 
  for(ii=1; ii<=ct; ii++) {
    if(spl[ii]=="ill") {has_ill=1; has_dip=has_dip+1;}
    else {
      split(spl[ii], spl1, ".");      
      if(spl1[1]=="svim" || spl1[1]=="pbsv" || spl1[1]=="sniffles") has_lr=has_lr+1;
      if(spl1[1]=="sv" || spl1[1]=="pav" || spl1[1]=="paf" || spl1[1]=="svimasm") has_ass=has_ass+1;
      if(spl1[2]=="mat") has_mat=has_mat+1;
      if(spl1[2]=="pat") has_pat=has_pat+1;
      if(spl1[2]=="diploid") has_dip=has_dip+1;
    }
  }
  if($10>0 || $13>0.5) hits_str=1;
  if($11>0 || $14>0.5) hits_segdup=1;
  if($12>0 || $15>0.5) hits_sat=1;
  print $0, hits_str, hits_segdup, hits_sat, has_ill, has_lr, has_ass, has_mat, has_pat, has_dip;
}' | cut -f -9,16- | awk '{if($1!~/chr[UMXYE]/ && $1!~/random/) print $0}' \
| cut -f 5- | gzip -c > all_samples.0811.annot.2.summ.bed.gz
