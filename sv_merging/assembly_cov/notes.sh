#!/bin/bash

mkdir logs


## get contig ids from each assembly
while read sample mat pat
do
  mkdir -p  $sample/logs
  LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
  LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub -G compute-ccdg \
  -oo logs/$sample.%J.log \
  -q ccdg -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
   "bash ./get_contigs.sh $sample $mat $pat"
done < sample.map.v2
#


#copy contig to reference alignments out of cromwell dirs


while read sample shard
do
  pth=/storage1/fs1/ccdg/Active/analysis/abelhj/pangenome/cromwell/ca_round2/cromwell-executions/CompetitiveGenotyping/87ef843c-5739-42ad-9c09-95332e1226c0/call-call_variants/shard-$shard/CallAssemblyVariants
  bam1=$pth/*/call-align_contig1_to_ref/execution/aligned.bam
  bam2=$pth/*/call-align_contig2_to_ref/execution/aligned.bam
  LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q ccdg \
  -oo logs/svtools.%J.log -M 16000000 -R'select[mem>16000] rusage[mem=16000]' \
  -a 'docker(apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca)' \
  "cp  $bam2 $sample/$sample.2.bam"
  LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q ccdg \
  -oo logs/svtools.%J.log -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' \
  -a 'docker(apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca)' \
  "cp  $bam2.bai $sample/$sample.2.bam.bai"
LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q ccdg \
  -oo logs/svtools.%J.log -M 16000000 -R'select[mem>16000] rusage[mem=16000]' \
  -a 'docker(apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca)' \
  "cp  $bam1 $sample/$sample.1.bam"
  LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q ccdg \
  -oo logs/svtools.%J.log -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' \
  -a 'docker(apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca)' \
  "cp  $bam1.bai $sample/$sample.1.bam.bai"

done < <(cat ../sample.map.new | cut -f -2)


#split bams by contig


while read sample
do
 bam1=$sample/$sample.1.bam
 bam2=$sample/$sample.2.bam
 cont1=$sample/$sample.mat.contigs.txt
 cont2=$sample/$sample.pat.contigs.txt

 mkdir -p $sample/by_contig/mat $sample/by_contig/pat
 while read contig
 do
 echo "./subset.sh $bam1  $contig $sample/by_contig/mat"
 done < $cont1

 while read contig
 do
   echo "./subset.sh $bam2  $contig $sample/by_contig/pat"
 done < $cont2
 
done < <(cat sample.map.v2 | cut -f 1) > tosubset.sh



cat tosubset.sh | awk 'BEGIN{OFS="\t"}{print $0, int(NR/50);}' > subsets.batched.txt

mkdir batches

while read batch
do
cat subsets.batched.txt \
| awk -v batch="$batch" 'BEGIN{FS="\t"}{if($2==batch) print $1;}' > batches/batch$batch.sh

chmod 755  batches/batch$batch.sh
LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q ccdg \
 -oo logs/svtools.%J.log -M 16000000 -R'select[mem>16000] rusage[mem=16000]' \
 -g /abelhj/subset2/$sample.2 \
 -a 'docker(apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca)' bash -lc \
  "bash batches/batch$batch.sh"

done < <(cat subsets.batched.txt | cut -f 2 | sort -u) 



#calculate coverage depth of each contig

mkdir batches/cov

while read sample
do
  bam1=$sample/$sample.1.bam
  bam2=$sample/$sample.2.bam
  cont1=$sample/$sample.mat.contigs.txt
  cont2=$sample/$sample.pat.contigs.txt

  mkdir -p $sample/by_contig/mat/cov $sample/by_contig/pat/cov
  base1=`basename $bam1 .bam`
  base2=`basename $bam2 .bam`
  while read contig
  do
    inbam=$sample/by_contig/mat/$base1.$contig.bam
    outcov=$sample/by_contig/mat/cov/$base1.$contig.bedg
    printf "/opt/hall-lab/bin/bedtools genomecov  -bg -ibam %s | sort -k1,1V -k2,3n > %s\n" $inbam $outcov
  done < $cont1 > batches/cov/$sample.mat.sh
  
  
  while read contig
  do
    inbam=$sample/by_contig/pat/$base2.$contig.bam
    outcov=$sample/by_contig/pat/cov/$base2.$contig.bedg
    printf "/opt/hall-lab/bin/bedtools genomecov  -bg -ibam %s | sort -k1,1V -k2,3n > %s\n" $inbam $outcov
  done < $cont2 > batches/cov/$sample.pat.sh
  
  chmod 755 batches/cov/$sample.pat.sh batches/cov/$sample.mat.sh
  
  LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q ccdg \
  -a 'docker(halllab/vcf_bed_utils@sha256:09c18a5827d67891792ffc110627c7fa05b2262df4b91d6967ad6e544f41e8ec)' \
  -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' \
  -oo logs/cov.$sample.%J.log bash -lc \
  "bash batches/cov/$sample.pat.sh"
  
  LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q ccdg \
  -a 'docker(halllab/vcf_bed_utils@sha256:09c18a5827d67891792ffc110627c7fa05b2262df4b91d6967ad6e544f41e8ec)' \
  -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' \
  -oo logs/cov.$sample.%J.log bash -lc \
  "bash batches/cov/$sample.mat.sh"


done < <(cat sample.map.v2 | cut -f 1 ) 


# calculate genome coverage by all contigs
# exclude file is all regions with >4x total coverage or covered by more than one contig

while read sample 
do
 str=`ls $sample/by_contig/mat/cov/*.bedg `
 /opt/hall-lab/bin/bedtools unionbedg -i $str | /opt/hall-lab/htslib-1.9/bin/bgzip -c > $sample/by_contig/mat/$sample.mat.union.bedg.gz

  str=`ls $sample/by_contig/pat/cov/*.bedg `
  /opt/hall-lab/bin/bedtools unionbedg -i $str | /opt/hall-lab/htslib-1.9/bin/bgzip -c > $sample/by_contig/pat/$sample.pat.union.bedg.gz
  for pp in mat pat
  do
    zcat $sample/by_contig/$pp/$sample.$pp.union.bedg.gz \
    | awk 'BEGIN{OFS="\t"}{
      for(jj=4; jj<=NF; jj++) {
         gt0=0; 
         if($jj>0) gt0=1; 
         if (gt0==1) print $1, $2, $3, $jj, gt0;
       }
     }' | sort -k 1,1V -k2,3n \
     | /opt/hall-lab/bin/bedtools groupby -g 1,2,3 -c 4,5 -o sum,sum | gzip -c > $sample/by_contig/$pp/$sample.$pp.summary.bedg.gz
   done

  for pp in mat pat
   do
     zcat $sample/by_contig/$pp/$sample.$pp.summary.bedg.gz \
     | awk '{if($4>3 || $5>1) print $0;}' \
     | sort -k1,1V -k2,3n | cut -f -3 \
     | /opt/hall-lab/bin/bedtools merge -i stdin -d 10000 > $sample/$sample.$pp.exclude.bed
   done
done < <(cat sample.map.new | cut -f 1)

