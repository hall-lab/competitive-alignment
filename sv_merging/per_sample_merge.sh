#!/bin/bash

sample=$1
outdir=$2

echo $outdir


#filter hall-lab SV calls. exclude variants where either breakpoint or >50% of outer span hits exclude region

for hap in mat pat
do
  vcf1=$outdir/$sample/$hap/sv.$sample.$hap.unfiltered.vcf
  LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub  -G compute-ccdg -q general  \
  -oo $outdir/$sample/logs/tovcf.%J.log \
   -a 'docker(halllab/vcf_bed_utils@sha256:09c18a5827d67891792ffc110627c7fa05b2262df4b91d6967ad6e544f41e8ec)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
   "bash ./sv_filter.sh $sample $hap $vcf1  $outdir"
done

#round 1 of svtools sort and merge

LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub -G compute-ccdg -K \
-oo $outdir/$sample/logs/annot.%J.log \
-q general -a  'docker(halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
"bash ./sort_merge.1.sh $sample $outdir"

#merge in illumina calls
LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub -K -G compute-ccdg \
-oo $outdir/$sample/logs/annot.%J.log \
-q general -a  'docker(halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
"bash ./sort_merge.2.sh $sample  $outdir"


LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
LSF_DOCKER_VOLUMES="${HOME}:${HOME} /scratch1/fs1/ccdg:/scratch1/fs1/ccdg /storage1/fs1/ccdg/Active:/storage1/fs1/ccdg/Active" bsub -G compute-ccdg -K \
-oo $outdir/$sample/logs/annot.%J.log -q general -a 'docker(halllab/vcf_bed_utils@sha256:09c18a5827d67891792ffc110627c7fa05b2262df4b91d6967ad6e544f41e8ec)' -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' bash -lc \
"bash ./annot_summarize.2.sh $sample $outdir "



#




