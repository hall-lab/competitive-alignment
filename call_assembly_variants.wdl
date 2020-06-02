version 1.0

workflow CallAssemblyVariants {
    input {
        String assembly_name
        File contigs1
        File contigs2
        File ref
        File ref_index
    }

    call align_contigs as align_contig1_to_ref {
        input:
            contigs=contigs1,
            ref=ref,
            assembly_name=assembly_name,
    }

    call align_contigs as align_contig2_to_ref {
        input:
            contigs=contigs2,
            ref=ref,
            assembly_name=assembly_name,
    }

    call align_contigs align_contigs_to_each_other {
        input:
            contigs=contigs1,
            ref=contigs2,
            assembly_name=assembly_name
    }

    call call_small_variants as call_small_variants1_ref {
        input:
            alignment=align_contigs1_to_ref.bam,
            ref=ref,
            assembly_name=assembly_name
    }

    call call_small_variants as call_small_variants2_ref {
        input:
            alignment=align_contigs2_to_ref.bam,
            ref=ref,
            assembly_name=assembly_name
    }
    
    call call_small_variants as call_small_variants_self {
        input:
            alignment=align_contigs_to_each_other.bam,
            ref=contigs2,
            assembly_name=assembly_name
    }

    call call_sv as call_sv1_ref {
        input:
            alignment=align_contigs1_to_ref.bam,
            contigs=contigs1,
            ref=ref,
            assembly_name=assembly_name
    }

    call call_sv as call_sv2_ref {
        input:
            alignment=align_contigs2_to_ref.bam,
            contigs=contigs2,
            ref=ref,
            assembly_name=assembly_name
    }

    call call_sv as call_sv_self {
        input:
            alignment=align_contigs_to_each_other.bam,
            contigs=contigs1,
            ref=contigs2,
            assembly_name=assembly_name
    }

    call combine_small_variants {
        input:
            small_variants1_ref = call_small_variants1_ref.vcf,
            small_variants1_ref_index = call_small_variants1_ref.vcf_index,
            small_variants2_ref = call_small_variants2_ref.vcf,
            small_variants2_ref_index = call_small_variants2_ref.vcf_index,
            small_variants_self = call_small_variants_self.vcf,
            small_variants_self_index = call_small_variants_self.vcf_index,
            contigs1 = contigs1,
            contigs2 = contigs2,
            ref = ref
    }

    call combine_sv {
        input:
            sv_ref1 = call_sv1_ref.bedpe,
            sv_ref2 = call_sv2_ref.bedpe,
            sv_self = call_sv_self.bedpe,
            contigs1 = contigs1,
            contigs2 = contigs2,
            ref = ref
    }

    output {
        small_variants = combine_small_variants.fasta
        sv = combine_sv.fasta
    }
}

task combine_sv {
    input {
        File sv_ref1
        File sv_ref2
        File sv_self
        File contigs1
        File contigs2
        File ref
    }
    command <<<
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File fasta = ~{assembly_name}.sv.combined.fasta
    }
}

task combine_small_variants {
    input {
        File small_variants1_ref
        File small_variants1_ref_index
        File small_variants2_ref
        File small_variants2_ref_index
        File small_variants_self
        File small_variants_self_index
        File contigs1
        File contigs2
        File ref
    }
    command <<<
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File fasta = ~{assembly_name}.small_variants.combined.fasta
    }
}

task call_sv {
    input {
        File alignment
        File contigs
        File ref
        String assembly_name
    }
    command <<<
        set -exo pipefail
        SAMTOOLS=/opt/hall-lab/samtools-1.9/bin/samtools
        PYTHON=/opt/hall-lab/python-2.7.15/bin/python
        SPLIT_TO_BEDPE=/opt/hall-lab/scripts/splitReadSamToBedpe
        BEDPE_TO_BKPTS=/opt/hall-lab/scripts/splitterToBreakpoint
        SVTOOLS=/opt/hall-lab/python-2.7.15/bin/svtools
        PERL=/usr/bin/perl
        REARRANGE_BREAKPOINTS=/opt/hall-lab/scripts/rearrange_breakpoints.pl
        GREP=/bin/grep
        ADD_ALIGNMENT_GAP_INFO=/opt/hall-lab/scripts/add_alignment_gap_info.pl

        $SAMTOOLS sort -n T ~{assembly_name}.tmp -O bam ~{alignment} > ~{assembly_name}.namesorted.bam
        $SAMTOOLS view -h -F 4 ~{assembly_name}.namesorted.bam | $PYTHON $SPLIT_TO_BEDPE -i stdin > ~{assembly_name}.bedpe
        $PYTHON $BEDPE_TO_BKPTS -i ~{assembly_name}.bedpe -f ~{assembly_name} -q ~{contigs} -e ~{ref} > ~{assembly_name}.breakpoints.bedpe
        $SVTOOLS bedpesort ~{assembly_name}.breakpoints.bedpe | $PERL $REARRANGE_BREAKPOINTS > ~{assembly_name}.breakpoints.sorted.bedpe
        cat <($GREP "^#" ~{assembly_name}.breakpoints.sorted.bedpe) <(paste <($GREP -v "^#" ~{assembly_name}.breakpoints.sorted.bedpe | cut -f 1-6) <(paste -d : <($GREP -v "^#" ~{assembly_name}.breakpoints.sorted.bedpe | cut -f 7) <($GREP -v "^#" ~{assembly_name}.breakpoints.sorted.bedpe | cut -f 19 | sed 's/.*SVLEN=/SVLEN=/' | sed 's/;.*//')) <($GREP -v "^#" ~{assembly_name}.breakpoints.sorted.bedpe | cut -f 8-)) | $PERL $ADD_ALIGNMENT_GAP_INFO > ~{assembly_name}.breakpoints.sorted.fixed.bedpe
    mv ~{assembly_name}.breakpoints.sorted.fixed.bedpe ~{assembly_name}.breakpoints.sorted.bedpe
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File bedpe = ~{assembly_name}.breakpoints.sorted.bedpe
    }
}

task call_small_variants {
    input {
        File alignment
        File ref
        String assembly_name
    }
    command <<<
        set -exo pipefail
        SAMTOOLS=/opt/hall-lab/samtools-1.9/bin/samtools
        PAFTOOLS=/opt/hall-lab/minimap2/misc/paftools.js
        K8=/opt/hall-lab/minimap2/k8
        BGZIP=/opt/hall-lab/htslib-1.9/bin/bgzip
        PYTHON=/opt/hall-lab/python-2.7.15/bin/python
        VAR_TO_VCF=/opt/hall-lab/scripts/varToVcf.py
        SVTOOLS=/opt/hall-lab/python-2.7.15/bin/svtools
        PERL=/usr/bin/perl
        GENOTYPE_VCF=/opt/hall-lab/scripts/vcfToGenotyped.pl
        TABIX=/opt/hall-lab/htslib-1.9/bin/tabix
        $SAMTOOLS view -h ~{alignment} | $K8 $PAFTOOLS sam2paf - | sort -k6,6 -k8,8n | $K8 $PAFTOOLS call -l 1 -L 1 -q 0 - | grep "^V" | sort -V | $BGZIP -c > ~{assembly_name}.loose.vcf.txt.gz
        $PYTHON $VAR_TO_VCF -i <(zcat ~{assembly_name}.loose.var.txt.gz) -r ~{ref} -s ~{assembly_name} -o ~{assembly_name}.loose.vcf
        $SVTOOLS vcfsort ~{assembly_name}.loose.vcf | $PERL $GENOTYPE_VCF | $BGZIP -c > ~{assembly_name}.loose.genotyped.vcf.gz
        $TABIX -f -p vcf $OUTPUT_PREFIX.loose.genotyped.vcf.gz
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File vcf = ~{assembly_name}.loose.genotyped.vcf.gz
        File vcf_index = ~{assembly_name}.loose.genotyped.vcf.gz.tbi
    }
}

task align_contigs {
    input {
            File contigs
            File ref
            String assembly_name
    }
    command <<<
        set -exo pipefail
        MINIMAP2=/opt/hall-lab/minimap2/minimap2
        SAMTOOLS=/opt/hall-lab/samtools-1.9/bin/samtools
        $MINIMAP2 -ax asm5 -L --cs ~{ref} ~{contigs} | $SAMTOOLS sort -T ~{assembly_name}.tmp -O bam - > ~{assembly_name}.bam
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File bam = ~{assembly_name}.bam
    }
}
