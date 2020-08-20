version 1.0
import "convert_to_fasta.wdl" as convert_to_fasta
import "get_read_support.wdl" as get_read_support

workflow CallAssemblyVariants {
    input {
        String assembly_name
        File contigs1
        File contigs2
        File ref
        File ref_index
        String ref_name
        File fastq_list
    }

    call align_contigs as align_contig1_to_ref {
        input:
            contigs=contigs1,
            ref=ref
    }

    call align_contigs as align_contig2_to_ref {
        input:
            contigs=contigs2,
            ref=ref
    }

    call align_contigs as align_contigs_to_each_other {
        input:
            contigs=contigs1,
            ref=contigs2
    }

    call index_fasta as index_contigs1 {
        input:
            fasta=contigs1
    }

    call index_fasta as index_contigs2 {
        input:
            fasta=contigs2
    }

    call call_small_variants as call_small_variants1_ref {
        input:
            alignment=align_contig1_to_ref.bam,
            ref=ref,
            ref_index=ref_index,
            assembly_name=assembly_name,
            id_prefix="ref1"
    }

    call call_small_variants as call_small_variants2_ref {
        input:
            alignment=align_contig2_to_ref.bam,
            ref=ref,
            ref_index=ref_index,
            assembly_name=assembly_name,
            id_prefix="ref2"
    }

    call call_small_variants as call_small_variants_self {
        input:
            alignment=align_contigs_to_each_other.bam,
            ref=index_contigs2.unzipped_fasta,
            ref_index=index_contigs2.fasta_index,
            assembly_name=assembly_name,
            id_prefix="self"
    }

    call call_sv as call_sv1_ref {
        input:
            alignment=align_contig1_to_ref.bam,
            contigs=index_contigs1.unzipped_fasta,
            contigs_index=index_contigs1.fasta_index,
            ref=ref,
            ref_index=ref_index,
            assembly_name=assembly_name
    }

    call call_sv as call_sv2_ref {
        input:
            alignment=align_contig2_to_ref.bam,
            contigs=index_contigs2.unzipped_fasta,
            contigs_index=index_contigs2.fasta_index,
            ref=ref,
            ref_index=ref_index,
            assembly_name=assembly_name
    }

    call call_sv as call_sv_self {
        input:
            alignment=align_contigs_to_each_other.bam,
            contigs=index_contigs1.unzipped_fasta,
            contigs_index=index_contigs1.fasta_index,
            ref=index_contigs2.unzipped_fasta,
            ref_index=index_contigs2.fasta_index,
            assembly_name=assembly_name
    }

    call get_read_support.GetReadSupport as get_read_support{
        input:
            ref=ref,
            contigs1=contigs1,
            contigs2=contigs2,
            reads_list=fastq_list,
            small_variants_ref_contigs1=call_small_variants1_ref.vcf,
            small_variants_contigs1_ref=call_small_variants1_ref.inverse_vcf,
            small_variants_ref_contigs2=call_small_variants2_ref.vcf,
            small_variants_contigs2_ref=call_small_variants2_ref.inverse_vcf,
            small_variants_contigs1_contigs2=call_small_variants_self.vcf,
            small_variants_contigs2_contigs1=call_small_variants_self.inverse_vcf
    }

    call combine_small_variants_vcf {
        input:
            small_variants1_ref = call_small_variants1_ref.vcf,
            small_variants1_ref_by_query = call_small_variants1_ref.inverse_vcf,
            small_variants2_ref = call_small_variants2_ref.vcf,
            small_variants2_ref_by_query = call_small_variants2_ref.inverse_vcf,
            small_variants_self = call_small_variants_self.vcf,
            small_variants_self_by_query = call_small_variants_self.inverse_vcf,
            small_variants1_ref_index = call_small_variants1_ref.vcf_index,
            small_variants1_ref_by_query_index = call_small_variants1_ref.inverse_vcf_index,
            small_variants2_ref_index = call_small_variants2_ref.vcf_index,
            small_variants2_ref_by_query_index = call_small_variants2_ref.inverse_vcf_index,
            small_variants_self_index = call_small_variants_self.vcf_index,
            small_variants_self_by_query_index = call_small_variants_self.inverse_vcf_index
    }

    call convert_to_fasta.ConvertToFasta as convert_ref {
        input:
            vcf=combine_small_variants_vcf.combined_vcf_ref,
            vcf_index=combine_small_variants_vcf.combined_vcf_ref_index,
            ref=ref,
            ref_index=ref_index,
            ref_name=ref_name
    }

    call convert_to_fasta.ConvertToFasta as convert_self {
        input:
            vcf=combine_small_variants_vcf.remaining_vcf_contigs,
            vcf_index=combine_small_variants_vcf.remaining_vcf_contigs_index,
            ref=index_contigs2.unzipped_fasta,
            ref_index=index_contigs2.fasta_index,
            ref_name=assembly_name
    }

    call combine_small_variants {
        input:
            small_variants_ref = convert_ref.marker_fasta,
            small_variants_self = convert_self.marker_fasta,
            small_variants_ref_marker_positions = convert_ref.marker_positions,
            small_variants_self_marker_positions = convert_self.marker_positions
    }

    #call combine_sv {
    #    input:
    #        sv_ref1 = call_sv1_ref.bedpe,
    #        sv_ref2 = call_sv2_ref.bedpe,
    #        sv_self = call_sv_self.bedpe,
    #        contigs1 = contigs1,
    #        contigs2 = contigs2,
    #        ref = ref
    #}

    output {
        File small_variants = combine_small_variants.fasta
        File small_variants_marker_positions = combine_small_variants.marker_positions
        File small_variant_support_ref_contigs1 = get_read_support.small_variant_support_ref_contigs1
        File small_variant_support_ref_contigs2 = get_read_support.small_variant_support_ref_contigs2
        File small_variant_support_contigs1_contigs2 = get_read_support.small_variant_support_contigs1_contigs2
    #   File sv = combine_sv.fasta
    }
}

task index_fasta {
    input {
        File fasta
    }
    command <<<
        set -exo pipefail
        SAMTOOLS=/opt/hall-lab/samtools-1.9/bin/samtools
        zcat ~{fasta} > unzipped.fa
        $SAMTOOLS faidx unzipped.fa
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/analyze_assemblies@sha256:f1e125e1261e163ef790ff11f3f2749f71427bf8f71ca95dd1894ce7c6c804eb"
    }
    output {
        File unzipped_fasta = "unzipped.fa"
        File fasta_index = "unzipped.fa.fai"
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
    #TODO
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:f1e125e1261e163ef790ff11f3f2749f71427bf8f71ca95dd1894ce7c6c804eb"
    }
    output {
        File fasta = "sv.combined.fasta"
    }
}

task combine_small_variants_vcf {
    input {
        File small_variants1_ref
        File small_variants1_ref_by_query
        File small_variants2_ref
        File small_variants2_ref_by_query
        File small_variants_self
        File small_variants_self_by_query
        File small_variants1_ref_index
        File small_variants1_ref_by_query_index
        File small_variants2_ref_index
        File small_variants2_ref_by_query_index
        File small_variants_self_index
        File small_variants_self_by_query_index
    }
    command <<<
        set -exo pipefail
        BCFTOOLS=/opt/hall-lab/bcftools-1.9/bin/bcftools
        TABIX=/opt/hall-lab/htslib-1.9/bin/tabix
        mkdir tmp1 tmp2 tmp3 tmp4
        $BCFTOOLS isec -p tmp1 ~{small_variants_self} ~{small_variants1_ref_by_query} && $BCFTOOLS query -f '%ID\n' tmp1/0002.vcf > self1_in_ref
        $BCFTOOLS isec -p tmp2 ~{small_variants_self_by_query} ~{small_variants2_ref_by_query} && $BCFTOOLS query -f '%ID\n' tmp2/0002.vcf > self2_in_ref
        cat self1_in_ref self2_in_ref > exclude_ids_ref
        $BCFTOOLS isec -p tmp3 ~{small_variants1_ref_by_query} ~{small_variants_self} && $BCFTOOLS query -f '%ID\n' tmp3/0002.vcf > ref_in_self1
        $BCFTOOLS isec -p tmp4 ~{small_variants2_ref_by_query} ~{small_variants_self_by_query} && $BCFTOOLS query -f '%ID' tmp4/0002.vcf > ref_in_self2
        cat ref_in_self1 ref_in_self2 > exclude_ids_self

        $BCFTOOLS concat -a -d all ~{small_variants1_ref} ~{small_variants2_ref} | $BCFTOOLS view -i '%ID!=@exclude_ids_ref' -o small_variants.combined.vcf.gz -O z
        $BCFTOOLS view -o small_variants.contigs.vcf.gz -O z -i '%ID!=@exclude_ids_self' ~{small_variants_self}
        $TABIX -f -p vcf small_variants.combined.vcf.gz
        $TABIX -f -p vcf small_variants.contigs.vcf.gz
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/bcftools_samtools@sha256:49dc9b0b9419281b87df4a9faf9ca6681725317108bca66ba37f9bd1d86e9de2"
    }
    output {
        File combined_vcf_ref = "small_variants.combined.vcf.gz"
        File remaining_vcf_contigs = "small_variants.contigs.vcf.gz"
        File combined_vcf_ref_index = "small_variants.combined.vcf.gz.tbi"
        File remaining_vcf_contigs_index = "small_variants.contigs.vcf.gz.tbi"
    }
}

task combine_small_variants {
    input {
        File small_variants_ref
        File small_variants_self
        File small_variants_ref_marker_positions
        File small_variants_self_marker_positions
    }
    command <<<
        set -exo pipefail
        PYTHON=/opt/hall-lab/python-2.7.15/bin/python
        FIND_DUPS=/storage1/fs1/ccdg/Active/analysis/ref_grant/assembly_analysis_20200220/multiple_competitive_alignment/find_duplicate_markers.py #TODO
        #combine fasta files and sort by sequence
        cat ~{small_variants_ref} ~{small_variants_self} | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep -v "^$" | paste - - - - | awk -v OFS="\t" -v FS="\t" '{if($2<$4) {print($2, $4, $1, $3)} else{print($4,$2,$3,$1)}}' | sort | awk -v OFS="\n" -v FS="\t" '{print($3,$1,$4,$2)}' > tmp
        #find duplicate markers
        $PYTHON $FIND_DUPS -i tmp > small_variants.combined.fasta
        cat ~{small_variants_ref_marker_positions} ~{small_variants_self_marker_positions} | sort -u > small_variants.marker_positions.txt
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:5cbac56b15b739783c37d2a92261bef138d5bae3e99171557df06d3e39cb485a"
    }
    output {
        File fasta = "small_variants.combined.fasta"
        File marker_positions = "small_variants.marker_positions.txt"
    }
}

task call_sv {
    input {
        File alignment
        File contigs
        File contigs_index
        File ref
        File ref_index
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

        ln -s ~{contigs} contigs.fa
        ln -s ~{contigs_index} contigs.fa.fai
        ln -s ~{ref} ref.fa
        ln -s ~{ref_index} ref.fa.fai
        $SAMTOOLS sort -n -T tmp -O bam ~{alignment} > namesorted.bam
        $SAMTOOLS view -h -F 4 namesorted.bam | $PYTHON $SPLIT_TO_BEDPE -i stdin > split.bedpe
        $PYTHON $BEDPE_TO_BKPTS -i split.bedpe -f ~{assembly_name} -q contigs.fa -e ref.fa > breakpoints.bedpe
        $SVTOOLS bedpesort breakpoints.bedpe | $PERL $REARRANGE_BREAKPOINTS > breakpoints.sorted.bedpe
        cat <($GREP "^#" breakpoints.sorted.bedpe) <(paste <($GREP -v "^#" breakpoints.sorted.bedpe | cut -f 1-6) <(paste -d : <($GREP -v "^#" breakpoints.sorted.bedpe | cut -f 7) <($GREP -v "^#" breakpoints.sorted.bedpe | cut -f 19 | sed 's/.*SVLEN=/SVLEN=/' | sed 's/;.*//')) <($GREP -v "^#" breakpoints.sorted.bedpe | cut -f 8-)) | $PERL $ADD_ALIGNMENT_GAP_INFO > breakpoints.sorted.fixed.bedpe
    mv breakpoints.sorted.fixed.bedpe breakpoints.sorted.bedpe
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:f1e125e1261e163ef790ff11f3f2749f71427bf8f71ca95dd1894ce7c6c804eb"
    }
    output {
        File bedpe = "breakpoints.sorted.bedpe"
    }
}

task call_small_variants {
    input {
        File alignment
        File ref
        File ref_index
        String assembly_name
        String id_prefix
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
        ln -s ~{ref} ref.fa
        ln -s ~{ref_index} ref.fa.fai
        $SAMTOOLS view -h ~{alignment} | $K8 $PAFTOOLS sam2paf - | sort -k6,6 -k8,8n | $K8 $PAFTOOLS call -l 1 -L 1 -q 0 - | grep "^V" | sort -V | $BGZIP -c > loose.var.txt.gz
        $PYTHON $VAR_TO_VCF -i <(zcat loose.var.txt.gz) -r ref.fa -s ~{assembly_name} -o loose.vcf -p ~{id_prefix}
        $SVTOOLS vcfsort loose.vcf | $PERL $GENOTYPE_VCF | $BGZIP -c > loose.genotyped.vcf.gz
        $TABIX -f -p vcf loose.genotyped.vcf.gz
        $SVTOOLS vcfsort loose.vcf.2.vcf | $PERL $GENOTYPE_VCF | $BGZIP -c > loose2.genotyped.vcf.gz
        $TABIX -f -p vcf loose2.genotyped.vcf.gz
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:5cbac56b15b739783c37d2a92261bef138d5bae3e99171557df06d3e39cb485a"
    }
    output {
        File vcf = "loose.genotyped.vcf.gz"
        File vcf_index = "loose.genotyped.vcf.gz.tbi"
        File inverse_vcf = "loose2.genotyped.vcf.gz"
        File inverse_vcf_index = "loose2.genotyped.vcf.gz.tbi"
    }
}

task align_contigs {
    input {
            File contigs
            File ref
    }
    command <<<
        set -exo pipefail
        MINIMAP2=/opt/hall-lab/minimap2/minimap2
        SAMTOOLS=/opt/hall-lab/samtools-1.9/bin/samtools
        $MINIMAP2 -ax asm5 -L --cs ~{ref} ~{contigs} | $SAMTOOLS sort -T tmp -O bam - > aligned.bam
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:f1e125e1261e163ef790ff11f3f2749f71427bf8f71ca95dd1894ce7c6c804eb"
    }
    output {
        File bam = "aligned.bam"
    }
}
