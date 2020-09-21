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
        File segdup_bed
        File str_bed
        #File fastq_list
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

    #call get_read_support.GetReadSupport as get_read_support{
    #    input:
    #        ref=ref,
    #        contigs1=contigs1,
    #        contigs2=contigs2,
    #        reads_list=fastq_list,
    #        small_variants_ref_contigs1=call_small_variants1_ref.vcf,
    #        small_variants_contigs1_ref=call_small_variants1_ref.inverse_vcf,
    #        small_variants_ref_contigs2=call_small_variants2_ref.vcf,
    #        small_variants_contigs2_ref=call_small_variants2_ref.inverse_vcf,
    #        small_variants_contigs1_contigs2=call_small_variants_self.vcf,
    #        small_variants_contigs2_contigs1=call_small_variants_self.inverse_vcf
    #}

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

    #call convert_to_fasta.ConvertToFasta as convert_ref {
    #    input:
    #        vcf=combine_small_variants_vcf.combined_vcf_ref,
    #        vcf_index=combine_small_variants_vcf.combined_vcf_ref_index,
    #        ref=ref,
    #        ref_index=ref_index,
    #        ref_name=ref_name
    #}

    #call convert_to_fasta.ConvertToFasta as convert_self {
    #    input:
    #        vcf=combine_small_variants_vcf.remaining_vcf_contigs,
    #        vcf_index=combine_small_variants_vcf.remaining_vcf_contigs_index,
    #        ref=index_contigs2.unzipped_fasta,
    #        ref_index=index_contigs2.fasta_index,
    #        ref_name=assembly_name
    #}

    #call combine_small_variants {
    #    input:
    #        small_variants_ref = convert_ref.marker_fasta,
    #        small_variants_self = convert_self.marker_fasta,
    #        small_variants_ref_marker_positions = convert_ref.marker_positions,
    #        small_variants_self_marker_positions = convert_self.marker_positions
    #}

    call combine_sv {
        input:
            sv_ref1 = call_sv1_ref.bedpe,
            sv_ref2 = call_sv2_ref.bedpe,
    }

    call count_variants {
        input:
            small_variants_vcf = combine_small_variants_vcf.combined_vcf_ref,
            sv_bedpe = combine_sv.bedpe,
            segdup_bed = segdup_bed,
            str_bed = str_bed,
            assembly_name = assembly_name
    }

    call count_self_variants {
        input:
            novel_vcf_self = combine_small_variants_vcf.novel_vcf_self,
            known_vcf_self = combine_small_variants_vcf.known_vcf_self,
            unique_vcf_ref = combine_small_variants_vcf.unique_vcf_ref,
            nonunique_vcf_ref = combine_small_variants_vcf.nonunique_vcf_ref,
            novel_vcf_self_index = combine_small_variants_vcf.novel_vcf_self_index,
            known_vcf_self_index = combine_small_variants_vcf.known_vcf_self_index,
            unique_vcf_ref_index = combine_small_variants_vcf.unique_vcf_ref_index,
            nonunique_vcf_ref_index = combine_small_variants_vcf.nonunique_vcf_ref_index,
            segdup_bed = segdup_bed,
            str_bed = str_bed,
            assembly_name = assembly_name
    }

    output {
        File sv_ref1 = call_sv1_ref.bedpe
        File sv_ref2 = call_sv2_ref.bedpe
        File sv_self = call_sv_self.bedpe
        File small_variants_ref1 = call_small_variants1_ref.vcf
        File small_variants_ref2 = call_small_variants2_ref.vcf
        File small_variants_self = call_small_variants_self.vcf
        File small_variants_ref_combined = combine_small_variants_vcf.combined_vcf_ref
        File small_variants_ref_combined_index = combine_small_variants_vcf.combined_vcf_ref_index
        File sv_combined = combine_sv.bedpe
        File counts = count_variants.counts
        File self_counts = count_self_variants.counts
    #    File small_variants_marker_positions = combine_small_variants.marker_positions
    #    File small_variant_support_ref_contigs1 = get_read_support.small_variant_support_ref_contigs1
    #    File small_variant_support_ref_contigs2 = get_read_support.small_variant_support_ref_contigs2
    #    File small_variant_support_contigs1_contigs2 = get_read_support.small_variant_support_contigs1_contigs2
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
        docker: "apregier/analyze_assemblies@sha256:54669591da03e517f61097f93f8eac512368ae503954276b0149b13ebae0aec4"
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
    }
    command <<<
        set -exo pipefail
        BEDTOOLS=/opt/hall-lab/bedtools
        SVTOOLS=/opt/hall-lab/python-2.7.15/bin/svtools
        $BEDTOOLS pairtopair -type both -a ~{sv_ref1} -b ~{sv_ref2} -is -slop 50 | $SVTOOLS bedpesort | cut -f 1-6,11 | uniq | awk '{print $s ".homalt." NR}' > ref1_ref2_homozygous.bedpe
        $BEDTOOLS pairtopair -type notboth -a ~{sv_ref1} -b ~{sv_ref2} -is -slop 50 | $SVTOOLS bedpesort | cut -f 1-6,11 | uniq | awk '{print $s ".ref1." NR}' > ref1_het.bedpe
        $BEDTOOLS pairtopair -type notboth -b ~{sv_ref1} -a ~{sv_ref2} -is -slop 50 | $SVTOOLS bedpesort | cut -f 1-6,11 | uniq | awk '{print $s ".ref2." NR}' > ref2_het.bedpe
        $SVTOOLS bedpesort <(cat ref1_ref2_homozygous.bedpe ref1_het.bedpe ref2_het.bedpe) > ref1_ref2.bedpe
    >>>
    runtime {
        memory: "16G"
        docker: "apregier/analyze_assemblies@sha256:54669591da03e517f61097f93f8eac512368ae503954276b0149b13ebae0aec4"
    }
    output {
        File bedpe = "ref1_ref2.bedpe"
    }
}

task count_self_variants {
    input {
        File novel_vcf_self
        File known_vcf_self
        File unique_vcf_ref
        File nonunique_vcf_ref
        File novel_vcf_self_index
        File known_vcf_self_index
        File unique_vcf_ref_index
        File nonunique_vcf_ref_index
        File segdup_bed
        File str_bed
        String assembly_name
    }
    command <<<
        BEDTOOLS=/opt/hall-lab/bedtools
        GREP=/bin/grep
        STR_BED=str.bed
        SEGDUP_BED=segdup.bed
        cat ~{str_bed} | grep -v "^track" | sed 's/^/chr/' > $STR_BED
        cat ~{segdup_bed} | grep -v "^track" | sed 's/^/chr/' > $SEGDUP_BED
        rm -f counts.txt
        touch counts.txt
        zcat ~{novel_vcf_self} | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/self\tunique\tall\t/' >> counts.txt
        zcat ~{known_vcf_self} | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/self\tnonunique\tall\t/' >> counts.txt
        zcat ~{nonunique_vcf_ref} | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/ref\tnonunique\tall\t/' >> counts.txt
        zcat ~{unique_vcf_ref} | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/ref\tunique\tall\t/' >> counts.txt
        $BEDTOOLS intersect -v -a <(cat <(zcat ~{nonunique_vcf_ref} | $GREP "^#") <($BEDTOOLS intersect -v -a ~{nonunique_vcf_ref} -b $STR_BED)) -b $SEGDUP_BED | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/ref\tnonunique\tnonRep\t/' >> counts.txt
        $BEDTOOLS intersect -v -a <(cat <(zcat ~{unique_vcf_ref} | $GREP "^#") <($BEDTOOLS intersect -v -a ~{unique_vcf_ref} -b $STR_BED)) -b $SEGDUP_BED | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/ref\tunique\tnonRep\t/' >> counts.txt
        cat counts.txt | sed 's/^/~{assembly_name}\t/' > tmp
        mv tmp counts.txt
    >>>
    runtime {
        memory: "16G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File counts = "counts.txt"
    }
}

task count_variants {
    input {
        File small_variants_vcf
        File sv_bedpe
        File segdup_bed
        File str_bed
        String assembly_name
    }
    command <<<
        set -exo pipefail
        BEDTOOLS=/opt/hall-lab/bedtools
        PYTHON=/opt/hall-lab/python-2.7.15/bin/python
        VCFTOBEDPE=/opt/hall-lab/scripts/vcfToBedpe.py
        SVLENGTHS=/opt/hall-lab/scripts/sv_lengths.py
        GREP=/bin/grep
        STR_BED=str.bed
        SEGDUP_BED=segdup.bed
        cat ~{str_bed} | grep -v "^track" | sed 's/^/chr/' > $STR_BED
        cat ~{segdup_bed} | grep -v "^track" | sed 's/^/chr/' > $SEGDUP_BED
        $BEDTOOLS pairtobed -a ~{sv_bedpe} -b $STR_BED -type either > sv.str.bedpe
        $BEDTOOLS pairtobed -a ~{sv_bedpe} -b $SEGDUP_BED -type either > sv.allSegDup.bedpe
        $BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a ~{sv_bedpe} -b $STR_BED -type neither) -b $SEGDUP_BED -type either > sv.segDup.bedpe
        $BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a ~{sv_bedpe} -b $STR_BED -type neither) -b $SEGDUP_BED -type neither > sv.nonRep.bedpe
        rm -f counts.txt
        cat ~{sv_bedpe} | cut -f 7 | cut -f 1-2 -d . | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/sv\tall\t/' | sed 's/\./\t/' >> counts.txt
        cat sv.str.bedpe | cut -f 1-7 | uniq | cut -f 7 | cut -f 1-2 -d . | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/sv\tstr\t/' | sed 's/\./\t/' >> counts.txt
        cat sv.segDup.bedpe | cut -f 1-7 | uniq | cut -f 7 | cut -f 1-2 -d . | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/sv\tsegDup\t/' | sed 's/\./\t/' >> counts.txt
        cat sv.nonRep.bedpe | cut -f 7 | cut -f 1-2 -d . | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/sv\tnonRep\t/' | sed 's/\./\t/' >> counts.txt
        $PYTHON $VCFTOBEDPE -i ~{small_variants_vcf} -o indels.bedpe -m 1 -M 49
        #$PYTHON $SVLENTHS -v ~{small_variants_vcf} -o indel_lengths.txt
        $PYTHON $VCFTOBEDPE -i ~{small_variants_vcf} -o large_indels.bedpe -m 50
        $BEDTOOLS pairtobed -a large_indels.bedpe -b $STR_BED -type either > large_indels.str.bedpe
        $BEDTOOLS pairtobed -a large_indels.bedpe -b $SEGDUP_BED -type either > large_indels.allSegDup.bedpe
        $BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a large_indels.bedpe -b $STR_BED -type neither) -b $SEGDUP_BED -type either > large_indels.segDup.bedpe
        $BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a large_indels.bedpe -b $STR_BED -type neither) -b $SEGDUP_BED -type neither > large_indels.nonRep.bedpe
        paste <(grep -v GENOTYPE large_indels.bedpe | cut -f 8) <(grep -v GENOTYPE large_indels.bedpe | cut -f 7) | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/indels\tall\t/' >> counts.txt
        paste <(grep -v GENOTYPE large_indels.str.bedpe | cut -f 1-8 | uniq | cut -f 8) <(grep -v GENOTYPE large_indels.str.bedpe | cut -f 1-8 | uniq | cut -f 7) | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/indels\tstr\t/' >> counts.txt
        paste <(grep -v GENOTYPE large_indels.segDup.bedpe | cut -f 1-8 | uniq | cut -f 8) <(grep -v GENOTYPE large_indels.segDup.bedpe | cut -f 1-8 | uniq | cut -f 7) | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/indels\tsegDup\t/' >> counts.txt
        paste <(grep -v GENOTYPE large_indels.nonRep.bedpe | cut -f 8) <(grep -v GENOTYPE large_indels.nonRep.bedpe | cut -f 7) | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/indels\tnonRep\t/' >> counts.txt
        $BEDTOOLS pairtobed -a indels.bedpe -b $STR_BED -type either > indels.str.bedpe
        $BEDTOOLS pairtobed -a indels.bedpe -b $SEGDUP_BED -type either > indels.allSegDup.bedpe
        $BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a indels.bedpe -b $STR_BED -type neither) -b $SEGDUP_BED -type either > indels.segDup.bedpe
        $BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a indels.bedpe -b $STR_BED -type neither) -b $SEGDUP_BED -type neither > indels.nonRep.bedpe
        paste <(grep -v GENOTYPE indels.bedpe | cut -f 8) <(grep -v GENOTYPE indels.bedpe | cut -f 7) | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/small_indels\tall\t/' >> counts.txt
        paste <(grep -v GENOTYPE indels.str.bedpe | cut -f 1-8 | uniq | cut -f 8) <(grep -v GENOTYPE indels.str.bedpe | cut -f 1-8 | uniq | cut -f 7) | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/small_indels\tstr\t/' >> counts.txt
        paste <(grep -v GENOTYPE indels.segDup.bedpe | cut -f 1-8 | uniq | cut -f 8) <(grep -v GENOTYPE indels.segDup.bedpe | cut -f 1-8 | uniq | cut -f 7) | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/small_indels\tsegDup\t/' >> counts.txt
        paste <(grep -v GENOTYPE indels.nonRep.bedpe | cut -f 1-8 | uniq | cut -f 8) <(grep -v GENOTYPE indels.nonRep.bedpe | cut -f 1-8 | uniq | cut -f 7) | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/small_indels\tnonRep\t/' >> counts.txt
        zcat ~{small_variants_vcf} | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/snps\tall\t/' | awk 'BEGIN {OFS = "\t" ;}{print $1,$2,$3,"SNP",$4}' >> counts.txt
        $BEDTOOLS intersect -u -a ~{small_variants_vcf} -b $STR_BED | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/snps\tstr\t/' | awk 'BEGIN {OFS = "\t" ;}{print $1,$2,$3,"SNP",$4}' >> counts.txt
        $BEDTOOLS intersect -u -a <(cat <(zcat ~{small_variants_vcf} | $GREP "^#") <($BEDTOOLS intersect -v -a ~{small_variants_vcf} -b $STR_BED)) -b $SEGDUP_BED | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/snps\tsegDup\t/' | awk 'BEGIN {OFS = "\t" ;}{print $1,$2,$3,"SNP",$4}' >> counts.txt
        $BEDTOOLS intersect -v -a <(cat <(zcat ~{small_variants_vcf} | $GREP "^#") <($BEDTOOLS intersect -v -a ~{small_variants_vcf} -b $STR_BED)) -b $SEGDUP_BED | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | cut -f 10 | sort | uniq -c | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t//' | sed 's/^/snps\tnonRep\t/' | awk 'BEGIN {OFS = "\t" ;}{print $1,$2,$3,"SNP",$4}' >> counts.txt
        cat counts.txt | sed 's/^/~{assembly_name}\t/' > tmp
        mv tmp counts.txt
    >>>
    runtime {
        memory: "16G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File counts = "counts.txt"
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
        BGZIP=/opt/hall-lab/htslib-1.9/bin/bgzip
        mkdir tmp1 tmp2 tmp3 tmp4
        $BCFTOOLS isec -c none -p tmp1 ~{small_variants1_ref} ~{small_variants2_ref}
        cat tmp1/0003.vcf | sed '/^chr/ s|...$|1/1|g' | $BGZIP -c > ref1_ref2_homozygous.vcf.gz
        cat tmp1/0000.vcf | sed '/^chr/ s/...$/1|0/g' | $BGZIP -c > ref1_het.vcf.gz
        cat tmp1/0001.vcf | sed '/^chr/ s/...$/0|1/g' | $BGZIP -c > ref2_het.vcf.gz
        $TABIX -fp vcf ref1_ref2_homozygous.vcf.gz
        $TABIX -fp vcf ref1_het.vcf.gz
        $TABIX -fp vcf ref2_het.vcf.gz
        $BCFTOOLS concat -a -d all ref1_ref2_homozygous.vcf.gz ref1_het.vcf.gz ref2_het.vcf.gz -o small_variants.combined.vcf.gz -O z
        $TABIX -fp vcf small_variants.combined.vcf.gz

        $BCFTOOLS isec -c snps -p tmp2 ~{small_variants_self} ~{small_variants2_ref_by_query}
        $BCFTOOLS query -f '%ID\n' tmp2/0002.vcf > self_in_ref2
        $BCFTOOLS query -f '%ID\n' tmp2/0003.vcf > ref2_in_self
        $BCFTOOLS isec -c snps -p tmp3 ~{small_variants_self_by_query} ~{small_variants1_ref_by_query}
        $BCFTOOLS query -f '%ID\n' tmp3/0002.vcf > self_in_ref1
        $BCFTOOLS query -f '%ID\n' tmp3/0003.vcf > ref1_in_self
        cat self_in_ref2 self_in_ref1 > known_from_self
        cat ref2_in_self ref1_in_self > known_from_ref
        $BCFTOOLS view -i '%ID!=@known_from_self' ~{small_variants_self} | $BGZIP -c > self_novel_small_variants.vcf.gz
        $BCFTOOLS view -i '%ID==@known_from_self' ~{small_variants_self} | $BGZIP -c > self_known_small_variants.vcf.gz
        $BCFTOOLS view -i '%ID!=@known_from_ref' small_variants.combined.vcf.gz | $BGZIP -c > ref_unique_small_variants.vcf.gz
        $BCFTOOLS view -i '%ID==@known_from_ref' small_variants.combined.vcf.gz | $BGZIP -c > ref_nonunique_small_variants.vcf.gz
        $TABIX -fp vcf self_novel_small_variants.vcf.gz
        $TABIX -fp vcf self_known_small_variants.vcf.gz
        $TABIX -fp vcf ref_unique_small_variants.vcf.gz
        $TABIX -fp vcf ref_nonunique_small_variants.vcf.gz
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/bcftools_samtools@sha256:49dc9b0b9419281b87df4a9faf9ca6681725317108bca66ba37f9bd1d86e9de2"
    }
    output {
        File combined_vcf_ref = "small_variants.combined.vcf.gz"
        File combined_vcf_ref_index = "small_variants.combined.vcf.gz.tbi"
        File novel_vcf_self = "self_novel_small_variants.vcf.gz"
        File known_vcf_self = "self_known_small_variants.vcf.gz"
        File unique_vcf_ref = "ref_unique_small_variants.vcf.gz"
        File nonunique_vcf_ref = "ref_nonunique_small_variants.vcf.gz"
        File novel_vcf_self_index = "self_novel_small_variants.vcf.gz.tbi"
        File known_vcf_self_index = "self_known_small_variants.vcf.gz.tbi"
        File unique_vcf_ref_index = "ref_unique_small_variants.vcf.gz.tbi"
        File nonunique_vcf_ref_index = "ref_nonunique_small_variants.vcf.gz.tbi"
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
        FIND_DUPS=/opt/hall-lab/scripts/find_duplicate_markers.py
        #combine fasta files and sort by sequence
        cat ~{small_variants_ref} ~{small_variants_self} | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep -v "^$" | paste - - - - | awk -v OFS="\t" -v FS="\t" '{if($2<$4) {print($2, $4, $1, $3)} else{print($4,$2,$3,$1)}}' | sort | awk -v OFS="\n" -v FS="\t" '{print($3,$1,$4,$2)}' > tmp
        #find duplicate markers
        $PYTHON $FIND_DUPS -i tmp > small_variants.combined.fasta
        cat ~{small_variants_ref_marker_positions} ~{small_variants_self_marker_positions} | sort -u > small_variants.marker_positions.txt
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
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
        docker: "apregier/analyze_assemblies@sha256:54669591da03e517f61097f93f8eac512368ae503954276b0149b13ebae0aec4"
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
        docker: "apregier/analyze_assemblies@sha256:54669591da03e517f61097f93f8eac512368ae503954276b0149b13ebae0aec4"
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
        docker: "apregier/analyze_assemblies@sha256:54669591da03e517f61097f93f8eac512368ae503954276b0149b13ebae0aec4"
    }
    output {
        File bam = "aligned.bam"
    }
}
