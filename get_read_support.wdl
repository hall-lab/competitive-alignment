version 1.0

workflow GetReadSupport {
    input {
        File ref
        File contigs1
        File contigs2
        File reads_list
        File small_variants_ref_contigs1
        File small_variants_contigs1_ref
        File small_variants_ref_contigs2
        File small_variants_contigs2_ref
        File small_variants_contigs1_contigs2
        File small_variants_contigs2_contigs1
    }
    Array[File] reads = read_lines(reads_list)
    scatter (reads_file in reads) {
        call align_source_reads as align_ref {
            input:
                ref=ref,
                reads=reads_file,
        }

        call align_source_reads as align_contigs1 {
            input:
                ref=contigs1,
                reads=reads_file
        }

        call align_source_reads as align_contigs2 {
            input:
                ref=contigs2,
                reads=reads_file
        }

        call get_read_support_small as get_read_support_small_ref_contigs1 {
            input:
                alignment=align_ref.alignment,
                variants=small_variants_ref_contigs1
        }

        call get_read_support_small as get_read_support_small_contigs1_ref {
            input:
                alignment=align_contigs1.alignment,
                variants=small_variants_contigs1_ref
        }

        call collate_read_support as collate_small_variant_support_ref_contigs1 {
            input:
                contig_support=get_read_support_small_ref_contigs1.read_support,
                ref_support=get_read_support_small_contigs1_ref.read_support
        }

        call get_read_support_small as get_read_support_small_ref_contigs2 {
            input:
                alignment=align_ref.alignment,
                variants=small_variants_ref_contigs2
        }

        call get_read_support_small as get_read_support_small_contigs2_ref {
            input:
                alignment=align_contigs2.alignment,
                variants=small_variants_contigs2_ref
        }

        call collate_read_support as collate_small_variant_support_ref_contigs2 {
            input:
                contig_support=get_read_support_small_ref_contigs2.read_support,
                ref_support=get_read_support_small_contigs2_ref.read_support
        }
        call get_read_support_small as get_read_support_small_contigs1_contigs2 {
            input:
                alignment=align_contigs2.alignment,
                variants=small_variants_contigs1_contigs2
        }

        call get_read_support_small as get_read_support_small_contigs2_contigs1 {
            input:
                alignment=align_contigs1.alignment,
                variants=small_variants_contigs2_contigs1
        }

        call collate_read_support as collate_small_variant_support_contigs1_contigs2 {
            input:
                contig_support=get_read_support_small_contigs1_contigs2.read_support,
                ref_support=get_read_support_small_contigs2_contigs1.read_support
        }
    }
    call combine_read_support as combine_small_variant_support_ref_contigs1 {
        input:
            read_support = collate_small_variant_support_ref_contigs1.collated_support
    }
    call combine_read_support as combine_small_variant_support_ref_contigs2 {
        input:
            read_support = collate_small_variant_support_ref_contigs2.collated_support
    }
    call combine_read_support as combine_small_variant_support_contigs1_contigs2 {
        input:
            read_support = collate_small_variant_support_contigs1_contigs2.collated_support
    }

    call correspond_variant_ids {
        input:
            ref1_vcf = small_variants_ref_contigs1,
            ref1_by_query_vcf = small_variants_contigs1_ref,
            ref2_vcf = small_variants_ref_contigs2,
            ref2_by_query_vcf = small_variants_contigs2_ref,
            self_vcf = small_variants_ref_contigs1_contigs2,
            self_by_query_vcf = small_variants_contigs2_contigs1
    }

    output {
        File small_variant_support_ref_contigs1 = combine_small_variant_support_ref_contigs1.total_read_support
        File small_variant_support_ref_contigs2 = combine_small_variant_support_ref_contigs2.total_read_support
        File small_variant_support_contigs1_contigs2 = combine_small_variant_support_contigs1_contigs2.total_read_support
    }
}

task correspond_variant_ids {
    input {
        File ref1_vcf
        File ref1_by_query_vcf
        File ref2_vcf
        File ref2_by_query_vcf
        File self_vcf
        File self_by_query_vcf
    }
    command <<<
        set -exo pipefail
        join -j2 <(zcat ~{ref1_vcf} | grep -v "^#" | sed 's/	/_/' | cut -f 1,2 | sort -k2) <(zcat ~{ref1_by_query_vcf} | grep -v "^#" | sed 's/	/_/' | cut -f 1,2 | sort -k2) > ref1.txt
        join -j2 <(zcat ~{ref2_vcf} | grep -v "^#" | sed 's/	/_/' | cut -f 1,2 | sort -k2) <(zcat ~{ref2_by_query_vcf} | grep -v "^#" | sed 's/	/_/' | cut -f 1,2 | sort -k2) > ref2.txt
        join -j2 <(zcat ~{self_vcf} | grep -v "^#" | sed 's/v   /_/' | cut -f 1,2 | sort -k2) <(zcat ~{self_by_query_vcf} | grep -v "^#" | sed 's/	/_/' | cut -f 1,2 | sort -k2) > self.txt
        join -1 3 -2 3 <(sort -k3 ref1.txt) <(sort -k3 self.txt) > ref1_self.txt
        join -1 3 -2 2 <(sort -k3 ref2.txt) <(sort -k2 self.txt) > ref2_self.txt
        join -1 2 -2 2 <(sort -k2 ref1.txt) <(sort -k2 ref2.txt) > ref1_ref2.txt
    >>>
    runtime {
        memory: "16G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File ref1_self = "ref1_self.txt"
        File ref2_self = "ref2_self.txt"
        File ref1_ref2 = "ref1_ref2.txt"
    }
}

task align_source_reads {
    input {
        File ref
        File reads
    }
    command <<<
        set -exo pipefail
        MINIMAP2=/opt/hall-lab/minimap2/minimap2
        $MINIMAP2 -t 8 -I 8G -x map-pb --cs ~{ref} ~{reads} > alignment.paf
    >>>
    runtime {
        memory: "32G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File alignment="alignment.paf"
    }
}

task get_read_support_small {
    input {
        File alignment
        File variants
    }
    command <<<
        set -exo pipefail
        PYTHON=/opt/hall-lab/python-2.7.15/bin/python
        READ_COUNT_SCRIPT=/opt/hall-lab/scripts/read_support_from_paf.py
        $PYTHON $READ_COUNT_SCRIPT -a ~{alignment} -v ~{variants} > support.txt
    >>>
    runtime {
        memory: "32G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File read_support="support.txt"
    }
}

task collate_read_support {
    input {
        File contig_support
        File ref_support
    }
    command <<<
        set -exo pipefail
        PYTHON=/opt/hall-lab/python-2.7.15/bin/python
        BEDTOOLS=/opt/hall-lab/bedtools
        FILTER_SUPPORT_SCRIPT=/opt/hall-lab/scripts/filter_support.py
        join <(cat ~{contig_support} | sed 's/\t/ /' | sort) <(cat ~{ref_support} | sed 's/\t/ /' | sort) | sed 's/ /\t/' > collated.txt
        $PYTHON $FILTER_SUPPORT_SCRIPT -i collated.txt | sort | $BEDTOOLS groupby -g 1 -c 3,6 -o sum,sum > support.txt
    >>>
    runtime {
        memory: "32G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File collated_support="support.txt"
    }
}

task combine_read_support {
    input {
        Array[File] read_support
    }
    command <<<
        set -exo pipefail
        BEDTOOLS=/opt/hall-lab/bedtools
        cat ~{sep=" " read_support} | sort | $BEDTOOLS groupby -g 1 -c 2,3 -o sum,sum > support.txt
    >>>
    runtime {
        memory: "32G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File total_read_support="support.txt"
    }
}
