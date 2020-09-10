version 1.0

workflow ConvertToFasta {
    input {
        File vcf
        File vcf_index
        File ref
        File ref_index
        File ref_name
    }

    call get_regions {
        input:
            vcf=vcf,
            vcf_index=vcf_index
    }

    scatter (region in get_regions.regions) {
        call convert {
            input:
                region_file=region,
                vcf=get_regions.filtered_vcf,
                vcf_index=get_regions.filtered_vcf_index,
                ref=ref,
                ref_index=ref_index,
                ref_name=ref_name
        }
    }

    call combine {
        input:
            split_fastas = convert.fasta,
            split_marker_positions = convert.marker_positions
    }

    output {
        File marker_fasta = combine.fasta
        File marker_positions = combine.marker_positions
    }
}

task get_regions {
    input {
        File vcf
        File vcf_index
    }
    command <<<
        set -exo pipefail
        BCFTOOLS=/opt/hall-lab/bcftools-1.9/bin/bcftools
        TABIX=/opt/hall-lab/htslib-1.9/bin/tabix
        PERL=/usr/bin/perl
        PRINT_REGIONS=/opt/hall-lab/scripts/printRegions2.pl
        $BCFTOOLS view -m2 -M2 -v snps ~{vcf} -o tmp.vcf.gz -O z
        $TABIX -fp vcf tmp.vcf.gz
        zcat tmp.vcf.gz | $BCFTOOLS query -f '%CHROM\t%POS\n' | $PERL $PRINT_REGIONS > regions.txt
        split -l 5000 -a 5 -d regions.txt regions.sub
        ls regions.sub* > split_regions_files.txt
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/bcftools_samtools@sha256:49dc9b0b9419281b87df4a9faf9ca6681725317108bca66ba37f9bd1d86e9de2"
    }
    output {
        Array[File] regions = read_lines("split_regions_files.txt")
        File filtered_vcf = "tmp.vcf.gz"
        File filtered_vcf_index = "tmp.vcf.gz.tbi"
    }
}

task convert {
    input {
        File region_file
        File vcf
        File vcf_index
        File ref
        File ref_index
        String ref_name
    }
    command <<<
        set -exo pipefail

        BCFTOOLS=/opt/hall-lab/bcftools-1.9/bin/bcftools
        SAMTOOLS=/opt/hall-lab/samtools-1.9/bin/samtools
        OUT="split.fasta"
        MARKERS="split.marker_positions"
        ln -s ~{ref} ref.fa
        ln -s ~{ref_index} ref.fa.fai
        regions=`cat ~{region_file}`
        for region in $regions; do
            $SAMTOOLS faidx ref.fa $region | $BCFTOOLS consensus ~{vcf} | sed 's/^>/>alt_~{ref_name}_/' >> $OUT
            $SAMTOOLS faidx ref.fa $region | sed 's/^>/>ref_~{ref_name}_/' >> $OUT
            echo "ref_~{ref_name}_$region	100	ref_~{ref_name}_$region	~{ref_name}_$region	r" >> $MARKERS
            echo "alt_~{ref_name}_$region	100	alt_~{ref_name}_$region	~{ref_name}_$region	a" >> $MARKERS ##TODO double check 100 or 101?
        done
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/bcftools_samtools@sha256:49dc9b0b9419281b87df4a9faf9ca6681725317108bca66ba37f9bd1d86e9de2"
    }
    output {
        File fasta = "split.fasta"
        File marker_positions = "split.marker_positions"
    }
}

task combine {
    input {
        Array[File] split_fastas
        Array[File] split_marker_positions
    }
    command <<<
        set -exo pipefail
        cat ~{sep=" " split_fastas} > combined.fasta
        cat ~{sep=" " split_marker_positions} > combined_marker_positions.txt
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/analyze_assemblies@sha256:cae6b31b36f8f49fcd1fcba8ba18107d4e0d7ad967600514d292423300c52425"
    }
    output {
        File fasta = "combined.fasta"
        File marker_positions = "combined_marker_positions.txt"
    }
}
