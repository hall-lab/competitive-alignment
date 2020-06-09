version 1.0

workflow ConvertToFasta {
    input {
        File vcf
        File vcf_index
        File query
        File ref
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
                query=query,
                ref=ref,
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
        PRINT_REGIONS=/storage1/fs1/ccdg/Active/analysis/ref_grant/assembly_analysis_20200220/multiple_competitive_alignment/printRegions2.pl #TODO
        $BCFTOOLS view -f PASS -g het -v snps ~{vcf} -o tmp.vcf.gz -O z
        $TABIX -fp vcf tmp.vcf.gz
        zcat tmp.vcf.gz | $BCFTOOLS query -f '%CHROM\t%POS\n' | $PERL $PRINT_REGIONS > regions.txt
        split -l 5000 -a 3 -d regions.txt regions.sub
        ls regions.sub* > split_regions_files.txt
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
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
        File query
        File ref
        String ref_name
    }
    command <<<
        set -exo pipefail

        BCFTOOLS=/opt/hall-lab/bcftools-1.9/bin/bcftools
        SAMTOOLS=/opt/hall-lab/samtools-1.9/bin/samtools
        OUT="split.fasta"
        MARKERS="split.marker_positions"
        regions=`cat ~{region_file}`
        for region in $regions; do
            $SAMTOOLS faidx ~{ref} $region | $BCFTOOLS consensus ~{vcf} | sed 's/^>/>alt_~{ref_name}_/' >> $OUT
            $SAMTOOLS faidx ~{ref} $region | sed 's/^>/>ref_~{ref_name}_/' >> $OUT
            echo "ref_~{ref_name}_$region	100	ref_~{ref_name}_$region	~{ref_name}_$region	r" >> $MARKERS
            echo "alt_~{ref_name}_$region	100	alt_~{ref_name}_$region	~{ref_name}_$region	a" >> $MARKERS ##TODO double check 100 or 101?
        done
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/bcftools_samtools:latest"
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
        cat ~{split_fastas} > combined.fasta
        cat ~{split_marker_positions} > combined_marker_positions.txt
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File fasta = "combined.fasta"
        File marker_positions = "combined_marker_positions.txt"
    }
}
