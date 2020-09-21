version 1.0

workflow GenotypeMarkers {
    input {
        String dataset_name
        File dataset_fastq
        File variant_fasta
        File marker_positions
    }

    call align_dataset {
        input:
            dataset_name=dataset_name,
            dataset_fastq=dataset_fastq,
            variant_fasta=variant_fasta
    }

    call genotype {
        input:
            alignment=align_dataset.alignment,
            marker_positions=marker_positions
    }
    output {
        File marker_counts = genotype.marker_counts
    }
}

task align_dataset {
    input {
        File dataset_fastq
        File variant_fasta
        String dataset_name
    }
    command <<<
        set -exo pipefail
        MINIMAP2=/opt/hall-lab/minimap2/minimap2
        $MINIMAP2 -x sr -N 5 --cs ~{variant_fasta} ~{dataset_fastq} > alignment.paf
    >>>
    runtime {
        memory: "32G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File alignment="alignment.paf"
    }
}

task genotype {
    input {
        File alignment
        File marker_positions
    }
    command <<<
        PYTHON=/opt/hall-lab/python-2.7.15/bin/python
        BEDTOOLS=/opt/hall-lab/bedtools
        MARKER_COUNTS_SCRIPT=/opt/hall-lab/scripts/marker_counts_from_paf.py
        $PYTHON $MARKER_COUNTS_SCRIPT -a ~{alignment} -m ~{marker_positions} > read_marker_pairs.txt
        grep non-ambiguous read_marker_pairs.txt | cut -f 6-8 | sort | $BEDTOOLS groupby -g 1 -c 2,3 -o sum,sum > marker_counts.txt
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
    }
    output {
        File marker_counts="marker_counts.txt"
    }
}
