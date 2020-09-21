version 1.0
import "call_assembly_variants.wdl" as call_variants
import "genotype_markers.wdl" as genotype

workflow CompetitiveGenotyping {
    input {
        File assembly_list
        File dataset_list
        File ref
        File ref_index
        String ref_name
        File segdup_bed
        File str_bed
    }
    Array[Array[File]] assemblies = read_tsv(assembly_list)
    Array[Array[File]] datasets = read_tsv(dataset_list)

    scatter (assembly in assemblies) {
        call call_variants.CallAssemblyVariants as call_variants {
            input:
            assembly_name=assembly[0],
            contigs1=assembly[1],
            contigs2=assembly[2],
            ref=ref,
            ref_index=ref_index,
            ref_name=ref_name,
            segdup_bed=segdup_bed,
            str_bed=str_bed
#            fastq_list=assembly[3]
        }
    }

#    call merge_variants {
#        input:
#            small_variants=call_variants.small_variants,
#            small_variants_marker_positions=call_variants.small_variants_marker_positions #,
#            #sv=call_variants.sv  #TODO
#    }

#    scatter (dataset in datasets) {
#        call genotype.GenotypeMarkers as genotype {
#            input:
#            dataset_name=dataset[0],
#            dataset_fastq=dataset[1],
#            variant_fasta=merge_variants.fasta_representation,
#            marker_positions=merge_variants.marker_positions
#        }
#    }

    output {
        Array[File] sv_ref1 = call_variants.sv_ref1
        Array[File] sv_ref2 = call_variants.sv_ref2
        Array[File] sv_self = call_variants.sv_self
        Array[File] small_variants_ref1 = call_variants.small_variants_ref1
        Array[File] small_variants_ref2 = call_variants.small_variants_ref2
        Array[File] small_variants_self = call_variants.small_variants_self
        Array[File] small_variants_ref_combined = call_variants.small_variants_ref_combined
        Array[File] small_variants_ref_combined_index = call_variants.small_variants_ref_combined_index
        Array[File] sv_combined = call_variants.sv_combined
#        Array[File] marker_counts = genotype.marker_counts
    }
}

task merge_variants {
    input {
        Array[File] small_variants
        Array[File] small_variants_marker_positions
        #Array[File] sv #TODO
    }
    command <<<
        PYTHON=/opt/hall-lab/python-2.7.15/bin/python
        FIND_DUPS=/opt/hall-lab/scripts/find_duplicate_markers.py
        cat ~{sep=" " small_variants} | paste - - - - | awk -v OFS="\t" -v FS="\t" '{if($2<$4) {print($2, $4, $1, $3)} else{print($4,$2,$3,$1)}}' | sort | awk -v OFS="\n" -v FS="\t" '{print($3,$1,$4,$2)}' > tmp
        $PYTHON $FIND_DUPS -i tmp > variants_merged.fasta
        cat ~{sep=" " small_variants_marker_positions} | sort -u > marker_positions.txt
    >>>
    runtime {
        docker: "apregier/analyze_assemblies@sha256:4cd67e009ae65820772265b572fc8cb9ce9e6e09228d1d73ee1f5d9118e91fca"
        memory: "64 GB"
    }
    output {
        File fasta_representation="variants_merged.fasta"
        File marker_positions="marker_positions.txt"
    }
}
