import paf
import vcf
import sys
import argparse

def run(alignment_file, variant_file):
    markers_on_contigs = {}
    with open(variant_file, 'r') as mp:
        vcf_reader = vcf.Reader(mp)
        for record in vcf_reader:
            if record.CHROM in markers_on_contigs:
                markers_on_contigs[record.CHROM].append([record.POS, record.ID])
            else:
                markers_on_contigs[record.CHROM] = [[record.POS, record.ID]]

    with open(alignment_file, 'r') as fp:
        counter = 0
        for line in fp:
            line = line.strip()
            pafLine = paf.Paf(line)
            if pafLine.rname in markers_on_contigs:
                for marker_position in markers_on_contigs[pafLine.rname]:
                    pos = int(marker_position[0])-1
                    if pafLine.isInAlignedRangeRef(pos):
                        distance = pafLine.edit_distance_for_region_in_paf_line(
                                         pos,
                                         pos+1)
                        base_marker_id = marker_position[1]
                        query_pos_set = pafLine.query_pos_for_ref_pos_in_paf(pos)
                        if len(query_pos_set) == 1:
                            query_pos = query_pos_set.pop()
                            print("\t".join(map(str, [base_marker_id, pafLine.qname, distance, query_pos])))

def main():

    usage = """%prog -a <alignment file> -v <variants file>
read_support_from_paf.py

Author: Allison Regier
Description: Given an alignment of long reads to contigs, and variants on the contigs, print the distance between read and contig at the variant position.
    """
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("-a", dest="alignmentFile", help="Alignment file (paf) reads to contigs", metavar="FILE")
    parser.add_argument("-v", dest="variantFile", help="Variant vcf.  Position of variant on contig specified by QNAME and QSTART in info field", metavar="FILE")
    args = parser.parse_args()

    if args.alignmentFile is None or args.variantFile is None:
        parser.print_help()
    else:
        run(args.alignmentFile, args.variantFile)

if __name__ == "__main__":
    main()
