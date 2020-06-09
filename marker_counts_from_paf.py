import paf
import sys
import argparse

class Marker:
    def __init__(self, distance, contig_name, position_on_read):
        self.distance = distance
        self.contig_name = contig_name
        self.position_on_read = position_on_read

class MarkerPair:
    def __init__(self, name):
        self.a = None
        self.r = None
        self.name = name

class Alignment:
    def __init__(self, allele_id, distance, contig_name, position_on_read):
        self.allele_id = allele_id
        self.distance = distance
        self.contig_name = contig_name
        self.position_on_read = position_on_read

class ReadAlignments:
    def __init__(self, name):
        self.markers = {}
        self.name = name

    def by_marker(self):
        markerPairs = set()
        for marker in self.markers:
            markerPair = MarkerPair(marker)
            for alignment in self.markers[marker]:
                markerObj = Marker(alignment.distance, alignment.contig_name, alignment.position_on_read)
                if alignment.allele_id == 'a':
                    markerPair.a = markerObj
                if alignment.allele_id == 'r':
                    markerPair.r = markerObj
            markerPairs.add(markerPair)
        return markerPairs

    def add_alignment(self, pafLine, markers_on_contigs):
        if pafLine.rname in markers_on_contigs:
            for marker_position in markers_on_contigs[pafLine.rname]:
                pos = int(marker_position[0])-1
                if pafLine.isInAlignedRangeRef(pos):
                    distance = pafLine.edit_distance_for_region_in_paf_line(
                                     pos,
                                     pos+1)
                    base_marker_id = marker_position[3]
                    query_pos_set = pafLine.query_pos_for_ref_pos_in_paf(pos)
                    if len(query_pos_set) == 1:
                        query_pos = query_pos_set.pop()
                        newAlignment = Alignment(marker_position[4], distance, pafLine.rname, query_pos)
                        if base_marker_id in self.markers:
                            self.markers[base_marker_id].add(newAlignment)
                        else:
                            self.markers[base_marker_id] = set([newAlignment])

    def print_alignments(self):
        markers_to_connect = set()
        markers = self.by_marker()
        for marker in markers:
            if marker.a != None and marker.r != None and marker.a.distance + marker.r.distance == 1:
                print("\t".join(map(str, ["Marker", "non-ambiguous", self.name, marker.a.contig_name, marker.r.contig_name, marker.name, marker.a.distance, marker.r.distance])))
            elif marker.a == None:
                print("\t".join(map(str, ["Marker", "no a", self.name, "NA", marker.r.contig_name, marker.name, "NA", marker.r.distance])))
            elif marker.r == None:
                print("\t".join(map(str, ["Marker", "no r", self.name, marker.a.contig_name, "NA", marker.name, marker.a.distance, "NA"])))
            else:
                print("\t".join(map(str, ["Marker", "distance!=1", self.name, marker.a.contig_name, marker.r.contig_name, marker.name, marker.a.distance, marker.r.distance])))

def run(alignment_file, marker_file):
    markers_on_contigs = {}
    with open(marker_file, 'r') as mp:
        for line in mp:
            line = line.strip()
            fields = line.split()
            if fields[0] in markers_on_contigs:
                markers_on_contigs[fields[0]].append([fields[1], fields[2], fields[3], fields[4]])
            else:
                markers_on_contigs[fields[0]] = [[fields[1], fields[2], fields[3], fields[4]]]

    read_alignments = None
    with open(alignment_file, 'r') as fp:
        counter = 0
        for line in fp:
            line = line.strip()
            pafLine = paf.Paf(line)
            if read_alignments == None:
                read_alignments = ReadAlignments(pafLine.qname)
            elif read_alignments.name != pafLine.qname:
                read_alignments.print_alignments()
                read_alignments = ReadAlignments(pafLine.qname)
            if pafLine.qlength > 20000:
                read_alignments.add_alignment(pafLine, markers_on_contigs)
    if read_alignments != None:
        read_alignments.print_alignments()

def main():

    usage = """%prog -a <alignment file> -m <positions on contigs file>
marker_counts_from_paf.py

Author: Allison Regier
Description: Given an alignment of long reads to contigs, and positions of SNP markers on those contigs, classify the read as supporting the ref or alt allele of each SNP marker.  Some filtering is performed - for a read to classify a marker, it needs to align to both ref and alt alleles, and it has to match one allele exactly at the SNP position, and one allele must be a mismatch.
For each pair of SNP markers covered by a single read, classify what phasing the read supports (a-a, a-r, r-a, or r-r) and emit the distance between the markers on the read.
    """
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("-a", dest="alignmentFile", help="Alignment file (paf)", metavar="FILE")
    parser.add_argument("-m", dest="markerFile", help="Positions of markers on contigs.  Fields are contig, position, marker", metavar="FILE")
    args = parser.parse_args()

    if args.alignmentFile is None or args.markerFile is None:
        parser.print_help()
    else:
        run(args.alignmentFile, args.markerFile)

if __name__ == "__main__":
    main()
