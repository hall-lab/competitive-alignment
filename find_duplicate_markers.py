#!/usr/bin/env python
from __future__ import print_function
import sys
from optparse import OptionParser

duplicates = {}
def add_duplicates(current_marker, previous_marker):
    sequence1 = "{}.{}".format(previous_marker[1], previous_marker[3])
    sequence2 = "{}.{}".format(previous_marker[3], previous_marker[1])
    if sequence1 in duplicates:
        duplicates[sequence1].add(current_marker[0])
        duplicates[sequence1].add(current_marker[2])
    elif sequence2 in duplicates:
        duplicates[sequence2].add(current_marker[0])
        duplicates[sequence2].add(current_marker[2])
    else:
        duplicates[sequence] = set([current_marker[0], current_marker[2]])

def process_fasta(fasta_file):
    previous_marker = []
    current_marker = []

    with open(fasta_file, "r") as fa:
        lines = []
        for line in fa:
            lines.append(str(line.strip()))
            if len(lines) >= 4:
                if current_marker == []:
                    current_marker = lines
                else:
                    previous_marker = current_marker
                    current_marker = lines
                    if current_marker[1] == previous_marker[1] and current_marker[3] == previous_marker[3]:
                        add_duplicates(current_marker, previous_marker)
                lines = []
        if len(lines) > 0:
            raise Exception("Leftover lines in fasta, make sure ref and alt fastas are paired for each marker")

    fasta_to_remove = set()
    for duplicate_set in duplicates.values():
        for duplicate in duplicate_set:
            fasta_to_remove.add(duplicate)
    print("{} duplicate sequences removed".format([str(len(fasta_to_remove))]), file=sys.stderr)
    with open(fasta_file, "r") as fa:
        lines = []
        for line in fa:
            lines.append(str(line.strip()))
            if len(lines) >= 2:
                if not (lines[0] in fasta_to_remove):
                    print(lines[0])
                    print(lines[1])
                lines = []
        if len(lines) > 0:
            raise Exception("Leftover lines in fasta, make sure fastas are complete")

def main():
    usage = """%prog -i <fasta>

Given a fasta file of paired marker sequences, sorted by sequence, remove duplicated sequences.
Author: Allison Regier
    """

    parser = OptionParser(usage)
    parser.add_option("-i", dest="fasta",
        help="a fasta file",
        metavar="FILE")
    (opts, args) = parser.parse_args()
    if opts.fasta is None:
        parser.print_help()
        print()
    else:
        process_fasta(opts.fasta)

if __name__ == "__main__":
    sys.exit(main())
