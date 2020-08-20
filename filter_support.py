import argparse

def run(input_file):
    with open(input_file, "r") as fp:
        for line in fp:
            line = line.strip()
            fields = line.split()
            if int(fields[2]) + int(fields[5]) == 1:
                print("\t".join(fields))

def main():

    usage = """%prog -i <input file>
filter_support.py

Author: Allison Regier
Description: Input file gives distance from an aligned read to alt and ref versions of variants.  Output only variants with unambiguous distance (0 1 or 1 0)
    """
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("-i", dest="inputFile", help="input file", metavar="FILE")
    args = parser.parse_args()

    if args.inputFile is None:
        parser.print_help()
    else:
        run(args.inputFile)

if __name__ == "__main__":
    main()
