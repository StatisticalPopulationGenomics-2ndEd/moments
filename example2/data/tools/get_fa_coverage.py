## write a .bed file recording the regions for which high-confidence assignments
## of ancestral state (denoted by capital letters) exist in an input .fa file.

import argparse
import gzip
import numpy as np


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fasta_file", required=True
    )    
    parser.add_argument(
        "-o", "--out_file", required=True
    )
    return parser.parse_args()


def load_fa_array(file):
    """
    Load the characters in a .fa file and return them as a numpy array.
    """
    lines = []
    line0 = None
    with open(file, "r") as fin:
        for line in fin:
            if ">" in line:
                line0 = line
                continue
            lines.append(line.rstrip("\n"))
    array = np.array([c for c in "".join(lines)])

    return line0, array


def boolmask_to_regions(boolmask):
    """
    Transform a boolean mask where 1 is 'masked out' to an array of mask regions
    """
    jumps = np.diff(np.concatenate(([1], boolmask, [1])))
    starts = np.where(jumps == -1)[0]
    ends = np.where(jumps == 1)[0]
    regions = np.stack([starts, ends], axis=1)

    return regions


def write_bedfile(file, regions, chromnum):
    """
    write a .bed file (supports .gz) from an array of mask regions.
    """
    openfunc = gzip.open if file.endswith(".gz") else open 
    with openfunc(file, "wb") as fout:
        for start, end in regions:
            lineb = f"{chromnum}\t{start}\t{end}\n".encode()
            fout.write(lineb)
    return


def main():

    args = get_args()
    fasta_file = args.fasta_file
    out_file = args.out_file
    line0, array = load_fa_array(fasta_file)
    # retrieving the chromosome number in this way will only work when the label
    # format is like this: >ANCESTOR_for_chromosome:GRCh38:20:1:64444167:1\n
    chromnum = f"chr{line0.split(':')[2]}"
    boolmask = np.ones(len(array), dtype=np.uint8)
    for nt in ("A", "T", "C", "G"):
        boolmask[array == nt] = 0
    regions = boolmask_to_regions(boolmask)
    write_bedfile(out_file, regions, chromnum)

    return


if __name__ == "__main__":
    main()
