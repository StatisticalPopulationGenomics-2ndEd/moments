## Write a file recording the estimated ancestral nucleotide state from a .fa
## file for every site in an input .vcf file. 
## output file has columns `chrom` `pos` `state`

## This script will throw an error if any .vcf sites are not assigned with 
## high confidence (denoted by capitalized nucleotide codes e.g. A) in .fa file.
## It is assumed that the input .vcf has already been masked to account for this

import argparse
import gzip 
import numpy as np

from lib import *


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--vcf_file", required=True
    )
    parser.add_argument(
        "-f", "--fasta_file", required=True
    )
    parser.add_argument(
        "-o", "--out_file", required=True, help=".tab or .tab.gz format"
    )

    return parser.parse_args()


def main():

    args = get_args()
    vcf_file = args.vcf_file
    fasta_file = args.fasta_file
    out_file = args.out_file
    fa_array = load_fa_array(fasta_file)[1]
    positions = []
    states = []
    openfunc = gzip.open if vcf_file.endswith(".gz") else open 
    with openfunc(vcf_file, "rb") as fin:
        for lineb in fin:
            line = lineb.decode()
            if "#" in line:
                continue
            chromnum, pos = line.split()[:2]
            # .fa file treated as 0-indexed
            position0 = int(pos) - 1
            state = fa_array[position0]
            if state not in ("A", "T", "G", "C"):
                raise ValueError(
                    f"pos {position0} lacks high-confidence state assignment"
                )
            # use 1-indexed position for output
            positions.append(pos)
            states.append(state)
    write_tab_file(out_file, positions, states, chromnum)

    return


if __name__ == "__main__":
    main()
