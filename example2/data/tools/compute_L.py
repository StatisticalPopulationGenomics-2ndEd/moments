## compute the total length of sequence in several .bed files
## the parameter L is required to obtain the expected SFS

import numpy as np
import sys

from lib import *


def main():

    bed_files = sys.argv[1:]
    L_tot = 0
    for file in bed_files: 
        regions = read_bedfile(file)[0]
        raw_L = np.diff(regions, axis=1).sum()
        _regions = boolmask_to_regions(regions_to_boolmask(regions))
        L = int(np.diff(_regions, axis=1).sum())
        if raw_L != L:
            raise ValueError(f"file {file} has overlapping regions")
        print(f"L_{file} =\t{L}")
        L_tot += L
    print(f"L_tot =\t{L_tot}")

    return


if __name__ == "__main__":
    main()
