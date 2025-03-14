## add together spectra from several saved .sfs files and write the sum

import argparse
from collections import defaultdict
import moments 
import numpy as np
import pickle 


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--sfs_files", type=str, required=True, nargs="*"
    )
    parser.add_argument(
        "-o", "--out_file", type=str, required=True
    )

    return parser.parse_args()


def main():

    args = get_args()
    sfs_files = args.sfs_files
    out_file = args.out_file
    to_pops = ["MSL", "CHS", "GBR", "Vindija"]
    sum_arr = 0
    for file in sfs_files:
        sfs = moments.Spectrum.from_file(file)
        pop_ids = sfs.pop_ids
        to_marg = [pop_ids.index(pop) for pop in pop_ids if pop not in to_pops]
        _sfs = sfs.marginalize(to_marg)
        sum_arr += np.asarray(_sfs)
    pop_ids = _sfs.pop_ids
    sum_sfs = moments.Spectrum(sum_arr, mask_corners=True, pop_ids=pop_ids)
    sum_sfs.to_file(out_file)

    return


if __name__ == "__main__":
    main()

