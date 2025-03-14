## perform mask arithmetic, with the following steps:
## take the intersection of masks flagged `--isec`
## take the union of masks flagged `--remove` [exonic/genic regions etc]
## extend `--remove` mask by `flank` bp
## remove sites in extended mask from intersected mask
## output remaining mask

## Because we use uint8 types to count overlaps, this script cannot handle more
## than 256 masks at once : )

import argparse 
import gzip
import numpy as np 

from lib import * 


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--isec", nargs="*", required=True
    )    
    parser.add_argument(
        "-s", "--subtract", nargs="*", default=[]
    )     
    parser.add_argument(
        "-flank", "--flank", type=float, default=None
    )  
    parser.add_argument(
        "-o", "--out_file", required=True
    )
    return parser.parse_args()


def main():

    args = get_args()
    isec_files = args.isec
    sub_files = args.subtract
    flank = int(args.flank) if args.flank is not None else None
    out_file = args.out_file
    chromnum = read_bedfile(isec_files[0])[1]
    isec_regions = [read_bedfile(file)[0] for file in isec_files]
    isec = intersect_regions(isec_regions)
    if len(sub_files) > 0:
        sub_regions = [read_bedfile(file)[0] for file in sub_files]
        union = union_regions(sub_regions)
        if flank is not None:   
            union = flank_regions(union, flank)
        regions = subtract_regions(isec, union)
    else:
        regions = isec
    write_bedfile(out_file, regions, chromnum)

    return


if __name__ == "__main__":
    main()
