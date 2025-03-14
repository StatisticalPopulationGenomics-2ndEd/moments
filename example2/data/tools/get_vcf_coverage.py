## write a .bed file recording regions covered in a .vcf file.
## intended to record coverage of GVCF files. 
## usage: python get_vcf_coverage.py -i {input.vcf.gz} -o {output.bed.gz}
## written 25-02-2025 

import argparse 
import gzip
import numpy as np 


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--in_file", required=True
    )    
    parser.add_argument(
        "-o", "--out_file", required=True
    )
    return parser.parse_args()


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
    in_file = args.in_file
    out_file = args.out_file
    i = 0
    start = None
    lastpos = None
    regions = []
    openfunc = gzip.open if in_file.endswith(".gz") else open
    with openfunc(in_file, "rb") as fin:
        for bline in fin:
            line = bline.decode()
            if line.startswith("#"):
                continue
            pos = int(line.split()[1])
            if i == 0:
                start = pos
            else:
                if pos - lastpos > 1:
                    regions.append([start, lastpos])
                    start = pos
            i += 1
            lastpos = pos
    regions.append([start, lastpos])
    region_arr = np.array(regions, dtype=np.int64)
    # recall that .bed file `start` positions are 0-indexed
    region_arr[:, 0] -= 1
    chromnum = line.split()[0]
    write_bedfile(out_file, region_arr, chromnum)

    return 


if __name__ == "__main__":
    main()

