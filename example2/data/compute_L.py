# compute L for each window and save the result

import numpy as np
import moments 
import pandas


bed_file = lambda chrom: f"bed_files/combined_mask_chr{chrom}.bed.gz"


# interval size
l = 30000000


def load_genome_file(fname):
    ret = dict()
    with open(fname, "r") as fin:
        for line in fin:
            chrom, L = line.split()
            ret[chrom] = int(L)
    return ret 


genome_file = "hg19.genome"
seq_lens = load_genome_file(genome_file)


# we loop over chromosomes and intervals
Ls = {"region": [], "interval": [], "L": []}
for chrom in range(1, 23):
    length = seq_lens[f"chr{chrom}"]
    intervals = [[x, y] for x, y 
                 in zip(range(1, length, l), range(l + 1, length + l, l))]
    for ii, interval in enumerate(intervals):
        L = moments.Parsing.compute_L(bed_file(chrom), interval=interval)
        Ls["region"].append(f"chrom_{chrom}.region_{ii}")
        Ls["interval"].append(interval)
        Ls["L"].append(L)


df = pandas.DataFrame(Ls)
df.to_csv("region_L_tbl.csv", index=False)

