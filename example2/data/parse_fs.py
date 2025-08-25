# compute the frequency spectra across all chromosomes and save them

import numpy as np
import moments 


pop_file = "populations.txt"
vcf_file = lambda chrom: f"vcf_files/annotated_variants_chr{chrom}.vcf.gz"
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
for chrom in range(1, 23):
    L = seq_lens[f"chr{chrom}"]
    intervals = [[x, y] for x, y in zip(range(1, L, l), range(l + 1, L + l, l))]
    for ii, interval in enumerate(intervals):
        try:
            fs = moments.Spectrum.from_vcf(
                vcf_file(chrom),
                bed_file=bed_file(chrom),
                pop_file=pop_file,
                pops=["MSL", "GBR", "Vindija"],
                use_AA=True,
                verbose=100000,
                interval=interval,
            )
            fs_proj = fs.project([20, 20, 2])
            fs_proj.to_file(f"spectra/3pop.projected.chrom_{chrom}.region_{ii}.fs")
        except:
            print(f"Empty interval chr{chrom} {interval}")

