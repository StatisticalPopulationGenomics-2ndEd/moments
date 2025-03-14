## parse the SFS from a .vcf file with ancestral state annotations, then save it

import argparse
from collections import defaultdict
import moments 
import numpy as np
import pickle 


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--vcf_file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--population_file", type=str, required=True
    )
    parser.add_argument(
        "-o", "--out_file", type=str, required=True
    )

    return parser.parse_args()


def compute_sfs(data_dict, pop_ids, sample_sizes, mask_corners=True):
    """
    copied from moments
    """
    Npops = len(pop_ids)
    fs = np.zeros(np.asarray(sample_sizes) + 1)
    for snp, snp_info in data_dict.items():
        # Skip SNPs that aren't biallelic.
        if len(snp_info["segregating"]) != 2:
            continue

        allele1, allele2 = snp_info["segregating"]
        if (
            "outgroup_allele" in snp_info
            and snp_info["outgroup_allele"] != "-"
            and snp_info["outgroup_allele"] in snp_info["segregating"]
        ):
            # Otherwise we need to check that it's a useful outgroup
            outgroup_allele = snp_info["outgroup_allele"]
        else:
            # If we're polarized and we didn't have good outgroup info, skip
            # this SNP.
            continue

        # Extract the allele calls for each population.
        allele1_calls = [snp_info["calls"][pop][0] for pop in pop_ids]
        allele2_calls = [snp_info["calls"][pop][1] for pop in pop_ids]
        # assume all chromosomes were called
        # How many chromosomes did we call successfully in each population?
        successful_calls = [
            a1 + a2 for (a1, a2) in zip(allele1_calls, allele2_calls)
        ]
        assert np.all(successful_calls == sample_sizes)

        # Which allele is derived (different from outgroup)?
        if allele1 == outgroup_allele:
            derived_calls = allele2_calls
        elif allele2 == outgroup_allele:
            derived_calls = allele1_calls
        fs[tuple(derived_calls)] += 1

    fsout = moments.Spectrum(fs, mask_corners=mask_corners, pop_ids=pop_ids)
    assert np.all(fsout >= 0)

    return fsout


def main():

    args = get_args()
    vcf_file = args.vcf_file
    pop_file = args.population_file
    out_file = args.out_file
    pop_dict = defaultdict(int)
    with open(pop_file, "r") as fin:
        for line in fin:
            pop = line.split()[1]
            pop_dict[pop] += 2
    pops = list(pop_dict.keys())
    nums = [pop_dict[pop] for pop in pops]
    data_dict = moments.Misc.make_data_dict_vcf(vcf_file, pop_file)
    sfs = compute_sfs(data_dict, pops, nums)
    sfs.to_file(out_file)

    return


if __name__ == "__main__":
    main()

