import demes
import msprime
import numpy as np
import moments
import pickle
import gzip

_pop_ids = ["popA", "popB"]

def run_sims():
    """
    Returns a list of frequency spectra, over some number of replicates.
    """
    g = demes.load("model.yaml")

    # sample sizes of 30 diploids from each population
    n = 30
    sample_sets = [msprime.SampleSet(n, _pop_ids[0]), msprime.SampleSet(n, _pop_ids[1])]
    demog = msprime.Demography.from_demes(g)

    # simulate 500 1Mb regions
    L = 1e6
    u = 1e-8
    r = 1e-8
    num_reps = 500

    tss = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demog,
        sequence_length=L,
        recombination_rate=r,
        num_replicates=num_reps,
        random_seed=42,
    )

    # get list of replicate spectra
    spectra = np.zeros((num_reps, 2*n+1, 2*n+1))
    for i, ts in enumerate(tss):
        mts = msprime.sim_mutations(ts, rate=u, random_seed=i+13)
        fs_rep = mts.allele_frequency_spectrum(
            sample_sets=[range(2*n), range(2*n, 4*n)],
            polarised=True,
            span_normalise=False
        )
        spectra[i] = fs_rep
    return spectra


def bootstrap_spectra(spectra):
    """
    TODO: complete this docstring
    """
    num_reps = len(spectra)
    bs_spectra = []
    for _ in range(num_reps):
        choices = np.random.choice(num_reps, num_reps, replace=True)
        fs_rep = np.sum(spectra[choices], axis=0)
        bs_spectra.append(moments.Spectrum(fs_rep, pop_ids=_pop_ids))
    return bs_spectra


def main():
    spectra = run_sims()
    fs = np.sum(spectra, axis=0)
    fs = moments.Spectrum(fs, pop_ids=_pop_ids)
    fs.tofile("data/data.fs")

    bs_spectra = bootstrap_spectra(spectra)
    with gzip.open("data/bootstrapped_data.pkl.gz", "wb+") as fout:
        pickle.dump(bs_spectra, fout)

if __name__ == "__main__":
    main()
