## Example 1: Reinferring simulated parameters

Here, we use msprime to simulate a 2-population split-with-migration model,
in which an ancestral parental population splits some time in the past, and
the descendent populations change sizes over time and exchange migrants.

### The ground-truth demographic model

The input model ([model.yaml](model.yaml)):
```YAML
description: split-with-migration model for a taxon with an average age of reproduction of three
generation_time: 3
time_units: years
demes:
- name: ancestral
  epochs:
  - {start_size: 5000, end_time: 2500}
  - {start_size: 8000, end_time: 1200}
- name: popA
  ancestors: [ancestral]
  epochs:
  - {start_size: 2000, end_size: 10000}
- name: popB
  ancestors: [ancestral]
  epochs:
  - {start_size: 6000, end_size: 1000}
migrations:
- source: popA
  dest: popB
  rate: 1e-4
- source: popB
  dest: popA
  rate: 3e-4
```

Using the following to load and plot the model:
```python
import demes, demesdraw, matplotlib.pylab as plt
g = demes.load("model.yaml")
fig = plt.figure(figsize=(5, 4))
ax = plt.subplot(1, 1, 1)
demesdraw.tubes(g, ax=ax)
fig.tight_layout()
fig.savefig("model.png")
```

![The input demes model](model.png)

### Simulating data with msprime

Sample 30 individuals from each population A and B.
We'll simulate 100 1M sequences, aggregating across replicates
to construct the SFS.
```python
import msprime

n = 30
sample_sets = [msprime.SampleSet(n, "popA"), msprime.SampleSet(n, "popB")]
demog = msprime.Demography.from_demes(g)

L = 1e6
u = 1e-8
r = 1e-8

tss = msprime.sim_ancestry(
    samples=sample_sets,
    demography=demog,
    sequence_length=L,
    recombination_rate=r,
    num_replicates=100,
    random_seed=42
)

fs = np.zeros((2*n+1, 2*n+1))
for i, ts in enumerate(tss):
    mts = msprime.sim_mutations(ts, rate=u, random_seed=i+13)
    fs_rep = mts.allele_frequency_spectrum(
        sample_sets=[range(2*n), range(2*n, 4*n)],
        polarised=True,
        span_normalise=False
    )
    fs += fs_rep
```

Visualize the marginal spectra
```python
import moments

fs = moments.Spectrum(fs, pop_ids=["popA", "popB"])

moments.Plotting.plot_1d_comp_Poisson(
    fs.marginalize([1]),
    fs.marginalize([0]),
    labels=fs.pop_ids,
    out="marginal_spectra.png",
)
```
![Marginal SFS](marginal_spectra.png)

By running `fs.Fst()`, we find an FST value of around 0.068.

### Running inference



