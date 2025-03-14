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
import numpy as np

n = 30
sample_sets = [msprime.SampleSet(n, "popA"), msprime.SampleSet(n, "popB")]
demog = msprime.Demography.from_demes(g)

# simulate 200 1Mb regions
L = 1e6
u = 1e-8
r = 1e-8
num_reps = 200

tss = msprime.sim_ancestry(
    samples=sample_sets,
    demography=demog,
    sequence_length=L,
    recombination_rate=r,
    num_replicates=num_reps,
    random_seed=42
)

# get list of replicate spectra
spectra = []
for i, ts in enumerate(tss):
    mts = msprime.sim_mutations(ts, rate=u, random_seed=i+13)
    fs_rep = mts.allele_frequency_spectrum(
        sample_sets=[range(2*n), range(2*n, 4*n)],
        polarised=True,
        span_normalise=False
    )
    spectra.append(fs_rep)

fs = np.sum(spectra, axis=0)
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

By running `fs.Fst()`, we find an FST value of between around 0.065 to 0.07.

### Running inference

We'll first fit a simpler (misspecified) model that doesn't allow for
size changes within populations, and assumes a symmetric migration rate.

Our initial model guess is given by [test1.yaml](test1.yaml):
```YAML
description: split-with-migration model, with constant sizes and symmetric migration
generation_time: 3
time_units: years
demes:
- name: ancestral
  epochs:
  - {start_size: 5000, end_time: 1500}
- name: popA
  ancestors: [ancestral]
  epochs:
  - {start_size: 5000}
- name: popB
  ancestors: [ancestral]
  epochs:
  - {start_size: 5000}
migrations:
- demes: [popA, popB]
  rate: 1e-4
```

To run inference using `moments.Demes.Inference`, we also need the parameters
options file, which may be defined in [test1.options.yaml](test1.options.yaml).
In this file, we specify which parameters to optimize, which values in the
Demes-specified model to fit, and impose any bounds on those parameters:
```YAML
parameters:
- name: Ne
  values:
  - demes:
      ancestral:
        epochs:
          0: start_size
  lower_bound: 500
  upper_bound: 50000
- name: NA
  values:
  - demes:
      popA:
        epochs:
          0: start_size
  lower_bound: 500
  upper_bound: 50000
- name: NB
  values:
  - demes:
      popB:
        epochs:
          0: start_size
  lower_bound: 500
  upper_bound: 50000
- name: T
  values:
  - demes:
      ancestral:
        epochs:
          0: end_time
  lower_bound: 0
  upper_bound: 50000
- name: m
  values:
  - migrations:
      0: rate
  lower_bound: 1e-8
  upper_bound: 1e-2
```

Optimization then can be performed as below:
```python
# The total mutation rate
U = u * L * num_reps
# Run the optimization
moments.Demes.Inference.optimize(
    "test1.yaml",
    "test1.options.yaml",
    fs,
    method="fmin",
    perturb=1,
    uL=U,
    output="test1.fit.yaml"
)
```

The output optimized demes model is stored in [test1.fit.yaml](test1.fit.yaml).
Visualizing this fit:
![test1.fit.png](test1.fit.png)
