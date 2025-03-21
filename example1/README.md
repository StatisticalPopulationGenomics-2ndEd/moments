## Example 1: Reinferring simulated parameters

Here, we use msprime to simulate a 2-population split-with-migration model,
in which an ancestral parental population splits some time in the past, and
the descendent populations change sizes over time and exchange migrants.

### The ground-truth demographic model

The input model ([model.yaml](models_and_fits/model.yaml)), which is loosely
based on inferences of human history (see example 2):
```YAML
generation_time: 29
time_units: years
demes:
- name: ancestral
  epochs:
  - {start_size: 15000, end_time: 300000}
  - {start_size: 25000, end_time: 75000}
- name: popA
  ancestors: [ancestral]
  epochs:
  - {start_size: 30000}
- name: popB
  ancestors: [ancestral]
  epochs:
  - {start_size: 1200, end_size: 14000}
migrations:
- demes: [popA, popB]
  rate: 5e-5
```

Using the following to load and plot the model:
```python
import demes, demesdraw, matplotlib.pylab as plt
import moments

g = demes.load("models_and_fits/model.yaml")
fig = plt.figure(figsize=(5, 4))
ax = plt.subplot(1, 1, 1)
demesdraw.tubes(g, ax=ax)
fig.tight_layout()
fig.savefig("model.png")
```

![The input demes model](model.png)

### Simulating data with msprime

We suppose we sample 30 individuals from each population A and B.
We'll simulate 500 1Mb sequences, aggregating across replicates
to construct the SFS and to approximate the process of creating
bootstrapped replicate spectra.


Visualize the marginal spectra:
```python
data = moments.Spectrum.from_file("data/data.fs")

moments.Plotting.plot_1d_comp_Poisson(
    data.marginalize([1]),
    data.marginalize([0]),
    labels=data.pop_ids,
    out="marginal_spectra.png",
)
```
![Marginal SFS](marginal_spectra.png)

By running `fs.Fst()`, we find an FST value of between around 0.065 to 0.07.

### Running inference

We'll first fit a two simpler (misspecified) models that don't allow for
size changes within the ancestral populations.

The first model assumes constant sizes in each of the three demes (the
ancestral population and the two sampled populations). Our initial model guess
is given by [model.misspec.yaml](models_and_fits/model.misspec.yaml):
```YAML
description: split-with-migration model, with constant sizes and symmetric migration
generation_time: 29
time_units: years
demes:
- name: ancestral
  epochs:
  - {start_size: 10000, end_time: 100000}
- name: popA
  ancestors: [ancestral]
  epochs:
  - {start_size: 20000}
- name: popB
  ancestors: [ancestral]
  epochs:
  - {start_size: 20000}
migrations:
- demes: [popA, popB]
  rate: 1e-4
```

To run inference using `moments.Demes.Inference`, we also need the parameters
options file, which here is defined in
[model.misspec.options.yaml](models_and_fits/model.misspec.options.yaml).
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
  upper_bound: 500000
- name: NA
  values:
  - demes:
      popA:
        epochs:
          0: start_size
  lower_bound: 500
  upper_bound: 500000
- name: NB
  values:
  - demes:
      popB:
        epochs:
          0: start_size
  lower_bound: 500
  upper_bound: 500000
- name: T
  values:
  - demes:
      ancestral:
        epochs:
          0: end_time
  lower_bound: 0
  upper_bound: 500000
- name: m
  values:
  - migrations:
      0: rate
  lower_bound: 1e-8
  upper_bound: 1e-2
```

Optimization then can be performed using the `moments.Demes.Inference.optimize`
function:
```python
import moments
import gzip
import pickle

# Load the data
data = moments.Spectrum.from_file("data/data.fs")

# Specify model paths
g_in = "models_and_fits/model.misspec.yaml"
options = "models_and_fits/model.misspec.options.yaml"
g_out = "models_and_fits/model.misspec.out.yaml"

# Parameters from simulated data
u = 1e-8
num_reps = 500
L = 1e6
U = u * L * num_reps


# Run fit
ret = moments.Demes.Inference.optimize(
    g_in,
    options,
    data,
    method="fmin",
    perturb=1,
    uL=U,
    output=g_out,
    overwrite=True
)

# ret stores parameter names, optimal values, and the log-likelihood
params, vals, ll = ret
print("Param\tBest fit value")
for p, v in zip(params, vals):
    print(f"{p}\t{v}")
```
This prints
```
Param	Best fit value
Ne	15324.101946214449
NA	28453.36696520237
NB	5922.663260471864
T	181743.50274735267
m	6.47082934291272e-05
```

The output optimized demes model is stored in
[model.misspec.out.yaml](models_and_fits/model.misspec.out.yaml).

We similarly perform inference for a model that allows for exponential size
changes in `popB` while assuming constant size changes in the ancestor, and
for the simulated model in which we reinfer parameters assuming no model
misspecification. The scripts to run these optimizations are found in this
directory.

To compare model fits between misspecified and correctly specified models,
we use `demesdraw` to plot output models. The plotting script is also found
in this directory.
![Fit models](model.fits.png)

### Confidence intervals
