## Example 2: Inferring a model with empirical data 

In this section, we use moments and its Demes interface to fit a model of archaic admixutre to empirical data. Our data set consists of two modern human populations- the Mende from Sierre Leone (MSL) and British from Great Britain (GBR)- and a single Neandertal individual (from Vindija cave in Croatia). We will model the out-of-Africa migration and subsequent Neandertal admixture with GBR and Vindija as proxies. 

### Estimating the SFS from sequence data

The tools we used to estimate the SFS are available in `data`. We lifted sequences for the Vindija Neandertal from genome build hg19 to hg38 to allow direct comparison to a set of recently-resequenced high-coverage 1000 Genomes sequences (Byrska-Bishop et al, 2022).

### Projecting the estimated SFS down

With close to 100 diploid genomes spread across two populations, the SFS we have estimated is rather large. We may wish to project it down to a smaller size to speed up optimization. We can also marginalize populations out of the SFS, removing them and further reducing its size.

```python
import moments
sfs = moments.Spectrum.from_file("spectra/full_MSL_GBR_Vindija")
sfs.pop_ids
>>> ['MSL', 'GBR', 'Vindija']
sfs.sample_sizes
>>> array([170, 182, 2])

# marginalize to MSL and project to a size of 10 diploids
sfs_msl = sfs.marginalize([1, 2])
sfs_msl_projected = sfs_msl.project([20])
sfs_msl_projected.to_file("spectra/MSL")

# project MSL and GBR down to 10 diploids each
sfs_projected = sfs.project([20, 20, 2])
sfs_projected.to_file("spectra/MSL_GBR_Vindija")
```

Let's plot a comparison between the projected spectra for MSL and GBR:

```python
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.subplot(111)
ax.semilogy(sfs_projected.marginalize([1,2]), "-o", ms=6, lw=1, mfc="w", label="MSL")
ax.semilogy(sfs_projected.marginalize([0,2]), "-o", ms=6, lw=1, mfc="w", label="GBR")
ax.set_xlabel("Derived allele frequency")
ax.set_ylabel("Count")
ax.set_xticks([0, 5, 10, 15, 20])
ax.legend()
plt.savefig("comp_1d_MSL_GBR.png", dpi=244)
```

![Marginal SFS for MSL and GBR](comp_1d_MSL_GBR.png)

Let's also plot a 2d spectrum for MSL and GBR using `moments.Plotting`:

```python
moments.Plotting.plot_single_2d_sfs(sfs_projected.marginalize([2]), out="comp_2d_MSL_GBR.png")
```

![Joint SFS for MSL and GBR](comp_2d_MSL_GBR.png)


### Fitting marginal demographies 

Before fitting the three-population model, we can estimate parameters for single-population models. This should help us to select reasonable topologies and initial parameter values before approaching more complex multi-population models. Let's begin with the MSL population and fit a 3-epoch model. We assume a generation time of 29 years.

Our model: ([MSL_model.yaml](models/MSL/MSL_model.yaml))
```YAML
description: 3-epoch model for the MSL population with piecewise-constant population sizes
time_units: years
generation_time: 29
demes: 
- name: MSL
  epochs: 
  - {start_size: 20000, end_time: 2e5}
  - {start_size: 25000, end_time: 2e4}
  - {start_size: 40000, end_time: 0}
```

We specify these parameters: ([options_MSL_model.yaml](models/MSL/options_MSL.yaml))
```YAML
parameters: 
- name: N_A 
  description: Ancestral population size
  lower_bound: 100
  upper_bound: 50000
  values: 
    - demes: 
        MSL: 
          epochs: 
            0: start_size
- name: T_AMH
  description: Time of first size expansion
  upper_bound: 1e6
  values: 
    - demes: 
        MSL: 
          epochs: 
            0: end_time
- name: N_AMH
  description: Population size following first expansion
  lower_bound: 100
  upper_bound: 100000
  values: 
    - demes: 
        MSL: 
          epochs: 
            1: start_size
- name: T_MSL
  description: Time of second size expansion
  lower_bound: 1e3
  values: 
    - demes: 
        MSL: 
          epochs: 
            1: end_time
- name: N_MSL
  description: Contemporary population size 
  lower_bound: 100
  upper_bound: 100000
  values: 
    - demes: 
        MSL: 
          epochs: 
            2: start_size
constraints: 
- params: [T_AMH, T_MSL]
  constraint: greater_than
```

Let's take $u=1.5\cdot10^{-8}$ as an estimate of the human nucleotide mutation rate per generation. The optimization function in `moments.Demes` takes the compound parameter `uL`, so we must also know the length of sequence from which we estimated the SFS (here, $L=960,914,001$). We may choose to fit either the folded or unfolded SFS; when fitting an unfolded SFS, we can incorporate error in the assigment of ancestral nucleotide states into our inference by simaltaneously estimating a 'misid' probability. 

```python
import moments

# define paths and parameters
graph_file = "models/MSL/MSL_model.yaml"
options_file = "models/MSL/options_MSL.yaml"
data_file = "spectra/MSL"
output = "models/MSL/MSL_model.misid_fit.yaml"
u = 1.5e-8
L = 960914001
misid_guess = 0.03

data = moments.Spectrum.from_file(data_file)

# fit using the Powell algorithm
param_names, fit_params, ll = moments.Demes.Inference.optimize(
    graph_file, 
    options_file, 
    data, 
    maxiter=10000,
    fit_ancestral_misid=True, 
    misid_guess=misid_guess,
    uL=u * L, 
    verbose=50,
    output=output, 
    method="powell"
)

# print results
print("Best-fit log-likelihood:", -ll)
print("Best-fit parameters:")
for name, val in zip(param_names, fit_params):
    print(f"{name}\t{val:.3}")
```

This code will print convergence messages and the final best-fit parameters. After fitting, we can create plots showing the size history we've inferred and a comparison of the empirical and expected SFS.

```python
import demes 
import demesdraw 
import matplotlib.pyplot as plt

p_misid = 0.0318905

# load data and obtain model expectations (adjusted by fit misid pr.)
graph = demes.load(graph_file)
model = moments.Demes.SFS(graph, samples={"MSL": 20}, u=u, L=L)
model = moments.Misc.flip_ancestral_misid(model, p_misid)

# plot size history
demesdraw.size_history(graph)
plt.savefig("MSL_model.misid_fit.size_history.png")
plt.close()

# plot SFS fit
moments.Plotting.plot_1d_comp_Poisson(
    model, data, residual="linear", out="MSL_model.misid_fit.comp_1d.png"
)
```

![Inferred MSL size history](models/MSL/MSL_model.misid_fit.size_history.png)
![Empirical and best-fit SFS for MSL](models/MSL/MSL_model.misid_fit.comp_1d.png)


### Fitting a joint demography


### Iteratively introducing complexity



### Computing confidence invervals


### References