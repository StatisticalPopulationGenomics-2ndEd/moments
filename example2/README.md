## Example 2: Inferring a model from empirical data 

In this section, we use moments and its Demes interface to infer a 3-population model of archaic admixture with empirical data. Our data set consists of samples from two modern human populations sequenced in the 1000 genomes project- the Mende from Sierre Leone (MSL) and British from Great Britain (GBR)- and a single Neandertal individual from Vindija cave in Croatia. We construct a simple model of the out-of-Africa expansion of modern humans and the following admixture with Neandertals.

### Estimating the SFS from sequence data

The tools we used to estimate the SFS are available in [data](data/). We lifted sequences for the Vindija Neandertal over from genome build hg19 to hg38 ([data/liftover/](data/liftover/)) to allow a direct comparison to a set of recently-resequenced high-coverage 1000 Genomes sequences (Byrska-Bishop et al, 2022). We apply the 1000 genomes "strict" mask, exclude exons and promoters plus 10-kilobase flanking regions around them, and ignore positions that lack high-confidence ancestral states. We also ignore sites which were not genotyped in the four existing high-coverage archaic human sequences. We are left with approximately 1 billion base pairs of sequence. The scripts we used to perform these operations and estimate the SFS are available in [data/parsing/](data/parsing/) and [data/tools](data/tools/).

### Projecting the estimated SFS down

With close to 100 diploid genomes spread across two populations, the SFS we have estimated is large. We may wish to project it down to a smaller size to speed up optimization and other operations. We can also marginalize populations out of the SFS, removing them and further reducing its size. First we load the SFS and inspect its attributes; 
```python
import moments

sfs = moments.Spectrum.from_file("spectra/full_MSL_GBR_Vindija")
print(sfs.pop_ids)
print(sfs.sample_sizes)
```
```
['MSL', 'GBR', 'Vindija']
array([170, 182, 2])
```
We can use projection and marginalization along with the `to_file` method to save one, two and three-dimensional versions of the large three-population SFS, so that we don't need to repeat these operations every time we fit a model. 
```python
# marginalize to MSL, project it to 10 diploids and save the result
sfs_msl = sfs.marginalize([1, 2])
sfs_msl_projected = sfs_msl.project([20])
sfs_msl_projected.to_file("spectra/MSL")

# project MSL and GBR down to 10 diploids and save the result
sfs_projected = sfs.project([20, 20, 2])
sfs_projected.to_file("spectra/MSL_GBR_Vindija")
```

Let's plot a comparison between the projected spectra for MSL and GBR. 
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
plt.savefig("comp_1d_MSL_GBR.png")
```

![Marginal SFS for MSL and GBR](comp_1d_MSL_GBR.png)

We can also plot a 2d spectrum as a heatmap using `moments.Plotting`:

```python
moments.Plotting.plot_single_2d_sfs(sfs_projected.marginalize([2]), out="comp_2d_MSL_GBR.png")
```

![Joint SFS for MSL and GBR](comp_2d_MSL_GBR.png)

### Fitting marginal demographies 

Before fitting the three-population model, it will be useful to estimate parameters for single-population models. This will help us to select reasonable topologies and initial parameter values when approaching more complex multi-population models. Let's begin with the MSL population and fit a model that proposes three epochs with constant effective population sizes. We assume a generation time of 29 years.

The model ([MSL_model.yaml](models/MSL/MSL_model.yaml)):
```YAML
description: 3-epoch model for the MSL population with piecewise-constant population sizes
time_units: years
generation_time: 29
demes: 
- name: MSL
  epochs: 
  - {start_size: 20000, end_time: 200000}
  - {start_size: 25000, end_time: 20000}
  - {start_size: 40000, end_time: 0}
```

The parameters ([options_MSL_model.yaml](models/MSL/options_MSL.yaml)):
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

Let $u=1.5\cdot10^{-8}$ be our estimate of the per-generation nucleotide mutation rate in humans. Here we will use a Poisson likelihood function rather than a multionomial one, and the object we will fit the model to will be SFS counts rather than frequencies. To do so, the optimization function in `moments.Demes` requires the compound parameter `uL`, so we must also know the length of the sequence from which we estimated the SFS. Here, the number of called sites that passed our filters was $L=960,914,001$. We may choose to fit either the folded or unfolded SFS; when fitting an unfolded SFS, we can incorporate error in the assigment of ancestral nucleotide states into our inference by simaltaneously estimating a `misid` probability.

```python
import moments

# define paths and parameters
graph_file = "models/MSL/MSL_model.yaml"
options_file = "models/MSL/options_MSL.yaml"
data_file = "spectra/MSL"
output = "models/MSL/MSL_model.misid_fit.yaml"
u = 1.5e-8
L = 960914001
U = u * L
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
    uL=U, 
    verbose=50,
    output=output, 
    method="powell"
)
```

We can print the results like so;
```python
# print results
print("Best-fit log-likelihood:", -ll)
print("Best-fit parameters:")
for name, val in zip(param_names, fit_params):
    print(f"{name}\t{val:.3}")
```
```
Best-fit log-likelihood: -356.03751402092166
Best-fit parameters:
N_A     1.43e+04
T_AMH   4.35e+05
N_AMH   2.61e+04
T_MSL   1.5e+04
N_MSL   2.1e+04
p_misid 0.0331
```

After fitting the model, we can create plots showing the size history we've inferred and a comparison of the empirical and expected SFS using some built-in `moments` plotting functions. 

```python
import demes 
import demesdraw 
import matplotlib.pyplot as plt

# p_misid was returned as the last element of `fit_params`
p_misid = fit_params[-1]

# load data and obtain model expectations (adjusted by fit misid pr.)
graph = demes.load(output)
model = moments.Demes.SFS(graph, samples={"MSL": 20}, u=u, L=L)
model = moments.Misc.flip_ancestral_misid(model, p_misid)

# plot size history
fig, ax = plt.subplot(111)
demesdraw.size_history(graph, ax=ax)
plt.savefig("MSL_model.misid_fit.size_history.png")
plt.close()

# plot model expectation against data
moments.Plotting.plot_1d_comp_Poisson(
    model, data, residual="linear", out="MSL_model.misid_fit.comp_1d.png"
)
```

![Inferred MSL size history](models/MSL/MSL_model.misid_fit.size_history.png)
![Empirical and best-fit SFS for MSL](models/MSL/MSL_model.misid_fit.comp_1d.png)


### Fitting a joint demography

We fit a one-population model incorporating a sharp contraction followed by exponential growth for the GBR population in the same manner ([models/GBR/GBR_model.misid_fit.yaml](models/GBR/GBR_model.misid_fit.yaml)). Using the parameters inferred in the marginal MSL and GBR models, we now construct a 2-deme model where MSL and GBR diverge ~60,000 years ago ([models/MSL_GBR/MSL_GBR_model.yaml](models/MSL_GBR/MSL_GBR_model.yaml)).
```YAML
description: Simple OOA model with MSL, GBR demes and ancestral expansion.
time_units: years
generation_time: 29
demes:
- name: A
  epochs: 
  - {end_time: 450000, start_size: 14000}
- name: AMH
  ancestors: [A]
  epochs:
  - {end_time: 60000, start_size: 24000}
- name: MSL 
  ancestors: [AMH]
  epochs:
  - {end_time: 0, start_size: 21000}
- name: GBR
  ancestors: [AMH]
  epochs: 
  - {end_time: 0, start_size: 2000, end_size: 21000}
migrations:
- demes: [MSL, GBR]
  rate: 1e-5
```

We parameterize every feature of the model ([models/MSL_GBR/options_MSL_GBR.yaml](models/MSL_GBR/options_MSL_GBR.yaml)). 
```YAML
parameters:
- name: N_A 
  description: Ancestral population size
  lower_bound: 100
  upper_bound: 100000
  values: 
    - demes: 
        A: 
          epochs: 
            0: start_size
- name: T_EXP
  decription: Time of ancestral size expansion
  upper_bound: 1e6
  values: 
    - demes: 
        A: 
          epochs: 
            0: end_time
- name: N_AMH
  description: Ancestral size following expansion
  lower_bound: 100
  upper_bound: 100000
  values: 
    - demes: 
        AMH: 
          epochs: 
            0: start_size
- name: T_OOA
  description: Divergence time of MSL, GBR
  lower_bound: 30e3
  values: 
    - demes: 
        AMH:  
          epochs: 
            0: end_time
- name: N_MSL
  description: Final size of MSL 
  lower_bound: 100
  upper_bound: 100000
  values: 
    - demes: 
        MSL: 
          epochs: 
            0: start_size
- name: N_OOA
  description: Initial size of GBR 
  lower_bound: 100
  upper_bound: 100000
  values: 
    - demes: 
        GBR: 
          epochs: 
            0: start_size
- name: N_GBR
  description: Final size of GBR
  lower_bound: 100
  upper_bound: 200000
  values: 
    - demes: 
        GBR: 
          epochs: 
            0: end_size
- name: m
  description: Symmetric migration rate between MSL, GBR
  lower_bound: 1e-8
  upper_bound: 1e-3
  values: 
    - migrations: 
        0: rate
constraints:
- params: [T_EXP, T_OOA]
  constraint: greater_than
```

This is fitted using
```python
graph_file = "models/MSL_GBR/MSL_GBR_model.yaml"
options_file = "models/MSL_GBR/options_MSL_GBR.yaml"
data_file = "spectra/MSL_GBR"
output = "models/MSL_GBR/MSL_GBR_model.misid_fit.yaml"

data = moments.Spectrum.from_file(data_file)

# fit using the Powell algorithm
param_names, fit_params, ll = moments.Demes.Inference.optimize(
    graph_file, 
    options_file, 
    data, 
    maxiter=10000,
    fit_ancestral_misid=True, 
    misid_guess=misid_guess,
    uL=U, 
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
```
Best-fit log-likelihood: -6960.3200203990145
Best-fit parameters:
N_A     1.45e+04
T_EXP   4.36e+05
N_AMH   2.35e+04
T_OOA   7.25e+04
N_MSL   2.76e+04
N_OOA   1.2e+03
N_GBR   1.41e+04
m       3.35e-05
p_misid 0.0308
```

We are left with some reasonable parameters and the fit to data plotted below;
```python
p_misid = fit_params[-1]
graph = demes.load(output)

model = moments.Demes.SFS(graph, samples={"MSL": 20, "GBR": 20}, u=u, L=L)
model = moments.Misc.flip_ancestral_misid(model, p_misid)
data = moments.Spectrum.from_file(data_file)

# plot size history
ax = plt.subplot(111)
demesdraw.size_history(graph, ax=ax)
plt.savefig("models/MSL_GBR/MSL_GBR_model.misid_fit.size_history.png")
plt.close()

# plot joint fit (2d) for MSL, GBR
moments.Plotting.plot_2d_comp_Poisson(
    model, 
    data, 
    residual="linear", 
    out="models/MSL_GBR/MSL_GBR_model.misid_fit.comp_2d.png"
)
```
![Size history for MSL_GBR](models/MSL_GBR/MSL_GBR_model.misid_fit.size_history.png)
![2d comparison for MSL_GBR](models/MSL_GBR/MSL_GBR_model.misid_fit.comp_2d.png)


### Introducing a third population

Now that we have a relatively well-fitting model for two modern human populations, we incorporate a Neandertal branch and fit the relevant parameters. It is useful to do this in multiple steps- here we first optimize a few new parameters (Neandertal/African modern human divergence time, effective size of Neandertal populations, divergence time of the sampled and introgressing Neandertals and admixture pulse proportion) to obtain a reasonable fit, then run a second round of optimization to further refine all the parameters we considered above. We fix the time of Neandertal admixture to 50 kya (see Sümer et al. 2025) and use a single size parameter $N_N$ for all Neandertal demes, as we may expect these parameters to be poorly constrained given our small Neandertal sample size. The initial graph is [models/MSL_GBR_Vindija_round1/MSL_GBR_Vindija_model.yaml](models/MSL_GBR_Vindija_round1/MSL_GBR_Vindija_model.yaml) 
```YAML
description: Model of Vindija Neandertal, MSL and GBR lineages with admixture.
time_units: years
generation_time: 29
demes:
- name: A
  epochs:
  - {end_time: 6e5, start_size: 14502.424743309206}
- name: N 
  ancestors: [A]
  epochs: 
  - {start_size: 3000, end_time: 125000}
- name: NI
  ancestors: [N]
  epochs: 
  - {start_size: 3000, end_time: 50000}
- name: Vindija
  ancestors: [N]
  epochs: 
  - {start_size: 3000, end_time: 55000}
- name: AMH
  ancestors: [A]
  epochs:
  - {end_time: 436166.1984650331, start_size: 14502.424743309206}
  - {end_time: 72526.00759273666, start_size: 23525.824400305177}
- name: MSL
  ancestors: [AMH]
  epochs:
  - {end_time: 0, start_size: 27591.878990713663}
- name: GBR
  ancestors: [AMH]
  epochs:
  - {end_time: 0, start_size: 1201.9212244444434, end_size: 14129.839474007915}
migrations:
- demes: [MSL, GBR]
  rate: 3.345205234558448e-05
pulses:
- sources: [NI]
  proportions: [0.02]
  dest: GBR
  time: 50000
```

and the parameters are specified in [models/MSL_GBR_Vindija_round1/options_MSL_GBR_Vindija.yaml](models/MSL_GBR_Vindija_round1/options_MSL_GBR_Vindija.yaml):
```YAML
parameters: 
- name: T_NMH
  description: Split time of Neandertal and AMH
  upper_bound: 1e6
  lower_bound: 440000
  values: 
    - demes: 
        A: 
          epochs: 
            0: end_time     
- name: N_N
  description: Effective size of Neandertal demes
  lower_bound: 100
  values: 
    - demes: 
        N: 
          epochs: 
            0: start_size 
        NI: 
          epochs: 
            0: start_size 
        Vindija: 
          epochs: 
            0: start_size               
- name: T_NI
  description: Divergence time of `NI` and Vindija
  lower_bound: 55000
  values: 
    - demes: 
        N: 
          epochs: 
            0: end_time 
- name: p
  desciprion: Admixture proportion from `NI` to GBR
  lower_bound: 1e-5
  upper_bound: 0.10
  values: 
    - pulses: 
        0:  
          proportions: 0
constraints:
- params: [T_NMH, T_NI]
  constraint: greater_than
```

After optimizing these parameters, we run a second round of optimization using the model inferred in the first round ([models/MSL_GBR_Vindija_round1/MSL_GBR_Vindija_model.misid_fit.yaml](models/MSL_GBR_Vindija_round1/MSL_GBR_Vindija_model.misid_fit.yaml)) as a base. Parameters are specified in [models/MSL_GBR_Vindija_round3/options_MSL_GBR_Vindija_round3.yaml](models/MSL_GBR_Vindija_round3/options_MSL_GBR_Vindija_round3.yaml). Because the divergence time $T_{NI}$ appears poorly constrained, going to its lower bound as seen in [models/MSL_GBR_Vindija_round2/MSL_GBR_Vindija_round2.misid_fit.yaml](models/MSL_GBR_Vindija_round2/MSL_GBR_Vindija_round2.misid_fit.yaml), we fix this parameter to a value supported by the literature (90 kya, see Prüfer et al. 2017).
```python
graph_file = "models/MSL_GBR_Vindija_round3/MSL_GBR_Vindija_round3.misid_fit.yaml"
options_file = "models/MSL_GBR_Vindija_round3/options_MSL_GBR_Vindija_round3.yaml"
data_file = "spectra/MSL_GBR_Vindija"
output = "models/MSL_GBR_Vindija_round3/MSL_GBR_Vindija_round3.misid_fit.yaml"

data = moments.Spectrum.from_file(data_file)

param_names, fit_params, ll = moments.Demes.Inference.optimize(
    graph_file, 
    options_file, 
    data, 
    maxiter=10000,
    fit_ancestral_misid=True, 
    misid_guess=misid_guess,
    uL=U, 
    verbose=1,
    output=output,
    method="lbfgsb"
)

# print results
print("Best-fit log-likelihood:", -ll)
print("Best-fit parameters:")
for name, p in zip(param_names, fit_params):
    print(f"{name}\t{p:.3}")
```
```
Best-fit log-likelihood: -24640.177263672063
Best-fit parameters:
N_A     1.54e+04
T_NMH   5.47e+05
N_N     2.77e+03
p       0.0271
T_EXP   3.79e+05
N_AMH   2.12e+04
T_OOA   7.77e+04
N_MSL   2.98e+04
N_OOA   1.2e+03
N_GBR   1.08e+04
m       4.36e-05
p_misid 0.0223
```

Our final best-fit model is plotted below.

![Best-fit model](models/MSL_GBR_Vindija_round3/MSL_GBR_Vindija_round3.misid_fit.tubes.png)
![Best-fit model fit](models/MSL_GBR_Vindija_round3/MSL_GBR_Vindija_round3.misid_fit.comp_3d.png)

### Computing confidence invervals

We may wish to quantify the uncertainty in our best-fit parameter values. We can do this using the `moments.Demes.Inference.uncerts` function. Here we use the default `FIM` method, which does not take bootstrapped data. As described in example 1, this method underestimates parameter uncertainties. We can then print the estimated 95% confidence intervals about our best-fit parameters. 
```python
p_misid = fit_params[-1]

std_errs = moments.Demes.Inference.uncerts(
    output,
    options_file,
    data,
    uL=U,
    method="FIM",
    fit_ancestral_misid=True,
    misid_fit=p_misid
)

# print results
print(r"95% confidence intervals:")
print("param\t2.5%\t97.5%")
for name, val, err in zip(param_names, fit_params, std_errs):
    print(f"{name}\t{val - 1.96 * err:.3}\t{val + 1.96 * err:.3}")
```
```
95% confidence intervals:
param   2.5%    97.5%
N_A     1.54e+04        1.54e+04
T_NMH   5.46e+05        5.48e+05
N_N     2.76e+03        2.78e+03
p       0.0268  0.0273
T_EXP   3.72e+05        3.86e+05
N_AMH   2.11e+04        2.13e+04
T_OOA   7.72e+04        7.83e+04
N_MSL   2.96e+04        3e+04
N_OOA   1.19e+03        1.22e+03
N_GBR   1.06e+04        1.1e+04
m       4.31e-05        4.41e-05
```

We can see that the `FIM` 95% confidence intervals are very tight. 

### References

Marta Byrska-Bishop et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. Cell 185, no. 18 3426-3440, 2022. 

Kay Prüfer et al. A high-coverage Neandertal genome from Vindija Cave in Croatia. Science 358, no. 6363 655-658, 2017.

Arev Sümer et al. Earliest modern human genomes constrain timing of Neanderthal admixture. Nature 638, no. 8051 711-717, 2025.