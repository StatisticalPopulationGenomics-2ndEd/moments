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

```python
import msprime

sample_sets = [msprime.SampleSet(), msprime.SampleSet()]
```
