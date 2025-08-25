# bootstrap over site frequency spectra. also take their sum.

import numpy as np
import moments 
import pickle
import pandas
import random
import os 


fnames = sorted(os.listdir("spectra/"))
region_fs = dict()
for fname in fnames: 
    region = ".".join(fname.split(".")[2:4])
    region_fs[region] = moments.Spectrum.from_file("spectra/" + fname)


L_tbl = pandas.read_csv("region_L_tbl.csv")
L_tot = np.sum(L_tbl["L"])


# resample bootstrap replicates.
n_reps = 100
n_samps = len(region_fs)
boot_reps = []
for ii in range(n_reps):
    draws = random.choices(list(region_fs.keys()), k=n_samps)
    boot_fs = 0
    boot_L = 0
    for region in draws:
        boot_fs += region_fs[region]
        boot_L += np.array(L_tbl[L_tbl["region"] == region]["L"])[0]
    boot_reps.append((boot_fs, boot_L))
with open("bootstrap_data.pkl", "wb") as fout:
    pickle.dump(boot_reps, fout)


# take the sum
sum_fs = 0
for region in region_fs:
    sum_fs += region_fs[region]
sum_fs.to_file("MSL_GBR_Vindija.fs")