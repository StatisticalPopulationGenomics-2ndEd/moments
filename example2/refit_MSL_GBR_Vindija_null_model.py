## Fit parameters to the Neandertal branch of a three-deme MSL, GBR, Vindija model
## Then compute confidence intervals and plot the fit

import demes, demesdraw
import moments
import matplotlib.pyplot as plt


# define model paths
graph_file =  "models/MSL_GBR_Vindija_null_model/MSL_GBR_Vindija_null_model.misid_fit.yaml"
options_file = "models/MSL_GBR_Vindija_null_model/MSL_GBR_Vindija_null_model_options.yaml"

# define data path and load data
data_file = "spectra/MSL_GBR_Vindija.fs"
data = moments.Spectrum.from_file(data_file)

# define output graph filename
output_file = "models/MSL_GBR_Vindija_null_model/MSL_GBR_Vindija_null_model.misid_refit.yaml"

# define parameters
u = 1.5e-8
L = 960914001
U = u * L
misid_guess = 0.0184

# fit using the Powell algorithm
param_names, fit_params, ll = moments.Demes.Inference.optimize(
    graph_file, 
    options_file, 
    data, 
    maxiter=10000,
    fit_ancestral_misid=True, 
    misid_guess=misid_guess,
    uL=U, 
    verbose=1,
    output=output_file, 
    method="powell",
    overwrite=True
)

# print results
print("Best-fit log-likelihood:", -ll)
print("Best-fit parameters:")
for name, value in zip(param_names, fit_params):
    print(f"{name}\t{value:.3}")

# retrieve p_misid
p_misid = fit_params[-1]

# compute confidence intervals using FIM 
log_file = "models/MSL_GBR_Vindija_null_model/MSL_GBR_Vindija_null_model.misid_refit.uncerts.txt" 
uncerts = moments.Demes.Inference.uncerts(
    output_file,
    options_file,
    data,
    uL=U,
    method="FIM",
    fit_ancestral_misid=True,
    misid_fit=p_misid,
    verbose=10,
    output=log_file,
    overwrite=True,
    eps=1e-2
)

# print confidence intervals
print("95% confidence intervals:")
print("param\t2.5%\t97.5%")
for name, value, err in zip(param_names, fit_params, uncerts):
    print(f"{name}\t{value - 1.96 * err:.3}\t{value + 1.96 * err:.3}")
