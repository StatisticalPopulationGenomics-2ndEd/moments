import moments
import gzip
import pickle

# Load the data
data = moments.Spectrum.from_file("data/data.fs")

# Specify model paths
g_in = "model.yaml"
options = "model.options.yaml"
g_out = "model.out.yaml"

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
    verbose=1,
    method="fmin",
    perturb=1,
    uL=U,
    output=g_out,
    overwrite=True
)

# Compute CIs
with gzip.open("data/bootstrapped_data.pkl.gz", "rb") as fin:
    bs_data = pickle.load(fin)

uncerts = moments.Demes.Inference.uncerts(
    g_out,
    options,
    data,
    bootstraps=bs_data,
    uL=U,
    verbose=1,
    method="GIM",
)
