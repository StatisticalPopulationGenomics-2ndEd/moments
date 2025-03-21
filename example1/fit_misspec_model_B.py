import moments
import gzip
import pickle

# Load the data
data = moments.Spectrum.from_file("data/data.fs")

# Specify model paths
g_in = "model.misspecB.yaml"
options = "model.misspecB.options.yaml"
g_out = "model.misspecB.out.yaml"

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
    overwrite=True,
    verbose=1,
)
# ret stores parameter names, optimal values, and the log-likelihood
params, vals, ll = ret
print("Param\tBest fit value")
for p, v in zip(params, vals):
    print(f"{p}\t{v}")

# Compute CIs
with gzip.open("data/bootstrapped_data.pkl.gz", "rb") as fin:
    bs_data = pickle.load(fin)

uncerts_FIM = moments.Demes.Inference.uncerts(
    g_out,
    options,
    data,
    uL=U,
    method="FIM",
)

uncerts_GIM = moments.Demes.Inference.uncerts(
    g_out,
    options,
    data,
    bootstraps=bs_data,
    uL=U,
    method="GIM",
)

print("Using FIM:")
print("Param\tBest fit value\t\tstd. err.")
for p, v, u in zip(params, vals, uncerts_FIM):
    print(f"{p}\t{v}\t{u}")

print("Using GIM")
print("Param\tBest fit value\t\tstd. err.")
for p, v, u in zip(params, vals, uncerts_GIM):
    print(f"{p}\t{v}\t{u}")

