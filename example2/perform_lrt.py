# Perform the LRT between MSL_GBR_Vindija_null_model and MSL_GBR_Vindija_step2

import demes
import moments
import numpy as np
import pickle


# parameters
u = 1.5e-8
L = 960914001
U = u * L
samples = {"MSL": 20, "GBR": 20, "Vindija": 2}


# models
model_0_file = "models/MSL_GBR_Vindija_null_model/MSL_GBR_Vindija_null_model.misid_fit.yaml"
misid_0 = 0.0184
model_1_file =  "models/MSL_GBR_Vindija_step2/MSL_GBR_Vindija_model_step2.misid_fit.yaml"
misid_1 = 0.0223


# data
data_file = "spectra/MSL_GBR_Vindija.fs"
data = moments.Spectrum.from_file(data_file)

bootstrap_data_file = "data/bootstrap_data.pkl"
with open(bootstrap_data_file, "rb") as fin:
    all_boot = pickle.load(fin)
all_boot = [x[0] for x in all_boot]


graph_0 = demes.load(model_0_file)
model_0 = moments.Demes.SFS(graph_0, samples=samples, u=u, L=L)
model_0 = moments.Misc.flip_ancestral_misid(model_0, misid_0)
ll_0 = moments.Inference.ll(model_0, data)

graph_1 = demes.load(model_1_file)
model_1 = moments.Demes.SFS(graph_1, samples=samples, u=u, L=L)
model_1 = moments.Misc.flip_ancestral_misid(model_1, misid_1)
ll_1 = moments.Inference.ll(model_1, data)


# load parameter values from nested model
options_0_file = "models/MSL_GBR_Vindija_null_model/MSL_GBR_Vindija_null_model_options.yaml"
options_0 = moments.Demes.Inference._get_params_dict(options_0_file)
builder_0 = demes.load(model_0_file).asdict()
pnames, p0, *_ = moments.Demes.Inference._set_up_params_and_bounds(options_0, builder_0)

nested_idx = np.array([3])
adj_params = np.insert(p0, 3, 0.0)


options_1_file = "models/MSL_GBR_Vindija_step2/MSL_GBR_Vindija_options_step2.yaml"
builder_1 = demes.load_asdict(model_1_file)
options_1 = moments.Demes.Inference._get_params_dict(options_1_file)
pnames_1, p1, *_ = moments.Demes.Inference._set_up_params_and_bounds(options_1, builder_1)

def model_func(params, ns):
    global builder_1
    builder_1 = moments.Demes.Inference._update_builder(builder_1, options_1, params)
    graph = demes.Graph.fromdict(builder_1)
    model = moments.Demes.SFS(graph, samples=samples, u=u, L=L)
    model = moments.Misc.flip_ancestral_misid(model, misid_1)
    return model


adj = moments.Godambe.LRT_adjust(
    model_func, all_boot, adj_params, data, nested_idx, multinom=False)


delta_ll = 2 * (ll_1 - ll_0)
print(f"ll_0 = {ll_0}")
print(f"ll_1 = {ll_1}")
print(f"delta_ll = {delta_ll}")
print(f"adjustment = {adj}")
delta_ll_adj = adj * delta_ll
print(f"delta_ll_adj = {delta_ll_adj}")
pval = moments.Godambe.sum_chi2_ppf(delta_ll_adj, (0.5, 0.5))
print(f"adjusted p-value = {pval}")
