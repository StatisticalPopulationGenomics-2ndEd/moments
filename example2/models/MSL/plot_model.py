## Make a plot of expectations from the model fit by ./fit_model.py

import demes
import demesdraw
import matplotlib.pyplot as plt
import moments 


def main():
    # define paths and parameters
    graph_file = "MSL_model.misid_fit.yaml"
    data_file = "../../spectra/MSL"
    u = 1.5e-8
    L = 960914001
    p_misid = 0.0318905

    # load data and obtain model expectations (adjusted by fit misid pr.)
    graph = demes.load(graph_file)
    model = moments.Demes.SFS(graph, samples={"MSL": 20}, u=u, L=L)
    model = moments.Misc.flip_ancestral_misid(model, p_misid)
    data = moments.Spectrum.from_file(data_file)

    # plot size history
    ax = plt.subplot(111)
    demesdraw.size_history(graph, ax=ax)
    plt.savefig("MSL_model.misid_fit.size_history.png")
    plt.close()

    # plot SFS fit
    moments.Plotting.plot_1d_comp_Poisson(
        model, data, residual="linear", out="MSL_model.misid_fit.comp_1d.png"
    )

    return 

    
main() 
