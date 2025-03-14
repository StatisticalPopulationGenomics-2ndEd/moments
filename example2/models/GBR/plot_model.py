## Make a plot of expectations from the model fit by ./fit_model.py

import demes
import demesdraw
import matplotlib.pyplot as plt
import moments 


def main():
    # define paths and parameters
    graph_file = "GBR_model.misid_fit.yaml"
    data_file = "../../spectra/GBR"
    u = 1.5e-8
    L = 960914001
    p_misid = 0.0576934

    graph = demes.load(graph_file)
    model = moments.Demes.SFS(graph, samples={"GBR": 20}, u=u, L=L)
    model = moments.Misc.flip_ancestral_misid(model, p_misid)
    data = moments.Spectrum.from_file(data_file)

    # plot size history
    demesdraw.size_history(graph)
    plt.savefig("GBR_model.misid_fit.size_history.png")
    plt.close()

    # plot SFS fit
    moments.Plotting.plot_1d_comp_Poisson(
        model, data, residual="linear", out="GBR_model.misid_fit.comp_1d.png"
    )

    return 

    
main() 
