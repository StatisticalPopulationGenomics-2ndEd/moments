## Make a plot of expectations from the model fit by ./fit_model_folded.py

import demes
import demesdraw
import matplotlib.pyplot as plt
import moments 


def main():
    # define paths and parameters
    graph_file = "MSL_model.folded_fit.yaml"
    data_file = "../../spectra/MSL"
    u = 1.5e-8
    L = 960914001

    # load data, obtain model expectations and fold them
    graph = demes.load(graph_file)
    model = moments.Demes.SFS(graph, samples={"MSL": 20}, u=u, L=L)
    model = model.fold()
    data = moments.Spectrum.from_file(data_file)
    data = data.fold()

    # plot size history
    demesdraw.size_history(graph)
    plt.savefig("MSL_model.folded_fit.size_history.png")
    plt.close()

    # plot SFS fit
    moments.Plotting.plot_1d_comp_Poisson(
        model, data, residual="linear", out="MSL_model.folded_fit.comp_1d.png"
    )

    return 

    
main() 
