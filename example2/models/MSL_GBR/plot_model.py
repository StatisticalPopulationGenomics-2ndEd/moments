## Make a plot of expectations from the model fit by ./fit_model.py

import demes
import demesdraw
import matplotlib.pyplot as plt
import moments 


def main():
    # define paths and parameters
    graph_file = "MSL_GBR_model.misid_fit.yaml"
    data_file = "../../spectra/MSL_GBR"
    u = 1.5e-8
    L = 960914001
    p_misid = 0.0308139

    # load data and obtain model expectations (adjusted by fit misid pr.)
    graph = demes.load(graph_file)
    model = moments.Demes.SFS(graph, samples={"MSL": 20, "GBR": 20}, u=u, L=L)
    model = moments.Misc.flip_ancestral_misid(model, p_misid)
    data = moments.Spectrum.from_file(data_file)

    # plot size history
    ax = plt.subplot(111)
    demesdraw.size_history(graph, ax=ax)
    plt.savefig("MSL_GBR_model.misid_fit.size_history.png")
    plt.close()

    # plot marginal fit for MSL
    moments.Plotting.plot_1d_comp_Poisson(
        model.marginalize([1]), 
        data.marginalize([1]), 
        residual="linear", 
        out="MSL_GBR_model.misid_fit.comp_1d_MSL.png"
    )

    # plot marginal fit for GBR
    moments.Plotting.plot_1d_comp_Poisson(
        model.marginalize([0]), 
        data.marginalize([0]), 
        residual="linear", 
        out="MSL_GBR_model.misid_fit.comp_1d_GBR.png"
    )

    # plot joint fit (2d) for MSL, GBR
    moments.Plotting.plot_2d_comp_Poisson(
        model, 
        data, 
        residual="linear", 
        out="MSL_GBR_model.misid_fit.comp_2d.png"
    )

    return 

    
main() 
