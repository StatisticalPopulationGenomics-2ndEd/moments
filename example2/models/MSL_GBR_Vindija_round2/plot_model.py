## Make model/data comparison plots for the model fit by ./fit_model.py

import demes
import demesdraw 
import matplotlib.pyplot as plt
import moments


def main():

    graph_file = "MSL_GBR_Vindija_round2.misid_fit.yaml"
    data_file = "../../spectra/MSL_GBR_Vindija"
    u = 1.5e-8
    L = 960914001
    p_misid = 0.0218

    graph = demes.load(graph_file)
    data = moments.Spectrum.from_file(data_file)
    model = moments.Demes.SFS(
        graph, 
        samples={"MSL": 20, "GBR": 20, "Vindija": 2}, 
        u=u,
        L=L
    )
    model = moments.Misc.flip_ancestral_misid(model, p_misid)

    print(moments.Inference.ll(model, data))

    # plot demography in linear and log time
    ax = plt.subplot(111)
    demesdraw.tubes(graph, log_time=True, ax=ax)
    plt.savefig("MSL_GBR_Vindija_round2.misid_fit.tubes_logtime.png")
    plt.close()
    ax = plt.subplot(111)
    demesdraw.tubes(graph, ax=ax)
    plt.savefig("MSL_GBR_Vindija_round2.misid_fit.tubes.png")
    plt.close()

    # plot marginal fit for GBR
    moments.Plotting.plot_1d_comp_Poisson(
        model.marginalize([0, 2]), 
        data.marginalize([0, 2]),
        residual="linear", 
        out="MSL_GBR_Vindija_round2.misid_fit.comp_1d_GBR.png"
    )
    # and for MSL
    moments.Plotting.plot_1d_comp_Poisson(
        model.marginalize([1, 2]), 
        data.marginalize([1, 2]),
        residual="linear", 
        out="MSL_GBR_Vindija_round2.misid_fit.comp_1d_MSL.png"
    )

    # plot 3d comparison
    moments.Plotting.plot_3d_comp_Poisson(
        model, 
        data, 
        residual="linear", 
        out="MSL_GBR_Vindija_round2.misid_fit.comp_3d.png"
    )

    return


if __name__ == "__main__":
    main()

