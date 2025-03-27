
import demes
import demesdraw
import matplotlib
import matplotlib.pylab as plt
import moments
import numpy as np


plt.rcParams["legend.title_fontsize"] = "xx-small"
matplotlib.rc("xtick", labelsize=5)
matplotlib.rc("ytick", labelsize=5)
matplotlib.rc("axes", labelsize=6)
matplotlib.rc("axes", titlesize=6)
matplotlib.rc("legend", fontsize=5)


def plot_2d_comp(ax, model, data):

    ax.semilogy(data, "-o", ms=4, lw=1, mfc="w", label="Data")
    ax.semilogy(model, "-o", ms=2, lw=1, label="Model")        
    ax.set_xlim(0, data.sample_sizes[0])
    ax.set_xlabel("Allele frequency")
    ax.set_ylabel("Count")
    ax.legend()

    return


def main():

    graph = demes.load(
        "models/MSL_GBR_Vindija_round3/MSL_GBR_Vindija_round3.misid_fit.yaml"
    )
    data = moments.Spectrum.from_file("spectra/MSL_GBR_Vindija")
    p_misid = 0.0223
    u = 1.5e-8
    L = 960914001
    model = moments.Demes.SFS(
        graph, samples={"MSL": 20, "GBR": 20, "Vindija": 2}, u=u, L=L
    )
    model = moments.Misc.flip_ancestral_misid(model, p_misid)
    resid = moments.Inference.Anscombe_Poisson_residual(model, data)

    grid = (4, 3)
    fig = plt.figure(1, figsize=(6.5, 5.5))
    ax0 = plt.subplot2grid(grid, (0, 0), rowspan=2, colspan=2)
    ax1 = plt.subplot2grid(grid, (0, 2), rowspan=1, colspan=1)
    ax2 = plt.subplot2grid(grid, (1, 2), rowspan=1, colspan=1)
    ax3 = plt.subplot2grid(grid, (2, 0), rowspan=1, colspan=1)
    ax4 = plt.subplot2grid(grid, (2, 1), rowspan=1, colspan=1)
    ax5 = plt.subplot2grid(grid, (2, 2), rowspan=1, colspan=1)
    ax6 = plt.subplot2grid(grid, (3, 0), rowspan=1, colspan=1)
    ax7 = plt.subplot2grid(grid, (3, 1), rowspan=1, colspan=1)
    ax8 = plt.subplot2grid(grid, (3, 2), rowspan=1, colspan=1)

    # plot the model
    # custom position mapping
    pos = {
        "Vindija": 5000, 
        "N": 12000,
        "NI": 19000,
        "A": 33000,
        "GBR": 45000,
        "AMH": 60000,
        "MSL": 90000
    }
    demesdraw.tubes(graph, ax=ax0, labels=None, positions=pos)
    label_me = ["Vindija", "NI", "GBR", "MSL"]
    ax0.set_xticks([pos[x] for x in label_me], label_me)
    ax0.set_title("Best-fit model")

    # plot 1d comparisons
    plot_2d_comp(ax1, model.marginalize([1,2]), data.marginalize([1,2]))
    ax1.set_title("Marginal (MSL)")
    plot_2d_comp(ax2, model.marginalize([0,2]), data.marginalize([0,2]))
    ax2.set_title("Marginal (GBR)")

    # plot 2d comparisons
    moments.Plotting.plot_single_2d_sfs(model.marginalize([2]), ax=ax3)
    ax3.set_title("Model (MSL, GBR)")
    moments.Plotting.plot_single_2d_sfs(model.marginalize([0]), ax=ax4)
    ax4.set_title("Model (MSL, Vindija)")
    moments.Plotting.plot_single_2d_sfs(model.marginalize([1]), ax=ax5)
    ax5.set_title("Model (GBR, Vindija)")

    # plot residuals
    moments.Plotting.plot_2d_resid(resid.marginalize([2]), ax=ax6, show=False)
    ax6.set_title("Residual (MSL, GBR)")
    moments.Plotting.plot_2d_resid(resid.marginalize([0]), ax=ax7, show=False)
    ax7.set_title("Residual (MSL, Vindija)")
    moments.Plotting.plot_2d_resid(resid.marginalize([1]), ax=ax8, show=False)
    ax8.set_title("Residual (GBR, Vindija)")

    fig.tight_layout()

    fig.text(0.05, 0.98, "A", fontsize=7, va="center", ha="center")
    fig.text(0.65, 0.98, "B", fontsize=7, va="center", ha="center")
    fig.text(0.65, 0.73, "C", fontsize=7, va="center", ha="center")
    fig.text(0.05, 0.50, "D", fontsize=7, va="center", ha="center")
    fig.text(0.35, 0.50, "E", fontsize=7, va="center", ha="center")
    fig.text(0.65, 0.50, "F", fontsize=7, va="center", ha="center")
    fig.text(0.05, 0.25, "G", fontsize=7, va="center", ha="center")
    fig.text(0.35, 0.25, "H", fontsize=7, va="center", ha="center")
    fig.text(0.65, 0.25, "I", fontsize=7, va="center", ha="center")

    fig.savefig("fig2.pdf")
    plt.show()

    return


main() 
