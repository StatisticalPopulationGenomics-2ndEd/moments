import demes, demesdraw
import matplotlib, matplotlib.pylab as plt
import numpy as np
import moments

# set font sizes
plt.rcParams["legend.title_fontsize"] = "xx-small"
matplotlib.rc("xtick", labelsize=5)
matplotlib.rc("ytick", labelsize=5)
matplotlib.rc("axes", labelsize=6)
matplotlib.rc("axes", titlesize=6)
matplotlib.rc("legend", fontsize=5)


def plot_demes_model(g, ax, inf_ratio=0.2, top_lim=None, ylabel=None, title=None):
    demesdraw.tubes(g, ax=ax, labels="xticks", inf_ratio=inf_ratio)
    if top_lim is not None:
        ax.set_ylim(top=top_lim)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

def plot_1d_sfs(data, ax):
    pass

def plot_2d_sfs(data, ax):
    pass

def plot_2d_resid(data, model, ax):
    pass


grid = (4, 9)

g = demes.load("models_and_fits/model.yaml")
g1 = demes.load("models_and_fits/model.misspec.out.yaml")
g2 = demes.load("models_and_fits/model.misspecB.out.yaml")
g3 = demes.load("models_and_fits/model.out.yaml")

data = moments.Spectrum.from_file("data/data.fs")

model1 = moments.Demes.SFS(g1, samples={"popA": 60, "popB": 60}, u=500*1e-8*1e6)
resid1 = moments.Inference.Anscombe_Poisson_residual(model1, data)

model2 = moments.Demes.SFS(g2, samples={"popA": 60, "popB": 60}, u=500*1e-8*1e6)
resid2 = moments.Inference.Anscombe_Poisson_residual(model2, data)

model3 = moments.Demes.SFS(g3, samples={"popA": 60, "popB": 60}, u=500*1e-8*1e6)
resid3 = moments.Inference.Anscombe_Poisson_residual(model3, data)

fig = plt.figure(1, figsize=(6.5, 5), dpi=300)
fig.clf()

ax1 = plt.subplot2grid(grid, (0, 0), rowspan=2, colspan=3)
ax2 = plt.subplot2grid(grid, (0, 3), colspan=2)
ax3 = plt.subplot2grid(grid, (1, 3), colspan=2)
ax4 = plt.subplot2grid(grid, (0, 5), colspan=2)
ax5 = plt.subplot2grid(grid, (1, 5), colspan=2)
ax6 = plt.subplot2grid(grid, (0, 7), colspan=2)
ax7 = plt.subplot2grid(grid, (1, 7), colspan=2)
ax8 = plt.subplot2grid(grid, (2, 0), rowspan=2, colspan=3)
ax9 = plt.subplot2grid(grid, (2, 3), rowspan=2, colspan=3)
ax10 = plt.subplot2grid(grid, (2, 6), rowspan=2, colspan=3)

# plot data
moments.Plotting.plot_1d_fs(data.marginalize([1]), ax=ax2, show=False, ms=1, lw=0.5)
ax2.set_title("Marginal SFS (popA)")
moments.Plotting.plot_1d_fs(data.marginalize([0]), ax=ax3, show=False, ms=1, lw=0.5)
ax3.set_title("Marginal SFS (popB)")
moments.Plotting.plot_single_2d_sfs(data, ax=ax4, show=False)
ax4.set_title("Joint SFS")

# plot residuals
moments.Plotting.plot_2d_resid(resid1, ax=ax6, show=False, resid_range=10)
ax6.set_title("Residual (panel H)")
moments.Plotting.plot_2d_resid(resid2, ax=ax5, show=False, resid_range=10)
ax5.set_title("Residual (panel I)")
moments.Plotting.plot_2d_resid(resid3, ax=ax7, show=False, resid_range=10)
ax7.set_title("Residual (panel J)")

# plot models
plot_demes_model(g, ax1, title="Simulated model")
plot_demes_model(g1, ax8, top_lim=ax1.get_ylim()[1], inf_ratio=0.8,
                 title="Misspec. model ($LL\\approx-29000$)")
plot_demes_model(g2, ax9, top_lim=ax1.get_ylim()[1], ylabel="", inf_ratio=0.8,
                 title="Misspec. model ($LL\\approx-21000$)")
plot_demes_model(g3, ax10, top_lim=ax1.get_ylim()[1], ylabel="",
                 title="Reinferred model ($LL\\approx-15000$)")

# adjust width of plotted deme sizes
w1 = ax8.get_xlim()[1] - ax8.get_xlim()[0]
w2 = ax9.get_xlim()[1] - ax9.get_xlim()[0]
w3 = ax10.get_xlim()[1] - ax10.get_xlim()[0]
m1 = np.mean(ax8.get_xlim())
m2 = np.mean(ax9.get_xlim())
m3 = np.mean(ax10.get_xlim())

w = np.max([w1, w2, w3])

ax8.set_xlim((m1 - w/2, m1 + w/2))
ax9.set_xlim((m2 - w/2, m2 + w/2))
ax10.set_xlim((m3 - w/2, m3 + w/2))

# annotate true split time
ax8.hlines(75000, ax8.get_xlim()[0] + 1000, ax8.get_xlim()[1] - 1000,
           colors=["k"], linestyles="dashed", lw=0.75)
ax8.text(ax8.get_xlim()[0] + 20000, 85000, "True split time", fontsize=5, va="center", ha="center")
ax9.hlines(75000, ax9.get_xlim()[0] + 1000, ax9.get_xlim()[1] - 1000,
           colors=["k"], linestyles="dashed", lw=0.75)
ax9.text(ax9.get_xlim()[0] + 15000, 85000, "True split time", fontsize=5, va="center", ha="center")
ax10.hlines(75000, ax10.get_xlim()[0] + 1000, ax10.get_xlim()[1] - 1000,
           colors=["k"], linestyles="dashed", lw=0.75)
ax10.text(ax10.get_xlim()[0] + 15000, 85000, "True split time", fontsize=5, va="center", ha="center")


fig.tight_layout()

fig.text(0.05, 0.98, "A", fontsize=7, va="center", ha="center")
fig.text(0.35, 0.98, "B", fontsize=7, va="center", ha="center")
fig.text(0.35, 0.73, "C", fontsize=7, va="center", ha="center")
fig.text(0.58, 0.98, "D", fontsize=7, va="center", ha="center")
fig.text(0.58, 0.73, "F", fontsize=7, va="center", ha="center")
fig.text(0.79, 0.98, "E", fontsize=7, va="center", ha="center")
fig.text(0.79, 0.73, "G", fontsize=7, va="center", ha="center")
fig.text(0.05, 0.48, "H", fontsize=7, va="center", ha="center")
fig.text(0.35, 0.48, "I", fontsize=7, va="center", ha="center")
fig.text(0.67, 0.48, "J", fontsize=7, va="center", ha="center")

fig.text(0.75, 0.95, "Count", fontsize=6, va="center", ha="center")

fig.savefig("fig1.pdf")
#fig.show()

