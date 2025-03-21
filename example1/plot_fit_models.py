import demes, demesdraw
import matplotlib.pylab as plt
import numpy as np

g1 = demes.load("models_and_fits/model.misspec.out.yaml")
g2 = demes.load("models_and_fits/model.misspecB.out.yaml")
g3 = demes.load("models_and_fits/model.out.yaml")

fig = plt.figure(1, figsize=(12, 4))
ax1 = plt.subplot(1, 3, 1)
ax2 = plt.subplot(1, 3, 2)
ax3 = plt.subplot(1, 3, 3)

demesdraw.tubes(g1, ax=ax1, inf_ratio=0.8, labels="xticks", scale_bar=True)
demesdraw.tubes(g2, ax=ax2, inf_ratio=0.8, labels="xticks", scale_bar=True)
demesdraw.tubes(g3, ax=ax3, labels="xticks", scale_bar=True)

w1 = ax1.get_xlim()[1] - ax1.get_xlim()[0]
w2 = ax2.get_xlim()[1] - ax2.get_xlim()[0]
w3 = ax3.get_xlim()[1] - ax3.get_xlim()[0]
m1 = np.mean(ax1.get_xlim()) 
m2 = np.mean(ax2.get_xlim())
m3 = np.mean(ax3.get_xlim()) 

w = np.max([w1, w2, w3])

ax1.set_xlim((m1 - w/2, m1 + w/2))
ax2.set_xlim((m2 - w/2, m2 + w/2))
ax3.set_xlim((m3 - w/2, m3 + w/2))

ax1.set_ylim(top=ax3.get_ylim()[1])
ax2.set_ylim(top=ax3.get_ylim()[1])

ax1.set_title("Misspecified ($LL\\approx-29000$)")
ax2.set_title("Misspecified ($LL\\approx-21000$)")
ax3.set_title("Re-inferred true model ($LL\\approx-15000$)")

ax1.hlines(75000, ax1.get_xlim()[0] + 1000, ax1.get_xlim()[1] - 1000,
           colors=["k"], linestyles="dashed", lw=1)
ax1.text(ax1.get_xlim()[0] + 20000, 85000, "True split time", fontsize=8, va="center", ha="center")
ax2.hlines(75000, ax2.get_xlim()[0] + 1000, ax2.get_xlim()[1] - 1000,
           colors=["k"], linestyles="dashed", lw=1)
ax2.text(ax2.get_xlim()[0] + 15000, 85000, "True split time", fontsize=8, va="center", ha="center")
ax3.hlines(75000, ax3.get_xlim()[0] + 1000, ax3.get_xlim()[1] - 1000,
           colors=["k"], linestyles="dashed", lw=1)
ax3.text(ax3.get_xlim()[0] + 15000, 85000, "True split time", fontsize=8, va="center", ha="center")


fig.tight_layout()
fig.savefig("model.fits.png", dpi=300)
plt.show()
