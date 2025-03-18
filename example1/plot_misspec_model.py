import demes, demesdraw
import matplotlib.pylab as plt
import numpy as np

g = demes.load("model.yaml")
g2 = demes.load("model.misspec.out.yaml")

fig = plt.figure(figsize=(10, 4))
ax1 = plt.subplot(1, 2, 1)
ax2 = plt.subplot(1, 2, 2)

demesdraw.tubes(g, ax=ax1)
demesdraw.tubes(g2, ax=ax2, inf_ratio=0.55)

w = ax1.get_xlim()[1] - ax1.get_xlim()[0]
m2 = np.mean(ax2.get_xlim())
ax2.set_xlim((m2 - w/2, m2 + w/2))
ax2.set_ylim(top=ax1.get_ylim()[1])

ax1.set_title("Simulated model")
ax2.set_title("Inferred model")

fig.tight_layout()
fig.savefig("model.misspec.fit.png")
plt.show()
