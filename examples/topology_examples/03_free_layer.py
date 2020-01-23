# import nest
import pylab as pl
import numpy as np
import nest.topology as tp

fig, ax = pl.subplots(1, figsize=(6, 5))

pos = [[np.random.uniform(-0.5, 0.5),
        np.random.uniform(-0.5, 0.5)]
       for i in range(500)]

layer = tp.CreateLayer({"positions": pos,
                        "elements": "iaf_psc_alpha"})

# tp.PlotLayer([layer1[0], layer2[0]], nodecolor="k", nodesize=50, fig=fig)
tp.PlotLayer(layer, nodecolor="k", nodesize=50, fig=fig)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect('equal', 'box')
fig.savefig("figs/03.png")
# pl.show()
