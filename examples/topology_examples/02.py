# import nest
import pylab as pl
import numpy as np
import nest.topology as tp

fig, ax = pl.subplots(1, figsize=(6, 5))

nc, nr = 5, 3
d = 0.1
layer = tp.CreateLayer({"columns": nc,
                        "rows": nr,
                        "elements": "iaf_psc_alpha",
                        "extent": [nc * d, nr * d],
                        "center": [nc* d / 2., 0.]})

# tp.PlotLayer([layer1[0], layer2[0]], nodecolor="k", nodesize=50, fig=fig)
tp.PlotLayer(layer, nodecolor="k", nodesize=50, fig=fig)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect('equal', 'box')
ax.set_xticks(np.arange(0,0.6,0.1))
fig.savefig("figs/02.png")
# pl.show()
