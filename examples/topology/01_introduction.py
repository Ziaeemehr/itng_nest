# import nest
import pylab as pl
import nest.topology as tp

fig = pl.figure(figsize=(6, 5))

layer = tp.CreateLayer({"rows": 5,
                        "columns": 5,
                        "extent": [2.0, 2.0], # use float
                        "center" : [0.0, 1.0],
                        "elements": "iaf_psc_alpha"})
tp.PlotLayer(layer, nodecolor="k", nodesize=50, fig=fig)

pl.xlabel("x")
pl.ylabel("y")
pl.xticks([-1, 0, 1])
pl.yticks([-1, 0, 1])
fig.savefig("figs/01.png")
# pl.show()
