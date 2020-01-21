# connex.py
#
'''
NEST Topology Module Example

Create two 30x30 layers of iaf_psc_alpha neurons,
connect with circular mask, flat probability,
visualize.

Divergent connection When creating a divergent connection, Topology 
visits each node in the source layer and selects target nodes from 
the target layer. Masks, kernels, and boundary conditions are applied 
in the target layer.
'''

import nest
import nest.topology as topo
import pylab

nest.ResetKernel()

# create two test layers
a = topo.CreateLayer({'columns': 30, 'rows': 30,
                      'extent': [3.0, 3.0],
                      'elements': 'iaf_psc_alpha',
                      "edge_wrap": False})
b = topo.CreateLayer({'columns': 30, 'rows': 30,
                      'extent': [3.0, 3.0],
                      'elements': 'iaf_psc_alpha',
                      'edge_wrap': False})
conndict = {'connection_type': 'divergent',
            'mask': {'circular': {'radius': 0.5}},
            'kernel': 0.5,
            'weights': {'uniform': {'min': 0.5, 'max': 2.0}},
            'delays': 1.0}
topo.ConnectLayers(a, b, conndict)

# plot targets of neurons in different grid locations
fig, ax = pylab.subplots(1)

# plot targets of two source neurons into same figure, with mask
for src_pos in [[15, 10], [0, 0]]:
    # obtain node id for center
    src = topo.GetElement(a, src_pos)
    topo.PlotTargets(src, b, mask=conndict['mask'], fig=pylab.gcf())

# beautify
ax.set_xticks(pylab.arange(-1.5, 1.55, 0.5))
ax.set_yticks(pylab.arange(-1.5, 1.55, 0.5))
ax.grid(True)
ax.axis([-2.0, 2.0, -2.0, 2.0])
ax.set_aspect('equal', 'box')
ax.set_title('Connection targets')

pylab.savefig('figs/connex.png')
