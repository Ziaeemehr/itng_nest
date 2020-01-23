# -*- coding: utf-8 -*-
#
# conncon_targets.py
#
# This file is part of NEST.
'''
NEST Topology Module Example

Create two 30x30 layers of iaf_psc_alpha neurons,
connect with convergent projection and rectangular mask,
visualize connections from source perspective.
'''

import pylab
import nest
import nest.topology as topo

nest.ResetKernel()

# create two test layers
a = topo.CreateLayer({'columns': 30,
                      'rows': 30,
                      'extent': [3.0, 3.0],
                      'elements': 'iaf_psc_alpha',
                      'edge_wrap': True})
b = topo.CreateLayer({'columns': 30,
                      'rows': 30,
                      'extent': [3.0, 3.0],
                      'elements': 'iaf_psc_alpha',
                      'edge_wrap': True})

conndict = {'connection_type': 'convergent',
            'mask': {'rectangular': {'lower_left': [-0.2, -0.5],
                                     'upper_right': [0.2, 0.5]}},
            'kernel': 0.5,
            'weights': {'uniform': {'min': 0.5, 'max': 2.0}},
            'delays': 1.0}
topo.ConnectLayers(a, b, conndict)

fig, ax = pylab.subplots(1)

# plot targets of two source neurons into same figure, with mask
for src_pos in [[15, 15], [0, 0]]:
    # obtain node id for center
    src = topo.GetElement(a, src_pos)
    topo.PlotTargets(src, b, mask=conndict['mask'], fig=fig)

# beautify
ax.set_xticks(pylab.arange(-1.5, 1.55, 0.5))
ax.set_yticks(pylab.arange(-1.5, 1.55, 0.5))
ax.grid(True)
ax.axis([-2.0, 2.0, -2.0, 2.0])
ax.set_aspect('equal', 'box')
ax.set_title('Connection targets')

pylab.savefig('figs/conncon_targets.png')
