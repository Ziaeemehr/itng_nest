# -*- coding: utf-8 -*-
#
# gaussex.py
#
# This file is part of NEST.
'''
NEST Topology Module Example

Create two layers of 30x30 elements and connect
them using a Gaussian probabilistic kernel, visualize.
'''

import pylab
import nest
import nest.topology as topo


nest.ResetKernel()

# create two test layers
a = topo.CreateLayer({'columns': 30, 'rows': 30, 'extent': [3.0, 3.0],
                      'elements': 'iaf_psc_alpha'})
b = topo.CreateLayer({'columns': 30, 'rows': 30, 'extent': [3.0, 3.0],
                      'elements': 'iaf_psc_alpha'})

conndict = {'connection_type': 'divergent',
            'mask': {'circular': {'radius': 3.0}},
            'kernel': {'gaussian': {'p_center': 1.0, 'sigma': 0.5}},
            'weights': 1.0,
            'delays': 1.0}
topo.ConnectLayers(a, b, conndict)

# plot targets of neurons in different grid locations
fig, ax = pylab.subplots(1)

# plot targets of two source neurons into same figure, with mask
# use different colors
for src_pos, color in [([15, 15], 'blue'), ([0, 0], 'green')]:

    # obtain node id for center
    src = topo.GetElement(a, src_pos)
    topo.PlotTargets(src, b,
                     mask=conndict['mask'],
                     kernel=conndict['kernel'],
                     src_color=color,
                     tgt_color=color,
                     mask_color=color,
                     kernel_color=color,
                     src_size=100,
                     fig=fig)

# beautify
ax.set_xticks(pylab.arange(-1.5, 1.55, 0.5))
ax.set_yticks(pylab.arange(-1.5, 1.55, 0.5))
ax.grid(True)
ax.axis([-2.0, 2.0, -2.0, 2.0])
ax.set_aspect('equal', 'box')
ax.set_title('Connection targets, Gaussian kernel')

pylab.savefig('figs/gaussex.png')
