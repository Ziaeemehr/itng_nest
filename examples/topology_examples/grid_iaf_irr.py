# -*- coding: utf-8 -*-
#
# grid_iaf_irr.py

'''
NEST Topology Module Example

Create layer of 12 freely placed iaf_psc_alpha neurons, visualize

BCCN Tutorial @ CNS*09
Hans Ekkehard Plesser, UMB
'''

import nest
import pylab
import random
import nest.topology as topo

nest.ResetKernel()

fig, ax = pylab.subplots(1)

# generate list of 12 (x,y) pairs
pos = [[random.uniform(-0.75, 0.75), random.uniform(-0.5, 0.5)]
       for j in range(12)]

l1 = topo.CreateLayer({'extent': [2., 1.5],
                       'positions': pos,
                       'elements': 'iaf_psc_alpha'})

nest.PrintNetwork()
nest.PrintNetwork(2)
nest.PrintNetwork(2, l1)

topo.PlotLayer(l1, nodesize=50, fig=pylab.gcf())

# beautify

ax.axis([-1.0, 1.0, -0.75, 0.75])
ax.set_aspect('equal', 'box')
ax.set_xticks((-0.75, -0.25, 0.25, 0.75))
ax.set_yticks((-0.5, 0, 0.5))
ax.grid(True)
ax.set_xlabel('4 Columns, Extent: 1.5')
ax.set_ylabel('2 Rows, Extent: 1.0')

pylab.savefig('figs/grid_iaf_irr.png')
