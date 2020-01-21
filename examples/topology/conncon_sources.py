# -*- coding: utf-8 -*-
#
# conncon_sources.py
#
# This file is part of NEST.
#
'''
NEST Topology Module Example

Create two 30x30 layers of iaf_psc_alpha neurons,
connect with convergent projection and rectangular mask,
visualize connection from target perspective.

BCCN Tutorial @ CNS*09
Hans Ekkehard Plesser, UMB
'''

import nest
import nest.topology as topo
import pylab


nest.ResetKernel()
nest.set_verbosity('M_WARNING')

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

topo.ConnectLayers(a, b, {'connection_type': 'convergent',
                          'mask': {'rectangular': {'lower_left': [-0.2, -0.5],
                                                   'upper_right': [0.2, 0.5]}},
                          'kernel': 0.5,
                          'weights': {'uniform': {'min': 0.5, 'max': 2.0}},
                          'delays': 1.0})

fig, ax = pylab.subplots(1)

# plot sources of neurons in different grid locations
for tgt_pos in [[15, 15], [0, 0]]:
    # obtain node id for center
    tgt = topo.GetElement(b, tgt_pos)

    # obtain list of outgoing connections for ctr
    # int() required to cast numpy.int64
    spos = tuple(zip(*[topo.GetPosition([int(conn[0])])[0] for conn in
                       nest.GetConnections(target=tgt)]))

    # scatter-plot
    pylab.scatter(spos[0], spos[1], 20, zorder=10)

    # mark sender position with transparent red circle
    ctrpos = pylab.array(topo.GetPosition(tgt)[0])
    pylab.gca().add_patch(pylab.Circle(ctrpos, radius=0.1, zorder=99,
                                       fc='r', alpha=0.4, ec='none'))

    # mark mask position with open red rectangle
    pylab.gca().add_patch(
        pylab.Rectangle(ctrpos - (0.2, 0.5), 0.4, 1.0, zorder=1,
                        fc='none', ec='r', lw=3))

# mark layer edge
pylab.gca().add_patch(pylab.Rectangle((-1.5, -1.5), 3.0, 3.0, zorder=1,
                                      fc='none', ec='k', lw=3))

# beautify
ax.set_xticks(pylab.arange(-1.5, 1.55, 0.5))
ax.set_yticks(pylab.arange(-1.5, 1.55, 0.5))
ax.grid(True)
ax.axis([-2.0, 2.0, -2.0, 2.0])
ax.set_aspect('equal', 'box')
ax.set_title('Connection sources')

pylab.savefig("figs/conncon.png")
