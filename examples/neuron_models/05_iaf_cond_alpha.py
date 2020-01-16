import nest
import pylab as pl
import numpy as np
import numpy
import matplotlib.gridspec as gridspec

nest.ResetKernel()
nest.SetKernelStatus(
    {'resolution': 0.01,
     "overwrite_files": True,
     "print_time": False,
     })
N = 2
neuron = nest.Create('iaf_cond_alpha',
                     params={"tau_syn_ex": 1.0,
                             "V_reset": -70.0})

multimeter = nest.Create("multimeter",
                         params={"interval": 0.1,
                                 "record_from":
                                 ["V_m", "g_ex", "g_in"],
                                 "withgid": True,
                                 })

s_ex = nest.Create("spike_generator",
                   params={"spike_times":
                           numpy.array([10.0, 20.0, 50.0])})
s_in = nest.Create("spike_generator",
                   params={"spike_times":
                           numpy.array([15.0, 25.0, 55.0])})

nest.Connect(s_ex, neuron, syn_spec={"weight": 40.0})
nest.Connect(s_in, neuron, syn_spec={"weight": -20.0})
nest.Connect(multimeter, neuron)

nest.Simulate(100.0)  # time in ms

fig = pl.figure(figsize=(6, 5))
gs1 = gridspec.GridSpec(2, 1, hspace=0.0)
axs = []
axs.append(fig.add_subplot(gs1[0]))
axs.append(fig.add_subplot(gs1[1], sharex=axs[0]))

dmm = nest.GetStatus(multimeter)[0]['events']
Vms = dmm["V_m"]
ts = dmm["times"]
axs[0].plot(ts, Vms, label="Voltage")
axs[1].plot(ts, dmm["g_ex"], label="g_ex")
axs[1].plot(ts, dmm["g_in"], label='g_in')

for i in range(2):
    axs[i].legend()
axs[1].set_xlabel("time (ms)")

pl.tight_layout()
pl.savefig("data/05.png", dpi=600)
# pl.show()
