import nest
import pylab as pl
import numpy as np
import matplotlib.gridspec as gridspec

nest.ResetKernel()
nest.SetKernelStatus(
    {'resolution': 0.01,
     "overwrite_files": True,
     "print_time": False,
     })
N = 2
neurons = nest.Create('iaf_cond_alpha', N)

# print 'I_e = ', nest.GetStatus(neurons, "I_e")
# print 'V_reset and V_the are ', nest.GetStatus(neurons, ["V_reset", "V_th"])
Vrest = nest.GetStatus([neurons[0]], "V_reset")[0]
Vth = nest.GetStatus([neurons[0]], "V_th")[0]

Vms = Vrest+(Vth-Vrest)*np.random.rand(len(neurons))
nest.SetStatus(neurons, "V_m", Vms)

# print nest.GetStatus(neurons, "V_m")


nest.SetStatus(neurons, [{'I_e': 300.0}, {'I_e': 0.0}])
print 'I_e = ', nest.GetStatus(neurons, "I_e")

multimeters = nest.Create("multimeter", N)

nest.SetStatus(multimeters, {"withtime": True,
                             "record_from": ["V_m",
                                             "g_ex",
                                             "g_in"]})
spikedetector = nest.Create("spike_detector",
                            params={"withgid": True,
                                    "withtime": True})

nest.Connect([neurons[0]], [neurons[1]], syn_spec={"weight": 20, "delay": 1.0})

# for i in range(2):
#     conn = nest.GetConnections([neurons[i]])
#     print nest.GetStatus(conn, ['source', 'target'])

nest.Connect(multimeters, neurons, 'one_to_one')
nest.Connect(neurons, spikedetector)

# print nest.GetDefaults('iaf_cond_alpha')['recordables']
params = nest.GetDefaults('iaf_psc_alpha')
# print params
# for i, j in zip(params.keys(), params.values()):
#     print i, j

# print("iaf_cond_alpha recordables: {0}".format(
#   nest.GetDefaults("iaf_cond_alpha")["recordables"]))

nest.Simulate(200.0)  # time in ms


fig = pl.figure(figsize=(6, 5))
gs1 = gridspec.GridSpec(4, 1, hspace=0.0)
axs = []
axs.append(fig.add_subplot(gs1[:3]))
axs.append(fig.add_subplot(gs1[3], sharex=axs[0]))

# fig, ax = pl.subplots(2, figsize=(6, 5), sharex=True)
for i in range(N):
    dmm = nest.GetStatus(multimeters)[i]
    Vms = dmm["events"]["V_m"]
    gex = dmm["events"][""]
    ts = dmm["events"]["times"]
    axs[0].plot(ts, Vms)

dSD = nest.GetStatus(spikedetector, keys='events')[0]
evs = dSD['senders']
ts = dSD["times"]
axs[1].plot(ts, evs, '.')
pl.tight_layout()
pl.savefig("data/04.pdf")
# pl.show()
