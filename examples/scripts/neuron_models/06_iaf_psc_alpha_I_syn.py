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
neurons = nest.Create('iaf_psc_alpha', N)

# print 'I_e = ', nest.GetStatus(neurons, "I_e")
# print 'V_reset and V_the are ', nest.GetStatus(neurons, ["V_reset", "V_th"])
Vrest = nest.GetStatus([neurons[0]], "V_reset")[0]
Vth = nest.GetStatus([neurons[0]], "V_th")[0]

Vms = Vrest+(Vth-Vrest)*np.random.rand(len(neurons))
nest.SetStatus(neurons, "V_m", Vms)

nest.SetStatus(neurons, [{'I_e': 376.0}, {'I_e': 0.0}])
print 'I_e = ', nest.GetStatus(neurons, "I_e")

multimeters = nest.Create("multimeter", N)

nest.SetStatus(multimeters, {"withtime": True,
                             "record_from": ["V_m",
                                             "I_syn_ex",
                                             "I_syn_in"]})

nest.Connect([neurons[0]], [neurons[1]],
             syn_spec={"weight": -20.0, "delay": 1.0})

nest.Connect(multimeters, neurons, 'one_to_one')
nest.Simulate(200.0)  # time in ms


fig = pl.figure(figsize=(6, 5))
gs1 = gridspec.GridSpec(2, 1, hspace=0.0)
axs = []
axs.append(fig.add_subplot(gs1[0]))
axs.append(fig.add_subplot(gs1[1], sharex=axs[0]))

color = ['k', 'r']
for i in range(N):
    dmm = nest.GetStatus(multimeters)[i]['events']
    Vms = dmm["V_m"]
    ts = dmm["times"]
    axs[0].plot(ts, Vms, lw=2, label=str(i+1), c=color[i])

for i in range(N):
    dmm = nest.GetStatus(multimeters)[1]['events']
    ts = dmm["times"]
    axs[1].plot(ts, dmm["I_syn_ex"], label='ex-'+str(i+1))
    axs[1].plot(ts, dmm["I_syn_in"], label='in-'+str(i+1))

# axs[0].set_title('two connected neurons')
axs[0].legend(loc='upper right')
axs[1].legend(loc='upper right')
axs[1].set_xlabel("Time (ms)")
axs[0].set_ylabel("Voltage(mV)")
axs[1].set_ylabel(r"$I_{syn_{ex}}$(pA)")
pl.tight_layout()
pl.savefig("data/06.png", dpi=300)
# pl.show()
