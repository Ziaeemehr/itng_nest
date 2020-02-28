import nest
import pylab as pl
import matplotlib

matplotlib.rc('font', size=14)


# nest.Install("testModuleInstallmodule")
nest.Install("nestml_module")
model = "traub_receptors_nestml"

nest.set_verbosity("M_WARNING")
nest.ResetKernel()


dt = 0.01
tfinal = 1000.0


nest.SetKernelStatus({
    "resolution": dt})

neuron1 = nest.Create(model, 1)
nest.SetStatus(neuron1, {'I_e': 100.0})

neuron2 = nest.Create(model)
nest.SetStatus(neuron2, {"AMPA_Tau_1": 0.1,
                         "AMPA_Tau_2": 2.4,
                         "AMPA_g_peak": 0.1,
                         })


multimeter = nest.Create("multimeter", 2)
nest.SetStatus([multimeter[0]], {"withtime": True,
                                 "record_from": ["V_m"],
                                 "interval": dt})


record_from = ["V_m", "I_syn_ampa", "I_syn_nmda", "I_syn_gaba_a", "I_syn_gaba_b"]

nest.SetStatus([multimeter[1]], {"withtime": True,
                                 "record_from": record_from,
                                 "interval": dt})
#! {'AMPA': 1, 'NMDA': 2, 'GABA_A': 3, 'GABA_B': 4}
# nest.Connect(neuron1, neuron2)
nest.Connect(neuron1, neuron2, syn_spec={"receptor_type": 1})  # AMPA
nest.Connect(neuron1, neuron2, syn_spec={"receptor_type": 2})  # NMDA
nest.Connect(neuron1, neuron2, syn_spec={"receptor_type": 3})  # GABAA
nest.Connect(neuron1, neuron2, syn_spec={"receptor_type": 4})  # GABAB

nest.Connect([multimeter[0]], neuron1, "one_to_one")
nest.Connect([multimeter[1]], neuron2)

spikedetector = nest.Create("spike_detector",
                            params={"withgid": True,
                                    "withtime": True})
nest.Connect(neuron1, spikedetector)
nest.Simulate(tfinal)

fig, ax = pl.subplots(3, figsize=(8, 6), sharex=True)


dmm = nest.GetStatus(multimeter)[1]
Voltages = dmm["events"]["V_m"]
tv = dmm["events"]["times"]

ax[0].plot(tv, Voltages, lw=2, label=str(2))

labels = ["ampa", "nmda", "gaba_a", "gaba_b"]
j = 0
for i in record_from[1:]:
    g = dmm["events"][i]
    ax[1].plot(tv, g, lw=2, label=labels[j])
    j +=1

dSD = nest.GetStatus(spikedetector, keys='events')[0]
Spikes = dSD['senders']
ts = dSD["times"]

ax[2].plot(ts, Spikes, 'k.')
ax[2].set_xlabel("Time [ms]")
ax[2].set_xlim(0, tfinal)
ax[2].set_ylabel("Spikes")

ax[0].set_title("recording from PSP")
ax[0].set_ylabel("v [ms]")
# ax[0].set_ylim(-68, -66)
ax[1].set_ylabel("I_syn")
ax[1].legend(frameon=False, loc="upper right")

pl.savefig("figs/ampa_coupled.png")
# pl.show()