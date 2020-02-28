import nest
import pylab as pl



# nest.Install("testModuleInstallmodule")
nest.Install("nestml_module")
model = "traub_psc_alpha_nestml"

dt = 0.01
tfinal = 100.0
nest.SetKernelStatus({
    "resolution": dt})

neuron = nest.Create(model)
parameters = nest.GetDefaults(model)
# for i in parameters:
#     print(i, parameters[i])
# exit(0)

nest.SetStatus(neuron, {'I_e': 150.0})
multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime": True,
                            "record_from": ["V_m"],
                            "interval": dt})
spikedetector = nest.Create("spike_detector",
                            params={"withgid": True,
                                    "withtime": True})
nest.Connect(multimeter, neuron)
nest.Connect(neuron, spikedetector)
nest.Simulate(tfinal)
dmm = nest.GetStatus(multimeter)[0]
Voltages = dmm["events"]["V_m"]
tv = dmm["events"]["times"]
dSD = nest.GetStatus(spikedetector, keys='events')[0]
Spikes = dSD['senders']
ts = dSD["times"]

fig, ax = pl.subplots(2, figsize=(8, 6), sharex=True)
ax[0].plot(tv, Voltages, lw=2, color="k")
ax[1].plot(ts, Spikes, 'ko')
ax[1].set_xlabel("Time [ms]")
ax[1].set_xlim(0, tfinal)
ax[1].set_ylabel("Spikes")
ax[0].set_ylabel("v [ms]")
pl.savefig("fig.png")
