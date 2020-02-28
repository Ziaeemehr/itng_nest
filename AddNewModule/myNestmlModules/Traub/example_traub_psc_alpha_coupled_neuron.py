import nest
import pylab as pl


# nest.Install("testModuleInstallmodule")
nest.Install("nestml_module")
model = "traub_psc_alpha_nestml"

nest.set_verbosity("M_WARNING")
nest.ResetKernel()


dt = 0.01
tfinal = 300.0

nest.SetKernelStatus({
    "resolution": dt})

neuron1 = nest.Create(model, 5)
nest.SetStatus(neuron1, {'I_e': 150.0})

neuron2 = nest.Create(model)

multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime": True,
                            "record_from": ["V_m"],
                            "interval": dt})

nest.Connect(neuron1, neuron2, syn_spec={"weight": 20,
                                         "delay": 1.0})
nest.Connect(multimeter, neuron2)

spikedetector = nest.Create("spike_detector",
                            params={"withgid": True,
                                    "withtime": True})
nest.Connect(neuron1, spikedetector)

nest.Simulate(tfinal)
dmm = nest.GetStatus(multimeter)[0]
Voltages = dmm["events"]["V_m"]
tv = dmm["events"]["times"]
dSD = nest.GetStatus(spikedetector, keys='events')[0]
Spikes = dSD['senders']
ts = dSD["times"]

fig, ax = pl.subplots(2, figsize=(8, 6), sharex=True)
ax[0].plot(tv, Voltages, lw=2, color="k")
ax[1].plot(ts, Spikes, 'k.')
ax[1].set_xlabel("Time [ms]")
ax[1].set_xlim(0, tfinal)
ax[0].set_title("recording from PSP", fontsize=14)
ax[1].set_ylabel("Spikes")
ax[0].set_ylabel("v [ms]")
pl.savefig("coupled.png")
