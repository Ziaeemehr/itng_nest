import nest
import pylab as pl

nest.Install("testModuleInstallmodule")
model = "izhikevich_nestml"

neuron = nest.Create(model)
parameters = nest.GetDefaults(model)
# for i in parameters:
#     print(i, parameters[i])
# exit(0)

# nest.SetStatus(neuron, {'I_ext': 630.0})
multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime": True,
                            "record_from": ["v"]})
spikedetector = nest.Create("spike_detector",
                            params={"withgid": True,
                                    "withtime": True})
nest.Connect(multimeter, neuron)
nest.Connect(neuron, spikedetector)
nest.Simulate(200.0)
dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["v"]
ts = dmm["events"]["times"]
# dSD = nest.GetStatus(spikedetector, keys='events')[0]
# evs = dSD['senders']
# ts = dSD["times"]

pl.plot(ts, Vms)
pl.show()
