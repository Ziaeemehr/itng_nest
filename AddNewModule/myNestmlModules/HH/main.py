import nest
import pylab as pl

nest.Install("testModuleInstallmodule")
model = "hh_psc_alpha_nestml"

neuron = nest.Create(model)
parameters = nest.GetDefaults(model)
# for i in parameters:
#     print(i, parameters[i])
# exit(0)

nest.SetStatus(neuron, {'I_e': 630.0})
multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime": True,
                            "record_from": ["V_m"]})
spikedetector = nest.Create("spike_detector",
                            params={"withgid": True,
                                    "withtime": True})
nest.Connect(multimeter, neuron)
nest.Connect(neuron, spikedetector)
nest.Simulate(200.0)
dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["V_m"]
ts = dmm["events"]["times"]
# dSD = nest.GetStatus(spikedetector, keys='events')[0]
# evs = dSD['senders']
# ts = dSD["times"]

pl.plot(ts, Vms)
pl.show()
