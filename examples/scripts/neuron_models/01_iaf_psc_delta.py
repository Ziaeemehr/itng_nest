import nest
import pylab as pl #to plot figures

nest.ResetKernel()
neuron = nest.Create('iaf_psc_delta')

print 'I_e = ', nest.GetStatus(neuron, "I_e")
print 'V_reset and V_the are ', nest.GetStatus(neuron, ["V_reset", "V_th"])


nest.SetStatus(neuron,{'I_e':375.1})
print 'I_e = ', nest.GetStatus(neuron, "I_e")

multimeter = nest.Create("multimeter")

nest.SetStatus(multimeter,{"withtime":True,
                           "record_from":["V_m"]})

spikedetector = nest.Create("spike_detector",
                            params={"withgid" : True,
                                    "withtime": True})

nest.Connect(multimeter, neuron)
nest.Connect(neuron, spikedetector)

nest.Simulate(1000.0) # time in ms


fig, ax = pl.subplots(2, figsize=(10,8), sharex=True)

dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["V_m"]
ts = dmm["events"]["times"]
ax[0].plot(ts,Vms)

dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD['senders']
ts = dSD["times"]
ax[1].plot(ts,evs,'.')
pl.savefig("data/01.pdf")
# pl.show()