# https://www.nest-simulator.org/py_sample/lin_rate_ipn_network/
import nest
import pylab
import numpy
import pylab

dt = 0.1  # the resolution in ms
T = 100.0  # Simulation time in ms

order = 50
NE = int(4 * order)  # number of excitatory neurons
NI = int(1 * order)  # number of inhibitory neurons
N = int(NE+NI)       # total number of neurons

# Definition of the connections
d_e = 5.   # delay of excitatory connections in ms
g = 5.0  # ratio inhibitory weight/excitatory weight
epsilon = 0.1  # connection probability
w = 0.1/numpy.sqrt(N)  # excitatory connection strength

KE = int(epsilon * NE)  # number of excitatory synapses per neuron (outdegree)
KI = int(epsilon * NI)  # number of inhibitory synapses per neuron (outdegree)
K_tot = int(KI + KE)  # total number of synapses per neuron
connection_rule = 'fixed_outdegree'  # connection rule

# Definition of the neuron model and its neuron parameters
neuron_model = 'lin_rate_ipn'  # neuron model
neuron_params = {'linear_summation': True,
                 # type of non-linearity (not affecting linear rate models)
                 'tau': 10.0,
                 # time constant of neuronal dynamics in ms
                 'mean': 2.0,
                 # mean of Gaussian white noise input
                 'std': 5.
                 # standard deviation of Gaussian white noise input
                 }

nest.ResetKernel()
nest.SetKernelStatus({"resolution": dt, "use_wfr": False,
                      "print_time": True,
                      "overwrite_files": True})

print("Building network")  

nest.SetDefaults(neuron_model, neuron_params)

n_e = nest.Create(neuron_model, NE)
n_i = nest.Create(neuron_model, NI)

mm = nest.Create('multimeter', params={'record_from': ['rate'],
                                       'interval': dt})

syn_e = {'weight': w, 'delay': d_e, 'model': 'rate_connection_delayed'}
syn_i = {'weight': -g*w, 'model': 'rate_connection_instantaneous'}
conn_e = {'rule': connection_rule, 'outdegree': KE}
conn_i = {'rule': connection_rule, 'outdegree': KI}

# Connect rate units
nest.Connect(n_e, n_e, conn_e, syn_e)
nest.Connect(n_i, n_i, conn_i, syn_i)
nest.Connect(n_e, n_i, conn_i, syn_e)
nest.Connect(n_i, n_e, conn_e, syn_i)

# Connect recording device to rate units
nest.Connect(mm, n_e+n_i)

nest.Simulate(T)

data = nest.GetStatus(mm)[0]['events']
rate_ex = data['rate'][numpy.where(data['senders'] == n_e[0])]
rate_in = data['rate'][numpy.where(data['senders'] == n_i[0])]
times = data['times'][numpy.where(data['senders'] == n_e[0])]

pylab.figure()
pylab.plot(times, rate_ex, label='excitatory')
pylab.plot(times, rate_in, label='inhibitory')
pylab.legend()
pylab.xlabel('time (ms)')
pylab.ylabel('rate (a.u.)')
pylab.savefig("f")
pylab.show()