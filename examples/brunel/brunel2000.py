# brunel_delta_nest.py


import nest 
import nest.raster_plot
import pylab as pl 
from time import time 
import numpy 

nest.ResetKernel()
n = 4  # number of threads
nest.SetKernelStatus({'local_num_threads': n})


startbuild = time()

# simulation parameters
dt = 0.1
simtime = 1000.0
delay = 1.5  # synaptic delay in ms

# Definition of the parameters crucial for asynchronous irregular firing
# of the neurons.

g = 5.0  # ratio inhibitory weight/excitatory weight
eta = 2.0  # external rate relative to threshold rate
epsilon = 0.1

# Definition of the number of neurons in the network and the number of
# neuron recorded from

order = 2500
NE = 4 * order  # number of excitatory neurons
NI = 1 * order  # number of inhibitory neurons
N_neurons = NE + NI  # number of neurons in total
N_rec = 50  # record from 50 neurons


# Definition of connectivity parameter

CE = int(epsilon * NE)  # number of excitatory synapses per neuron
CI = int(epsilon * NI)  # number of inhibitory synapses per neuron
C_tot = int(CI + CE)  # total number of synapses per neuron



# Initialization of the parameters of the integrate and fire neuron and
# the synapses. The parameter of the neuron are stored in a dictionary.

tauMem = 20.0  # time constant of membrane potential in ms
theta = 20.0  # membrane threshold potential in mV
neuron_params = {"C_m": 1.0,
                 "tau_m": tauMem,
                 "t_ref": 2.0,
                 "E_L": 0.0,
                 "V_reset": 0.0,
                 "V_m": 0.0,
                 "V_th": theta}
J = 0.1  # postsynaptic amplitude in mV
J_ex = J  # amplitude of excitatory postsynaptic potential
J_in = -g * J_ex  # amplitude of inhibitory postsynaptic potential


'''
Definition of threshold rate, which is the external rate needed to fix
the membrane potential around its threshold, the external firing rate
and the rate of the poisson generator which is multiplied by the
in-degree CE and converted to Hz by multiplication by 1000.
'''

nu_th = theta / (J * CE * tauMem)
nu_ex = eta * nu_th
p_rate = 1000.0 * nu_ex * CE


nest.SetKernelStatus({"resolution": dt, "print_time": True,
                      "overwrite_files": True})

print "Building network"


# Creating neurons and devices
# ------------------------------------------------------- #
nest.SetDefaults("iaf_psc_delta", neuron_params)
nest.SetDefaults("poisson_generator", {"rate": p_rate})

nodes_ex = nest.Create("iaf_psc_delta", NE)
nodes_in = nest.Create("iaf_psc_delta", NI)
noise = nest.Create("poisson_generator")
espikes = nest.Create("spike_detector")
ispikes = nest.Create("spike_detector")

nest.SetStatus(espikes, [{"label": "data/brunel-py-ex",
                          "withtime": True,
                          "withgid": True,
                          "to_file": True}])

nest.SetStatus(ispikes, [{"label": "data/brunel-py-in",
                          "withtime": True,
                          "withgid": True,
                          "to_file": True}])


print "Connecting devices"


# Connecting the network
# ------------------------------------------------------- #
nest.CopyModel("static_synapse", "excitatory",
               {"weight": J_ex, "delay": delay})
nest.CopyModel("static_synapse", "inhibitory",
               {"weight": J_in, "delay": delay})


nest.Connect(noise, nodes_ex, syn_spec="excitatory")
nest.Connect(noise, nodes_in, syn_spec="excitatory")

nest.Connect(nodes_ex[:N_rec], espikes, syn_spec="excitatory")
nest.Connect(nodes_in[:N_rec], ispikes, syn_spec="excitatory")

print "Connecting network"

print "Excitatory connections"

conn_params_ex = {'rule': 'fixed_indegree', 'indegree': CE}
nest.Connect(nodes_ex, nodes_ex + nodes_in, conn_params_ex, "excitatory")

print "Inhibitory connections"
conn_params_in = {'rule': 'fixed_indegree', 'indegree': CI}
nest.Connect(nodes_in, nodes_ex + nodes_in, conn_params_in, "inhibitory")

endbuild = time()

print "Simulating"

nest.Simulate(simtime)

endsimulate = time()


events_ex = nest.GetStatus(espikes, "n_events")[0]
events_in = nest.GetStatus(ispikes, "n_events")[0]


'''
Calculation of the average firing rate of the excitatory and the
inhibitory neurons by dividing the total number of recorded spikes by
the number of neurons recorded from and the simulation time. The
multiplication by 1000.0 converts the unit 1/ms to 1/s=Hz.
'''

rate_ex = events_ex / simtime * 1000.0 / N_rec
rate_in = events_in / simtime * 1000.0 / N_rec

'''
Reading out the number of connections established using the excitatory
and inhibitory synapse model. The numbers are summed up resulting in
the total number of synapses.
'''

num_synapses = (nest.GetDefaults("excitatory")["num_connections"] +
                nest.GetDefaults("inhibitory")["num_connections"])

'''
Establishing the time it took to build and simulate the network by
taking the difference of the pre-defined time variables.
'''

build_time = endbuild - startbuild
sim_time = endsimulate - endbuild

'''
Printing the network properties, firing rates and building times.
'''

print("Brunel network simulation (Python)")
print("Number of neurons : {0}".format(N_neurons))
print("Number of synapses: {0}".format(num_synapses))
print("       Exitatory  : {0}".format(int(CE * N_neurons) + N_neurons))
print("       Inhibitory : {0}".format(int(CI * N_neurons)))
print("Excitatory rate   : %.2f Hz" % rate_ex)
print("Inhibitory rate   : %.2f Hz" % rate_in)
print("Building time     : %.2f s" % build_time)
print("Simulation time   : %.2f s" % sim_time)

'''
Plot a raster of the excitatory neurons and a histogram.
'''

nest.raster_plot.from_device(espikes, hist=True)
pl.show()
