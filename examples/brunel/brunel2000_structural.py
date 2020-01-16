# brunel_delta_nest.py

# Abolfazl Ziaeemehr
# Institute for Advanced Studies in
# Basic Sciences (IASBS)
# tel: +98 3315 2148
# github.com/ziaeemehr
# a.ziaeemehr@gmail.com


import nest 
import nest.raster_plot
import pylab as pl 
from time import time 
import numpy 
import os




class Brunel2000(object):
    '''
    Implementation of the sparsely connected random network,
    described by Brunel (2000) J. Comp. Neurosci.
    Parameters are chosen for the asynchronous irregular
    state (AI).
    '''

    n_threads = 4  # number of threads

    dt = 0.1
    simtime = 1000.0
    delay = 1.5  # synaptic delay in ms

    g = 5.0  # ratio inhibitory weight/excitatory weight
    eta = 2.0  # external rate relative to threshold rate
    epsilon = 0.1

    order = 2500
    NE = 4 * order  # number of excitatory neurons
    NI = 1 * order  # number of inhibitory neurons
    N_neurons = NE + NI  # number of neurons in total
    N_rec = 50  # record from 50 neurons
    
    tau_m = 20.0  # time constant of membrane potential in ms
    theta = 20.0  # membrane threshold potential in mV
    neuron_params = {"C_m": 1.0,
                    "tau_m": tau_m,
                    "t_ref": 2.0,
                    "E_L": 0.0,
                    "V_reset": 0.0,
                    "V_m": 0.0,
                    "V_th": theta}
    J = 0.1  # postsynaptic amplitude in mV
    J_ex = J  # amplitude of excitatory postsynaptic potential
    J_in = -g * J_ex  # amplitude of inhibitory postsynaptic potential

    data_path = "data/"

    built = False       # True, if build() was called
    connected = False   # True, if connect() was called

    def __init__(self):
        self.name = self.__class__.__name__
        nest.ResetKernel()

        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)

        nest.SetKernelStatus({
                      "resolution": self.dt, 
                      "print_time": True,
                      "overwrite_files": True, 
                      "data_path": self.data_path,
                      "local_num_threads": self.n_threads})
        
        # Create and seed RNGs
        msd = 1000
        # master seed
        n_vp = nest.GetKernelStatus('total_num_virtual_procs')
        msdrange1 = range(msd, msd + n_vp)
        self.pyrngs = [numpy.random.RandomState(s) for s in msdrange1]
        msdrange2 = range(msd + n_vp + 1, msd + 1 + 2 * n_vp)
        nest.SetKernelStatus({'grng_seed': msd + n_vp,
                            'rng_seeds': msdrange2})        
        
    def calibrate(self):
        '''
        Compute all parameter dependent variables of the model.
        '''
        self.CE = int(self.epsilon * self.NE)  # number of excitatory synapses per neuron
        self.CI = int(self.epsilon * self.NI)  # number of inhibitory synapses per neuron
        C_tot = int(self.CI + self.CE)  # total number of synapses per neuron
        nu_th = self.theta / (self.J * self.CE * self.tau_m)
        nu_ex = self.eta * nu_th
        p_rate = 1000.0 * nu_ex * self.CE


        nest.SetDefaults("iaf_psc_delta", self.neuron_params)
        nest.SetDefaults("poisson_generator", {"rate": p_rate})


    def build(self):
        '''
        Create all nodes, used in the model.
        '''
        if self.built:
            return
        self.calibrate()

        self.nodes_ex = nest.Create("iaf_psc_delta", self.NE)


        self.nodes_in = nest.Create("iaf_psc_delta", self.NI)
        self.noise = nest.Create("poisson_generator")
        self.espikes = nest.Create("spike_detector")
        self.ispikes = nest.Create("spike_detector")

        nest.SetStatus(self.espikes, [{"label": "brunel-py-ex",
                                "withtime": True,
                                "withgid": True,
                                "to_file": True}])

        nest.SetStatus(self.ispikes, [{"label": "brunel-py-in",
                                "withtime": True,
                                "withgid": True,
                                "to_file": True}])
        node_info = nest.GetStatus(self.nodes_ex+self.nodes_in)
        local_nodes = [(ni['global_id'], ni['vp'])
                    for ni in node_info if ni['local']]
        for gid, vp in local_nodes:
            nest.SetStatus([gid], {'V_m': self.pyrngs[vp].uniform(-self.theta, self.theta)})
    
        self.built = True
    
    def connect(self):
        '''
        Connect all nodes in the model.
        '''

        if self.connected:
            return 
        if not self.built:
            self.build()

        nest.CopyModel("static_synapse", "excitatory")

        nest.CopyModel("static_synapse", "inhibitory",
                    {"weight": self.J_in, "delay": self.delay})

        nest.CopyModel("static_synapse", "excitatory_input",
                    {"weight": self.J_ex, "delay": self.delay})

        nest.Connect(self.noise, self.nodes_ex, syn_spec="excitatory_input")
        nest.Connect(self.noise, self.nodes_in, syn_spec="excitatory_input")

        nest.Connect(self.nodes_ex[:self.N_rec], self.espikes)
        nest.Connect(self.nodes_in[:self.N_rec], self.ispikes)

        print "Excitatory connections"

        conn_params_ex = {'rule': 'fixed_indegree', 'indegree': self.CE}

        nest.Connect(self.nodes_ex, self.nodes_ex + self.nodes_in, conn_params_ex,
                    {'model': 'excitatory',
                    'delay': self.delay,
                    'weight': {'distribution': 'uniform',
                                'low': 0.5*self.J_ex,
                                'high': 1.5*self.J_ex}})

        print "Inhibitory connections"
        conn_params_in = {'rule': 'fixed_indegree', 'indegree': self.CI}
        nest.Connect(self.nodes_in, self.nodes_ex + self.nodes_in, conn_params_in, "inhibitory")

        self.connected = True

    def run(self, simtime=300):
        '''
        Simulate the model for simtime milliseconds and print the firing 
        rates of the network during this period.
        '''
        if not self.connected:
            self.connect()
        nest.Simulate(simtime)

        events_ex = nest.GetStatus(self.espikes, "n_events")[0]
        events_in = nest.GetStatus(self.ispikes, "n_events")[0]

        rate_ex = events_ex / simtime * 1000.0 / self.N_rec
        rate_in = events_in / simtime * 1000.0 / self.N_rec

        num_synapses = (nest.GetDefaults("excitatory")["num_connections"] +
                        nest.GetDefaults("inhibitory")["num_connections"])

        print("Brunel network simulation (Python)")
        print("Number of neurons : {0}".format(self.N_neurons))
        print("Number of synapses: {0}".format(num_synapses))
        print("       Exitatory  : {0}".format(int(self.CE * self.N_neurons) + self.N_neurons))
        print("       Inhibitory : {0}".format(int(self.CI * self.N_neurons)))
        print("Excitatory rate   : %.2f Hz" % rate_ex)
        print("Inhibitory rate   : %.2f Hz" % rate_in)
        # print("Building time     : %.2f s" % build_time)
        # print("Simulation time   : %.2f s" % sim_time)

        nest.raster_plot.from_device(self.espikes, hist=True)
        pl.show()


net = Brunel2000()
net.run(1000.0)
