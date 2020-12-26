import numpy as np
import pylab as pl
import nest


nest.SetKernelStatus(
            {'resolution': self.dt,
             "data_path": self.data_path,
             "overwrite_files": True,
             "print_time": False,
             })
nest.set_verbosity("M_WARNING")

# Create neuron
neuron = nest.Create("iaf_psc_alpha")
nest.GetStatus(neuron, "I_e")
nest.GetStatus(neuron, ["V_reset", "V_th"])
nest.SetStatus(neuron, {"I_e": 376.0})

ndict = {"I_e": 200.0, "tau_m": 20.0}
neuronpop = nest.Create("iaf_psc_alpha", 100, params=ndict)
# params is a dict of list of dict


ndict = {"I_e": 200.0, "tau_m": 20.0}
# nest.SetDefaults(model, params)
nest.SetDefaults("iaf_psc_alpha", ndict)

print (nest.GetStatus(neuronpop[1:10], "I_e"))
print (nest.GetDefaults('iaf_psc_alpha')['tau_m'])


# If batches of neurons should be of the same model but using
# different parameters
edict = {"I_e": 200.0, "tau_m": 20.0}
# CopyModel(existing, new, params=None)
nest.CopyModel("iaf_psc_alpha", "exc_iaf_neuron")
nest.SetDefaults("exc_iaf_neuron", edict)
# # or in one step ----------------------------------- #
idict = {"I_e": 300.0}
nest.CopyModel("iaf_psc_alpha", "inh_iaf_neuron", params=idict)


# populations with an inhomogeneous set of parameters
# supply a list of dictionaries of the same length as
# the number of neurons (or synapses) created
parameter_list = [{"I_e": 200.0, "tau_m": 20.0},
                  {"I_e": 150.0, "tau_m": 30.0}]
epop3 = nest.Create("exc_iaf_neuron", 2, parameter_list)
# print nest.GetStatus(epop3, ['I_e', 'tau_m'])

# setting parameters for populations of neurons
Vth = -55.
Vrest = -70.
for neuron in epop1:
    nest.SetStatus([neuron], {"V_m": Vrest+(Vth-Vrest)*np.random.rand()})
# or
dVms = [{"V_m": Vrest+(Vth-Vrest)*np.random.rand()} for x in epop1]
nest.SetStatus(epop1, dVms)
# or
Vms = Vrest+(Vth-Vrest)*np.random.rand(len(epop1))
nest.SetStatus(epop1, "V_m", Vms)


# Specifying the behaviour of devices
# http://www.nest-simulator.org/helpindex/cc/RecordingDevice.html
# multimeter
multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime": True, "record_from": ["V_m"]})

recdict = {
    "to_memory": False,
    "to_file": True,
    "label": "epop_mp"}  # label is output filename

mm1 = nest.Create("multimeter", params=recdict)

# spikedetector
spikedetector = nest.Create("spike_detector", params={
                            "withgid": True, "withtime": True})

# Connect
# Connect(pre, post, conn_spec=None, syn_spec=None, model=None)
nest.Connect(multimeter, neuron)
nest.Connect(multimeter, pop2)
nest.Connect(neuron, spikedetector)
nest.Connect(neuron1, neuron2, syn_spec={"weight": 20, "delay": 1.0})
# default rule is all_to_all
nest.Connect(pop1, pop2, syn_spec={"weight": 20.0})
nest.Connect(pop1, pop2, 'one_to_one', syn_spec={'weight': 20.0, 'delay': 1.0})


# fixed_indegree
# fixed_outdegree
# fixed_total_number
# pairwise_bernoulli
epop1 = nest.Create("iaf_psc_alpha", 10)
ipop1 = nest.Create("iaf_psc_alpha", 10)
conn_dict_ex = {"rule": "fixed_indegree", "indegree": Ke}
conn_dict_in = {"rule": "fixed_indegree", "indegree": Ki}
conn_dict_in = {"rule": "fixed_outdegree", "outdegree": Ki}
conn_dict_in = {"rule": "fixed_total_number", "N": 10}
conn_dict_in = {"rule": "pairwise_bernoulli",
                "p": 0.1,
                "autapses": False,
                "multapses": False}

syn_dict_ex = {"delay": d, "weight": Je}
syn_dict_in = {"delay": d, "weight": Ji}
nest.Connect(epop1, ipop1, conn_dict_ex, syn_dict_ex)
nest.Connect(ipop1, epop1, conn_dict_in, syn_dict_in)
nest.Connect(ipop1, epop1,
             conn_spec=conn_dict_in,
             syn_spec=syn_dict_in)


n = 10
A = Create("iaf_psc_alpha", n)
B = Create("iaf_psc_alpha", n)
CopyModel("static_synapse", "excitatory", {"weight": 2.5, "delay": 0.5})
Connect(A, B, syn_spec="excitatory")


# distribution of synapse parameters
pop1 = nest.Create("iaf_psc_alpha", 50)
pop2 = nest.Create("iaf_psc_alpha", 50)

w_min = 0.5
w_max = 5.
syn_dict = {"model": "static_synapse",
            "weight": {"distribution": "uniform", "low": w_min, "high": w_max},
            "delay": 1.0}
syn_dict = {"model": "static_synapse",
            "weight": {"distribution": "normal", "mu": 0., "sigma": 1.},
            "delay": 1.0}
syn_dict = {"model": "static_synapse",
            "weight": {"distribution": "binomial", "n": 50, "p": .1},
            "delay": 1.0}
syn_dict = {"model": "static_synapse",
            "weight": {"distribution": "gamma", "order": 2., "scale": 2.},
            "delay": 1.0}
nest.Connect(pop1, pop2, "all_to_all", syn_dict)
conns = nest.GetConnections(pop1)
weights = []
for i in nest.GetStatus(conns, ['source', 'target', 'weight']):
    weights.append(i[2])


# Extracting and plotting data from devices
dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["V_m"]
ts = dmm["events"]["times"]
# pl.plot(ts,Vms)

# print nest.GetStatus(spikedetector)[0].keys()
# print nest.GetStatus(spikedetector)[0]['events']
dSD = nest.GetStatus(spikedetector, keys='events')[0]
evs = dSD['senders']
ts = dSD["times"]
# pl.plot(ts,evs,'.')

# get default values of models and synapses
nest.GetDefaults('iaf_psc_alpha')['recordables']
nest.GetDefaults('iaf_psc_alpha')


# poisson generator
noise_ex = nest.Create("poisson_generator")
noise_in = nest.Create("poisson_generator")
nest.SetStatus(noise_ex, {"rate": 80000.0})  # in [Hz]
nest.SetStatus(noise_in, {"rate": 15000.0})

# excitatory postsynaptic current of 1.2pA amplitude
syn_dict_ex = {"weight": 1.2}
# inhibitory postsynaptic current of -2pA amplitude
syn_dict_in = {"weight": -2.0}
nest.Connect(noise_ex, neuron, syn_spec=syn_dict_ex)
nest.Connect(noise_in, neuron, syn_spec=syn_dict_in)


# Efficient way to manage large data from recorders
mm = nest.Create('multimeter', params={'record_from': ['V_m', 'g_ex', 'g_in']})
nest.Simulate(2000)
res.to_pickle('data.pkl')
res.loc[123]
# would give you all lines in the dataframe with GID 123 as sender.
# Since 'senders' are set as dataframe index, this should be efficient
pd.read_csv()  # to load from disk
