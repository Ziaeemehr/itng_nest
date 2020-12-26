import nest
import pylab as pl


modelName = "ht_neuron"

dt = 0.01
tfinal = 100.0

nest.SetKernelStatus({
    "resolution": dt})

neurons = nest.Create(modelName, 5)
parameters = nest.GetDefaults(modelName)

# if 0:
#     for i in parameters:
#         print(i, parameters[i])

# nest.SetStatus(neurons, {'I_e': 380.0})
conn_dict = {'rule': 'fixed_indegree', 'indegree': 2}
nest.Connect(neurons, neurons, syn_spec={'receptor_type': 1},
conn_spec = conn_dict)

nest.Simulate(tfinal)
