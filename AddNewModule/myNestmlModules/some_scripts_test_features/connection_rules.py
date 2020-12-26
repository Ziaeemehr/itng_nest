import nest
import pylab as pl


# nest.Install("testModuleInstallmodule")
nest.Install("models_module")
modelName = "hill_tononi_nestml"

dt = 0.01
tfinal = 100.0

nest.SetKernelStatus({
    "resolution": dt})

neurons = nest.Create(modelName, 5)
parameters = nest.GetDefaults(modelName)

if 0:
    for i in parameters:
        print(i, parameters[i])

# nest.SetStatus(neurons, {'I_e': 380.0})
conn_dict = {'rule': 'fixed_indegree', 'indegree': 2,
             'receptor_type': 1}
nest.Connect(neurons, neurons, conn_dict)

nest.Simulate(tfinal)
