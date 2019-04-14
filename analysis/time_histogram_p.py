from elephant.spike_train_generation import homogeneous_poisson_process
from quantities import Hz, s, ms
from elephant.statistics import time_histogram

spiketrain = homogeneous_poisson_process(rate=10.0*Hz, t_start=0.0*s, t_stop=100.0*s)
binsize=5.0*ms

time_histogram( spiketrain, binsize, t_start=0.0, t_stop=100.0*s)
