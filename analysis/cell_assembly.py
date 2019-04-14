import matplotlib.pyplot as plt 
import elephant.conversion as conv
import elephant.spike_train_generation 
import quantities as pq
import numpy as np
import elephant.cell_assembly_detection as cad
np.random.seed(30)
# Generate correlated data and bin it with a binsize of 10ms 
sts = elephant.spike_train_generation.cpp( 
    rate=15*pq.Hz, A=[0]+[0.95]+[0]*4+[0.05], t_stop=10*pq.s)
binsize = 10*pq.ms 
spM = conv.BinnedSpikeTrain(sts, binsize=binsize)
# Call of the method 
patterns = cad.cell_assembly_detection(spM, maxlag=2)[0] 
# plotting
plt.figure() 
for neu in patterns['neurons']: 
    if neu == 0:
        plt.plot(
        patterns['times']*binsize, [neu]*len(patterns['times']),
        'ro', label='pattern')
    else: 
        plt.plot(
        patterns['times']*binsize, [neu] * len(patterns['times']),
        'ro')
# Raster plot of the data 
for st_idx, st in enumerate(sts): 
    if st_idx == 0:
        plt.plot(st.rescale(pq.ms), [st_idx] * len(st), 'k.',
            label='spikes')
    else: 
        plt.plot(st.rescale(pq.ms), [st_idx] * len(st), 'k.')

plt.ylim([-1, len(sts)])
plt.xlabel('time (ms)')
plt.ylabel('neurons ids')
plt.show()