import numpy as np
from matplotlib import pyplot as plt
import pyspike as spk


st1 = spk.generate_poisson_spikes(1.0, [0, 20])
st2 = spk.generate_poisson_spikes(1.0, [0, 20])
# spike_trains = spk.load_spike_trains_from_txt("PySpike_testdata.txt",
#                                               edges=(0, 4000))

# st1 = spike_trains[0]
# st2 = spike_trains[1]
d = spk.spike_directionality(st1, st2)

print "Spike Directionality of two Poissonian spike trains:", d

E = spk.spike_train_order_profile(st1, st2)
# E = spk.spike_train_order_profile(spike_trains)

plt.figure()
x, y = E.get_plottable_data()
plt.plot(x, y, '-ob')
plt.ylim(-1.1, 1.1)
plt.xlabel("t")
plt.ylabel("E")
plt.title("Spike Train Order Profile")
plt.savefig("fig/spk_order")
plt.close()

###### Optimize spike train order of 20 Random spike trains #######

M = 20

spike_trains = [spk.generate_poisson_spikes(1.0, [0, 100]) for m in xrange(M)]

F_init = spk.spike_train_order(spike_trains)

print "Initial Synfire Indicator for 20 Poissonian spike trains:", F_init

D_init = spk.spike_directionality_matrix(spike_trains)

phi, _ = spk.optimal_spike_train_sorting(spike_trains)

F_opt = spk.spike_train_order(spike_trains, indices=phi)

print "Synfire Indicator of optimized spike train sorting:", F_opt

D_opt = spk.permutate_matrix(D_init, phi)

plt.figure()
plt.imshow(D_init)
plt.title("Initial Directionality Matrix")

plt.figure()
plt.imshow(D_opt)
plt.title("Optimized Directionality Matrix")

plt.show()
