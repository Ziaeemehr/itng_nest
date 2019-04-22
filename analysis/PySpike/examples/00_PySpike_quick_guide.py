import pyspike as spk
import matplotlib.pyplot as plt



a = np.random.rand(10)*10
b = spk.SpikeTrain(a, [0,10], is_sorted=False)
print b.spikes

spike_trains = spk.load_spike_trains_from_txt("PySpike_testdata.txt",
                                              edges=(0, 4000))

# compute the two spike trains and multivariate ISI profile
f = spk.isi_profile(spike_trains[0], spike_trains[1])
f = spk.isi_profile(spike_trains)

# t = [900, 1100, 2000, 3100]
# print("ISI value at t =", t, ":", f(t))
# print("Average ISI distance:", f.avrg())

# compute the two spike trains and multivariate SPIKE profile
f = spk.spike_profile(spike_trains[0], spike_trains[1])
f = spk.spike_profile(spike_trains)

# t = [900, 1100, 2000, 3100]
# print("Multivariate SPIKE value at t =", t, ":", f(t))
# print("Average multivariate SPIKE distance:", f.avrg())


# plot the spike times
for (i, spike_train) in enumerate(spike_trains):
    plt.scatter(spike_train, i*np.ones_like(spike_train), marker='|')
    # print np.asarray(spike_train)

# profile of the first two spike trains
f = spk.isi_profile(spike_trains, indices=[0, 1])
x, y = f.get_plottable_data() 
# x = f.x
# y = f.y
plt.plot(x, np.abs(y), '--k', label="ISI-profile")


# spike distance profile
f = spk.spike_profile(spike_trains, indices=[0, 1])
f = spk.spike_profile(spike_trains[0], spike_trains[1])
f = spk.spike_profile(spike_trains)
x, y = f.get_plottable_data()
plt.plot(x, y, '-b', label="SPIKE-profile")

avrg = f.avrg()
print("Spike distance from average: %.8f" % avrg)

# compute average distance directly, should give the same result as above
spike_dist = spk.spike_distance(spike_trains)
print("Spike distance directly:     %.8f" % spike_dist)


# spike synchronization profile
f = spk.spike_sync_profile(spike_trains[0], spike_trains[1])
f = spk.spike_sync_profile(spike_trains)
x, y = f.get_plottable_data()
plt.plot(x, y, '-b', alpha=0.7, label="SPIKE-Sync profile")

x1, y1 = f.get_plottable_data(averaging_window_size=50)
plt.plot(x1, y1, '-k', lw=2.5, label="averaged SPIKE-Sync profile")
# The optional parameter averaging_window_size determines the size 
# of an averaging window to smoothen the profile. If this value is
# 0, no averaging is performed.
print("Average:", f.avrg())

# For the direct computation of the overall spike synchronization 
# value within some interval, the spike_sync() function can be used:
spike_sync = spk.spike_sync(spike_trains[0], spike_trains[1], interval=ival)


# Computes the peri-stimulus time histogram
# The PSTH is simply the histogram of merged spike events.
f_psth = spk.psth(spike_trains, bin_size=50.0)
x, y = f_psth.get_plottable_data()
plt.plot(x, y, '-k', alpha=1.0, label="PSTH")


# print("Number of spike trains: %d" % len(spike_trains))
# num_of_spikes = sum([len(spk) for spk in spike_trains])
# print("Number of spikes: %d" % num_of_spikes)

# spike_train_order.
# Generates a Poisson spike train with the given 
# rate in the given time interval
st1 = spk.generate_poisson_spikes(1.0, [0, 20])
st2 = spk.generate_poisson_spikes(1.0, [0, 20])

d = spk.spike_directionality(st1, st2)
E = spk.spike_train_order_profile(st1, st2)
x, y = E.get_plottable_data()
plt.plot(x, y, '-ob')
plt.xlabel("t")
plt.ylabel("E")
plt.title("Spike Train Order Profile")

# Optimize spike train order of 20 Random spike trains
spike_trains = [spk.generate_poisson_spikes(1.0, [0, 100]) for m in xrange(20)]
F_init = spk.spike_train_order(spike_trains)
print "Initial Synfire Indicator for 20 Poissonian spike trains:", F_init
D_init = spk.spike_directionality_matrix(spike_trains)
phi, _ = spk.optimal_spike_train_sorting(spike_trains)
F_opt = spk.spike_train_order(spike_trains, indices=phi)
print "Synfire Indicator of optimized spike train sorting:", F_opt
plt.imshow(D_init)
plt.title("Initial Directionality Matrix")
plt.imshow(D_opt)
plt.title("Optimized Directionality Matrix")


# merging of two spike trains
merged_spike_train = spk.merge_spike_trains([spike_trains[0], spike_trains[1]])
print (merged_spike_train.spikes) # <type 'numpy.ndarray'>

plt.plot(spike_trains[0], np.ones_like(spike_trains[0]), 'o')
plt.plot(spike_trains[1], np.ones_like(spike_trains[1]), 'x')
plt.plot(merged_spike_train.spikes,
         2*np.ones_like(merged_spike_train), 'o')


# compute the isi distance matrix of a set of spike trains.
isi_distance = spk.isi_distance_matrix(spike_trains)
plt.imshow(isi_distance, interpolation='none')

spike_distance = spk.spike_distance_matrix(spike_trains, interval=(0, 1000))
plt.imshow(spike_distance, interpolation='none')

spike_sync = spk.spike_sync_matrix(spike_trains, interval=(2000, 4000))
plt.imshow(spike_sync, interpolation='none')
