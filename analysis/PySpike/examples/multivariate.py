""" Example for the multivariate spike distance

Copyright 2014, Mario Mulansky <mario.mulansky@gmx.net>

"""
# from __future__ import print_function
import pyspike as spk
from time import clock


def time_diff_in_ms(start, end):
    """ Returns the time difference end-start in ms.
    """
    return (end-start)*1000


t_start = clock()

# load the data
time_loading = clock()
spike_trains = spk.load_spike_trains_from_txt("PySpike_testdata.txt",
                                              edges=(0, 4000))
t_loading = clock()

print("Number of spike trains: %d" % len(spike_trains))

num_of_spikes = sum([len(spk) for spk in spike_trains])
print("Number of spikes: %d" % num_of_spikes)

# calculate the multivariate spike distance
f = spk.spike_profile(spike_trains)

t_spike = clock()

# print the average
avrg = f.avrg()
print("Spike distance from average: %.8f" % avrg)

t_avrg = clock()

# compute average distance directly, should give the same result as above
spike_dist = spk.spike_distance(spike_trains)
print("Spike distance directly:     %.8f" % spike_dist)

t_dist = clock()

print("Loading:            %9.1f ms" % time_diff_in_ms(t_start, t_loading))
print("Computing profile:  %9.1f ms" % time_diff_in_ms(t_loading, t_spike))
print("Averaging:          %9.1f ms" % time_diff_in_ms(t_spike, t_avrg))
print("Computing distance: %9.1f ms" % time_diff_in_ms(t_avrg, t_dist))
print("Total:              %9.1f ms" % time_diff_in_ms(t_start, t_dist))
