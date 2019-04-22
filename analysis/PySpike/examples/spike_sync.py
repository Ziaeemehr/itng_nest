# from __future__ import print_function

import matplotlib.pyplot as plt
import pyspike as spk
from sys import exit

spike_trains = spk.load_spike_trains_from_txt("../test/SPIKE_Sync_Test.txt",
                                              edges=(0, 4000))

# print type(spike_trains), len(spike_trains), type(spike_trains[0])

plt.figure()
f = spk.spike_sync_profile(spike_trains[0], spike_trains[1])
x, y = f.get_plottable_data()
plt.plot(x, y, '--ok', label="SPIKE-SYNC profile")
# print(f.x)
# print(f.y)
# print(f.mp)
print("Average:", f.avrg())

f = spk.spike_profile(spike_trains[0], spike_trains[1])
x, y = f.get_plottable_data()

plt.plot(x, y, '-b', label="SPIKE-profile")

plt.axis([0, 4000, -0.1, 1.1])
plt.legend(loc="center right")
plt.savefig('fig/spike_sync_2.png')
plt.close()


# multivariant---------------------------------------------#
plt.figure()
plt.subplot(211)

f = spk.spike_sync_profile(spike_trains)
x, y = f.get_plottable_data()
plt.plot(x, y, '-b', alpha=0.7, label="SPIKE-Sync profile")

x1, y1 = f.get_plottable_data(averaging_window_size=50)
plt.plot(x1, y1, '-k', lw=2.5, label="averaged SPIKE-Sync profile")
plt.legend(loc='lower right')

f_direct = spk.spike_sync(spike_trains) #interval=[0,4000]

print "Average from spike_sync_profile :", f.avrg()
print "directly from spike_sync : ", f_direct


plt.subplot(212)
f_psth = spk.psth(spike_trains, bin_size=50.0)
x, y = f_psth.get_plottable_data()
plt.plot(x, y, '-k', alpha=1.0, label="PSTH")
plt.legend()


plt.savefig('fig/spike_sync_multi.png')
plt.close()

# plt.show()
