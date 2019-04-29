import scipy.io as sio
import pylab as pl 
from sys import exit
from scipy.signal import welch , filtfilt
from scipy.signal import butter, hilbert
import numpy as np 
import os


#----------------------------------------------#
def plot_lfp(time, lfp, ar_spikes, fname='f'):

	fig, (ax1, ax2) = pl.subplots(2, figsize=(8,6))

	f , PWelch_spec = welch(lfp ,fs=sf)
	filtered = filtering(lfp, sf, 10, 20, 5)

	ax1.plot(time, lfp, label="lfp")
	ax1.plot(time, filtered, label='filter')


	ax1.plot(ar_spikes, [0]*len(ar_spikes), 'ro', markersize=3, label='spike')

	ax1.set_xlabel("Time(ms)")
	ax1.set_ylabel(r"$LFP(\mu V)$")
	ax1.legend()

	ax2.set_yscale('log')
	ax2.plot(f, PWelch_spec)
	pl.savefig("fig/lfp/"+fname+'.png')
	pl.close()
#----------------------------------------------#
def power_spectrum_welch(time, lfp, ar_spikes):
	f , PWelch_spec = welch(lfp ,fs=sf)
	return f, PWelch_spec
#----------------------------------------------#

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def filtering(ar_signal, sampling_rate, lowcut, highcut, order=5):
	b, a = butter_bandpass(lowcut, highcut, sampling_rate, order)

	return filtfilt(b, a, ar_signal)
#----------------------------------------------#

PLOT_LFP_SIGNALS = True
PLOT_FILTERED = False


# loading data -------------------------------------------#
spikes = sio.loadmat('matfiles/spikes.mat')
spikes = spikes['spike_cell'][:,0]
lfp = sio.loadmat('matfiles/lfps.mat') #, squeez_me=True
time = lfp['time'][0]
lfpmat = lfp['lfp_matrix']
sf = lfp['sf'][0][0]
num_trial = 40


# ploting lfp signals ------------------------------------#
if PLOT_LFP_SIGNALS:
	for i in range(num_trial):
		plot_lfp(time, lfpmat[i,:], spikes[i][0], str("%02d"%i))
	



P_ave = 0.0
list_phase_spikes = []
for i in range(num_trial):
	f, P = power_spectrum_welch(time, lfpmat[i,:], spikes[i][0])
	P_ave += P

	filtered = filtering(lfpmat[i,:], sf, 10, 20, 5)
	H_fileterd = hilbert(filtered)
	H_norm = np.abs(H_fileterd)
	H_angle = np.angle(H_fileterd)
	indx_spike = ((spikes[i][0]-time[0])/1000*sf).astype(int)
	phase_spikes = H_angle[indx_spike].tolist()
	list_phase_spikes +=  phase_spikes

	if PLOT_FILTERED:
		fig , ax = pl.subplots(1, figsize=(10,6))
		ax.plot(time, filtered, label='filter')
		ax.plot(time, H_norm, label='norm_hilbert')
		ax.plot(time, H_angle, label='H_angle', c="k")
		ax.legend()
		# pl.ylim(-15,15)
		pl.savefig("./fig/f-"+str(i)+".png")
		pl.close()


l = np.asarray(list_phase_spikes)
n = len(list_phase_spikes)
plv = 1./float(n)* np.abs(np.sum(np.exp(1j*l)))
print "PLV = %g" % plv

pl.figure()
hist, bins = np.histogram(list_phase_spikes, bins=10)
width = (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
pl.bar(center, hist, align='center', width=width)
pl.savefig("fig/PLV.png")
pl.close()

P_ave /= float(num_trial)
fig2 = pl.figure()
pl.yscale('log')
pl.plot(f, P_ave)
pl.ylabel('Power A.U')
pl.xlabel("frequency(HZ)")
pl.savefig('./fig/ave.png')


















# print spikes[0][0]
# print spikes
# print len(spikes[0,0][0])
# print spikes[0,0][0]

# exit(0)

# print sf, len(time[0]), len(lfpmat)
# print type(lfpmat), lfpmat.shape


# pl.plot(PWelch_spec)
# pl.show()

# exit(0)
# print sf, len(time), lfpmat.shape
# print len(spikes), type(spikes), spikes.shape