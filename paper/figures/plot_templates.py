import obspy

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np

import pickle

params = {'xtick.labelsize':20, \
          'ytick.labelsize':20}
pylab.rcParams.update(params)
plt.figure(1, figsize=(24, 6))
namedir = '../../src/templates_2007_2009/080401.05.050/'
station = 'WDC'
channels = ['BHE', 'BHN', 'BHZ']
for (i, channel) in enumerate(channels):
    filename = namedir + station + '_' + channel + '.pkl'
    stack = pickle.load(open(filename, 'rb'))[0]
    dt = stack.stats.delta
    nt = stack.stats.npts
    t = dt * np.arange(0, nt)
    plt.subplot2grid((1, len(channels)), (0, i))
    plt.plot(t, stack.data, 'k')
    plt.xlim([15.0, 35.0])
    plt.ylim([-1.5, 1.5])
    plt.title('{} - {}'.format(station, channel), fontsize=20)
    plt.xlabel('Time (s)', fontsize=20)
plt.tight_layout()
plt.savefig('templates.png', format='png')
plt.close(1)
