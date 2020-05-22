import obspy

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle

def plot_new_templates(family, directory):
    """
    """
    params = {'xtick.labelsize':20, \
              'ytick.labelsize':20}
    pylab.rcParams.update(params)
    namedir = directory + '/' + family
    templates = os.listdir(namedir)
    for template in templates:
        name = template.split('.')[0]
        station = name.split('_')[0]
        channel = name.split('_')[1]
        filename = namedir + '/' + template
        stack = pickle.load(open(filename, 'rb'))[0]
        RMS = np.sqrt(np.mean(np.square(stack.data)))
        maximum = np.max(np.abs(stack.data))
        SNR = maximum / RMS
        plt.figure(1, figsize=(10, 5))
        dt = stack.stats.delta
        nt = stack.stats.npts
        t = dt * np.arange(0, nt)
        plt.plot(t, stack.data, 'k')
        plt.xlim([np.min(t), np.max(t)])
        plt.title('{} - {} - SNR = {:.2f}'.format(station, channel, SNR), fontsize=20)
        plt.xlabel('Time (s)', fontsize=20)
        plt.tight_layout()
        plt.savefig(namedir + '/' + station + '_' + channel + '.eps', format='eps')
        plt.close(1)

if __name__ == '__main__':

    # Set the parameters
    family = '080401.05.050'
    directory = 'templates_2007_2009'

    plot_new_templates(family, directory)
