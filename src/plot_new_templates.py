import obspy

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle

from math import cos, pi, sin, sqrt

def plot_new_templates(family, directory, timearrival, latitude, longitude, ie):
    """
    """
    if timearrival == True:
        data = pickle.load(open('tori.pkl', 'rb'))
#        tori = data[0]
#        sS = data[1]
#        sP = data[2]
        tori = data
        staloc = pd.read_csv('../data/Ducellier/stations_permanent.txt', \
            sep=r'\s{1,}', header=None, engine='python')
        staloc.columns = ['station', 'network', 'channels', 'location', \
            'server', 'latitude', 'longitude', 'time_on', 'time_off']
        a = 6378.136
        e = 0.006694470

    params = {'xtick.labelsize':20, \
              'ytick.labelsize':20}
    pylab.rcParams.update(params)
    namedir = directory + '/' + family
#    templates = os.listdir(namedir)
    templates = ['WDC_BHN.pkl']
    for template in templates:
        name = template.split('.')[0]
        station = name.split('_')[0]
        channel = name.split('_')[1]
        filename = namedir + '/' + template
        stack = pickle.load(open(filename, 'rb'))[0]
        RMS = np.sqrt(np.mean(np.square(stack.data)))
        maximum = np.max(np.abs(stack.data))
        SNR = maximum / RMS
        if timearrival == True:
            subdf = staloc.loc[staloc['station'] == station]
            distance = 0.0
            if (len(subdf) > 0):
                lat = subdf['latitude'].iloc[0]
                lon = subdf['longitude'].iloc[0]
                dx = (pi / 180.0) * a * cos(latitude * pi / 180.0) / sqrt(1.0 - e * e * \
                    sin(latitude * pi / 180.0) * sin(latitude * pi / 180.0))
                dy = (3.6 * pi / 648.0) * a * (1.0 - e * e) / ((1.0 - e * e * sin(latitude * \
                    pi / 180.0) * sin(latitude * pi / 180.0)) ** 1.5)
                x = dx * (lon - longitude)
                y = dy * (lat - latitude)
                distance = sqrt(x ** 2.0 + y ** 2.0)
#                tS = tori[ie] + distance * sS
#                tP = tori[ie] + distance * sP
#                tS = tori[ie] + distance * 0.25
#                tP = tori[ie] + distance * 0.125
                tS = 23.5 # for WDC and 16.5 for LRB
                tP = 16.5 # for WDC and 13.0 for LRB 
                plot_time = True
            else:
                plot_time = False
#        plt.figure(1, figsize=(10, 5))
        plt.figure(1, figsize=(6, 5))
        dt = stack.stats.delta
        nt = stack.stats.npts
        t = dt * np.arange(0, nt)
        if timearrival == True:
            if plot_time == True:
                plt.axvline(tS, linewidth=4, color='grey')
                plt.axvline(tP, linewidth=4, color='grey')
        plt.plot(t, stack.data, 'k')
#        plt.xlim([np.min(t), np.max(t)])
        plt.xlim([10.0, 30.0])
        if timearrival == True:
            plt.title('{} - {} - SNR = {:.2f} - distance = {:.2f} km'. \
                format(station, channel, SNR, distance), fontsize=20)
        plt.title('{} - {} - SNR = {:.2f}'.format(station, channel, SNR), fontsize=20)
        plt.xlabel('Time (s)', fontsize=20)
        plt.tight_layout()
        plt.savefig(namedir + '/' + station + '_' + channel + '.eps', format='eps')
        plt.close(1)

if __name__ == '__main__':

    # Set the parameters
    directory = 'templates_2007_2009'

    LFEloc = np.loadtxt('../data/Plourde_2015/templates_list.txt', \
        dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
        'eZ', 'nb'), \
             'formats': ('S13', 'S3', np.float, np.float, np.float, \
        np.float, np.float, np.int)}, \
        skiprows=1)

    for ie in range(0, 1): #len(LFEloc)):
        family = LFEloc[ie][0].decode('utf-8')
        latitude = LFEloc[ie][2]
        longitude = LFEloc[ie][3]
        plot_new_templates(family, directory, True, latitude, longitude, ie)
