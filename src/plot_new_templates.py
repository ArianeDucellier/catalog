import obspy

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle

from math import cos, pi, sin, sqrt

from get_stations import filter_stations

def plot_new_templates(family, directory, timearrival, latitude, longitude, \
    ie, tbegin, tend):
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
            'server', 'latitude', 'longitude', 'time_on', 'time_off', 'dt']
    else:
        # Get the network, channels, and location of the stations
        stations_BK = filter_stations('BK')
        stations_NC = filter_stations('NC')
        stations_PB = filter_stations('PB')
        stations_XQ = filter_stations('FAME')
        stations = pd.concat([stations_BK, stations_NC, stations_PB, \
            stations_XQ], ignore_index=True)
    a = 6378.136
    e = 0.006694470

    # Start figure    
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
        plot_time = False
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
                tS = tori[ie] + distance * 0.25
                tP = tori[ie] + distance * 0.125
#                tS = 23.5 # for WDC and 16.5 for LRB
#                tP = 16.5 # for WDC and 13.0 for LRB 
                plot_time = True
        else:
            subdf = stations.loc[stations['Stat'] == station]
            distance = 0.0
            if (len(subdf) > 0):
                lat = subdf['Lat'].iloc[0]
                lon = subdf['Lon'].iloc[0]
                dx = (pi / 180.0) * a * cos(latitude * pi / 180.0) / sqrt(1.0 - e * e * \
                    sin(latitude * pi / 180.0) * sin(latitude * pi / 180.0))
                dy = (3.6 * pi / 648.0) * a * (1.0 - e * e) / ((1.0 - e * e * sin(latitude * \
                    pi / 180.0) * sin(latitude * pi / 180.0)) ** 1.5)
                x = dx * (lon - longitude)
                y = dy * (lat - latitude)
                distance = sqrt(x ** 2.0 + y ** 2.0)
                
        plt.figure(1, figsize=(10, 5))
        dt = stack.stats.delta
        nt = stack.stats.npts
        t = dt * np.arange(0, nt)
        if timearrival == True:
            if plot_time == True:
                plt.axvline(tS, linewidth=4, color='grey')
                plt.axvline(tP, linewidth=4, color='grey')
#        plt.axvline(tbegin, linewidth=4, color='red')
#        plt.axvline(tend, linewidth=4, color='red')
        plt.plot(t, stack.data, 'k')
        plt.xlim([np.min(t), np.max(t)])
        if (len(subdf) > 0):
            plt.title('{} - {} - SNR = {:.2f} - distance = {:.2f} km'. \
                format(station, channel, SNR, distance), fontsize=20)
        else:
            plt.title('{} - {} - SNR = {:.2f}'.format(station, channel, SNR), fontsize=20)
        plt.xlabel('Time (s)', fontsize=20)
        plt.tight_layout()
        plt.savefig(namedir + '/' + station + '_' + channel + '.eps', format='eps')
        plt.close(1)

if __name__ == '__main__':

    # Set the parameters
    directory = 'templates_both'

    LFEloc = np.loadtxt('../data/Plourde_2015/templates_list.txt', \
        dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
        'eZ', 'nb'), \
             'formats': ('S13', 'S3', np.float, np.float, np.float, \
        np.float, np.float, np.int)}, \
        skiprows=1)

#    families = pd.read_csv('../data/Ducellier/families_permanent.txt', \
#        sep=r'\s{1,}', header=None, engine='python')
#    families.columns = ['family', 'stations', 'tbegin', 'tend']

    for ie in range(0, 1): #len(LFEloc)):
        family = LFEloc[ie][0].decode('utf-8')
#        for i in range(0, len(families)):
#            if families['family'][i] == family:
#                tbegin = families['tbegin'][i]
#                tend = families['tend'][i]
        tbegin = 0.0
        tend = 0.0
        latitude = LFEloc[ie][2]
        longitude = LFEloc[ie][3]
        plot_new_templates(family, directory, False, latitude, longitude, \
            ie, tbegin, tend)
