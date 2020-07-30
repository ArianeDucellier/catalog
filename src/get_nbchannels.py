"""
"""

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from datetime import datetime

def get_nb_channels(family, station_file, stations, Ystart, Yend):
    """
    """
    staloc = pd.read_csv(station_file, sep=r'\s{1,}', header=None, engine='python')
    staloc.columns = ['station', 'network', 'channels', 'location', \
        'server', 'latitude', 'longitude', 'Start time', 'End time', 'dt']
    startdate = pd.to_datetime(staloc['Start time'], format='%Y-%m-%d')
    enddate = staloc['End time'].str.replace('3000', '2025')
    staloc['End time'] = enddate
    enddate = pd.to_datetime(staloc['End time'], format='%Y-%m-%d')
    staloc['Start time'] = startdate
    staloc['End time'] = enddate

    nmonths = (Yend - Ystart + 1) * 12
    nchannels = np.zeros(nmonths)
    
    for year in range(Ystart, Yend + 1):
        for month in range(1, 13):
            mydate = datetime(year, month, 15)
            for station in stations:
                subset = staloc.loc[staloc['station'] == station]
                if (min(subset['Start time']) < mydate) & (max(subset['End time']) > mydate):
                    channels = subset['channels'].iloc[0]
                    mychannels = channels.split(',')
                    index = (year - Ystart) * 12 + month - 1
                    nchannels[index] = nchannels[index] + len(mychannels)

    plt.figure(1, figsize=(10, 6))
    params = {'xtick.labelsize':16,
              'ytick.labelsize':16}
    pylab.rcParams.update(params)
    time = Ystart + np.arange(0, nmonths) / 12.0
    plt.plot(time, nchannels)
    plt.xlabel('Year', fontsize=24)
    plt.ylabel('Number of channels', fontsize=24)
    plt.title('family ' + family, fontsize=30)
    plt.savefig('../data/Ducellier/timeline/' + family + '.eps', fromat='eps')
    plt.close(1)
                
if __name__ == '__main__':

    family_file = '../data/Ducellier/families_permanent.txt'
    station_file = '../data/Ducellier/stations_permanent.txt'
    Ystart = 1993
    Yend = 2013

    families = pd.read_csv(family_file, sep=r'\s{1,}', header=None, engine='python')
    families.columns = ['family', 'stations', 'duration']
    for i in range(0, len(families)):
        family = families['family'].iloc[i]
        stations = families['stations'].iloc[i].split(',')
        get_nb_channels(family, station_file, stations, Ystart, Yend)
