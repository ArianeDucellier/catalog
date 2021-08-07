"""
Script to plot tremor and LFEs as a function of latitude
"""

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

from math import floor

from date import ymdhms2day

# Time boundaries
xmin = [2004.39, 2004.81, 2005.66, 2006.36, 2007.48, 2008.25]
xmax = [2004.45, 2004.88, 2005.73, 2006.43, 2007.59, 2008.35]
years = [2004, 2004, 2005, 2006, 2007, 2008]

# Space boundaries
latmin = 39.4
latmax = 41.8

# Read tremor file
tremor = pd.read_csv('../tremor/tremor.xyddtt', sep=' ', header=None)
tremor.columns = ['longitude', 'latitude', 'downdip', 'alongstrike', 'year', 'epoch']
year_tremor = tremor['year']
lat_tremor = tremor['latitude']

# List of LFE families
templates = np.loadtxt('../Plourde_2015/templates_list.txt', \
    dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
    'eZ', 'nb'), \
         'formats': ('S13', 'S3', np.float, np.float, np.float, \
    np.float, np.float, np.int)}, \
    skiprows=1)

# Threshold for filtering the catalog
threshold = pd.read_csv('threshold_cc.txt', sep=r'\s{1,}', header=None, engine='python')
threshold.columns = ['family', 'threshold_FAME', 'threshold_perm']

# Start figure and plot tremor
params = {'legend.fontsize': 24, \
          'xtick.labelsize':24, \
          'ytick.labelsize':24}
pylab.rcParams.update(params)
plt.figure(1, figsize=(30, 10))

for event in range(0, 6):
    ax = plt.subplot2grid((1, 6), (0, event))
    plt.scatter(year_tremor - years[event], lat_tremor, c='r')

    # Loop on templates
    for i in range(0, np.shape(templates)[0]):

        # Plot only good data
        if threshold['threshold_perm'].iloc[i] > 0.0:

            # Get latitude
            latitude = templates[i][2]

            # Open LFE catalog
            namedir = 'catalogs/' + templates[i][0].astype(str)
            namefile = namedir + '/catalog_2004_2011.pkl'
            df = pickle.load(open(namefile, 'rb'))

            # Filter LFEs
            df = df.loc[df['cc'] * df['nchannel'] >= threshold['threshold_perm'].iloc[i]]

            time = np.arange(xmin[event], xmax[event], 1.0/365.5)
            nbLFEs = np.zeros(np.shape(time)[0])
            # Loop on LFEs
            for j in range(0, len(df)):
                myYear = df['year'].iloc[j]
                myMonth = df['month'].iloc[j]
                myDay = df['day'].iloc[j]
                myHour = df['hour'].iloc[j]
                myMinute = df['minute'].iloc[j]
                mySecond = int(floor(df['second'].iloc[j]))
                myMicrosecond = int(1000000.0 * (df['second'].iloc[j] - mySecond))
                t = ymdhms2day(myYear, myMonth, myDay, myHour, myMinute, mySecond)
                index = int((t - xmin[event]) * 365.5)
                if (index >= 0) and (index < len(nbLFEs)):
                    nbLFEs[index] = nbLFEs[index] + 1

            # Plot LFEs
            time = time[nbLFEs > 3]
            nbLFEs = nbLFEs[nbLFEs > 3]
            plt.scatter(time - years[event], np.repeat(latitude, np.shape(time)[0]), s=1 + nbLFEs * 5, c='k')

# Families
#plt.annotate('A', (xmax + 0.06 * (xmax - xmin), 40.09000), fontsize=24, color='grey')
#plt.annotate('B', (xmax + 0.06 * (xmax - xmin), 40.38000), fontsize=24, color='grey')
#plt.annotate('C', (xmax + 0.06 * (xmax - xmin), 40.65000), fontsize=24, color='grey')
#plt.annotate('D', (xmax + 0.06 * (xmax - xmin), 40.92000), fontsize=24, color='grey')
#plt.annotate('E', (xmax + 0.08 * (xmax - xmin), 40.97000), fontsize=24, color='grey')
#plt.annotate('F', (xmax + 0.06 * (xmax - xmin), 41.02000), fontsize=24, color='grey')
#plt.annotate('G', (xmax + 0.08 * (xmax - xmin), 41.10000), fontsize=24, color='grey')

    # End figure
    plt.xlim([xmin[event] - years[event], xmax[event] - years[event]]) # + 0.1 * (xmax - xmin)])
    plt.ylim([39.1, latmax])
    plt.xlabel('Time (years)', fontsize=24)
    if event == 0:
        plt.ylabel('Latitude', fontsize=24)
    plt.title('Year {:d}'.format(years[event]), fontsize=24)

plt.tight_layout()
plt.savefig('LFEdistribution_perm/zoom_ETS_events.eps', format='eps')
plt.close(1)
