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
xmin = 2004.0
xmax = 2012.0

# Space boundaries
lonmin = -123.1
lonmax = -122.2

# Start figure
params = {'legend.fontsize': 24, \
          'xtick.labelsize':24, \
          'ytick.labelsize':24}
pylab.rcParams.update(params)
plt.figure(1, figsize=(20, 10))

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

# Loop on templates
for i in [12, 13, 20, 23, 25]: #range(0, np.shape(templates)[0]):

    # Plot only good data
    if threshold['threshold_perm'].iloc[i] > 0.0:

        # Get latitude
        latitude = templates[i][2]
        longitude = templates[i][3]

        # Filter geography
        if (latitude >= 40.2) and (latitude <= 40.8):
            print(templates[i][0].astype(str))

            # Open LFE catalog
            namedir = 'catalogs/' + templates[i][0].astype(str)
            namefile = namedir + '/catalog_2004_2011.pkl'
            df = pickle.load(open(namefile, 'rb'))

            # Filter LFEs
            df = df.loc[df['cc'] * df['nchannel'] >= threshold['threshold_perm'].iloc[i]]

            time = np.arange(xmin, xmax, 1.0/365.5)
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
                index = int((t - xmin) * 365.5)
                if (index >= 0) and (index < len(nbLFEs)):
                    nbLFEs[index] = nbLFEs[index] + 1

            # Plot LFEs
            time = time[nbLFEs > 0]
            nbLFEs = nbLFEs[nbLFEs > 0]
            plt.scatter(time, np.repeat(longitude, np.shape(time)[0]), s=1 + nbLFEs * 5, c='k')

# End figure
plt.xlim([xmin, xmax])
plt.ylim([lonmin, lonmax])
plt.xlabel('Time (years)', fontsize=24)
plt.ylabel('Longitude', fontsize=24)
plt.title('LFEs', fontsize=24)
plt.tight_layout()
plt.savefig('LFEdistribution_perm/nbLFE_longitude.png', format='png')
plt.close(1)
