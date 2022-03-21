"""
This script looks at the catalog from southern Cascadia (2007-2009) and plots
for each family the daily number of LFEs in function of time
"""

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

from datetime import datetime, timedelta
from math import cos, floor, pi, sin, sqrt

params = {'legend.fontsize': 24, \
          'xtick.labelsize':24, \
          'ytick.labelsize':24}
pylab.rcParams.update(params)

# Earth's radius and ellipticity
#a = 6378.136
#e = 0.006694470

# Conversion lat-lon into kilometers
#lat0 = 40.281
#lon0 = -124.433
#dx = (pi / 180.0) * a * cos(lat0 * pi / 180.0) / sqrt(1.0 - e * e * \
#    sin(lat0 * pi / 180.0) * sin(lat0 * pi / 180.0))
#dy = (3.6 * pi / 648.0) * a * (1.0 - e * e) / ((1.0 - e * e * sin(lat0 * \
#    pi / 180.0) * sin(lat0 * pi / 180.0)) ** 1.5)

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

# Beginning and end of the period we are looking at
tbegin = datetime(2011, 3, 11, 5, 46, 24)
tend = datetime(2011, 3, 11, 11, 46, 24)
teq = datetime(2011, 3, 11, 5, 46, 24)

tmin = tbegin - teq
xmin = tmin.days * 24 + tmin.seconds / 3600.0
tmax = tend - teq
xmax = tmax.days * 24 + tmax.seconds / 3600.0

# Loop on templates
for i in range(0, np.shape(templates)[0]):

    # Plot only good data
    if threshold['threshold_perm'].iloc[i] > 0.0:

    # Distance to epicenter
#        latitude = templates[i][2]
#        longitude = templates[i][3]
#        x = (longitude - lon0) * dx
#        y = (latitude - lat0 ) * dy
#        dist = sqrt(x**2.0 + y**2.0)

        # Open LFE catalog
        namedir = 'catalogs/' + templates[i][0].astype(str)
        namefile = namedir + '/catalog_2004_2011.pkl'
        df = pickle.load(open(namefile, 'rb'))

        # Filter LFEs
        df = df.loc[df['cc'] * df['nchannel'] >= threshold['threshold_perm'].iloc[i]]

        # Create time series
        times = np.zeros(len(df))

        # Loop on LFEs
        for j in range(0, len(df)):
            myYear = df['year'].iloc[j]
            myMonth = df['month'].iloc[j]
            myDay = df['day'].iloc[j]
            myHour = df['hour'].iloc[j]
            myMinute = df['minute'].iloc[j]
            mySecond = int(floor(df['second'].iloc[j]))
            myMicrosecond = int(1000000.0 * (df['second'].iloc[j] - mySecond))
            t = datetime(myYear, myMonth, myDay, myHour, myMinute, mySecond, \
                myMicrosecond)
            dt = t - teq
            times[j] = dt.days * 24 + dt.seconds / 3600.0

        df['time'] = times

        # Plot figure
        plt.figure(1, figsize=(20, 10))
        plt.axvline(0, linewidth=2, color='red')
        plt.stem(df['time'], np.repeat(1, len(df)), 'k-', markerfmt=' ', basefmt=' ')
        plt.xlim([xmin, xmax])
        plt.xlabel('Time (hours) since 2011-03-09 at 05:46:24', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)
        plt.title('Family {}'.format(templates[i][0].astype(str)), fontsize=24)
#        plt.title('Family {} ({:4.2f} km from epicenter)'.format(templates[i][0].astype(str), dist), \
#            fontsize=24)
        plt.savefig('zoom_teq1_p6/' + templates[i][0].astype(str) + '.eps', format='eps')
        plt.close(1)
