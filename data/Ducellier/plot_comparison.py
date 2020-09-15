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
from math import floor

params = {'legend.fontsize': 24, \
          'xtick.labelsize':24, \
          'ytick.labelsize':24}
pylab.rcParams.update(params)

# List of LFE families
templates = np.loadtxt('../Plourde_2015/templates_list.txt', \
    dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
    'eZ', 'nb'), \
         'formats': ('S13', 'S3', np.float, np.float, np.float, \
    np.float, np.float, np.int)}, \
    skiprows=1)

# Beginning and end of the period we are looking at
tbegin = datetime(2007, 7, 1, 0, 0, 0)
tend = datetime(2009, 7, 1, 0, 0, 0)

# We construct the time series by counting the number of LFEs
# per one-day-long time window
window = 86400.0

# Length of the time series
dt = tend - tbegin
duration = dt.days * 86400.0 + dt.seconds + dt.microseconds * 0.000001
nw = int(duration / window)

# Loop on templates
for i in range(0, 1): #np.shape(templates)[0]):

    # Open LFE catalog
    namedir = 'catalogs/' + templates[i][0].astype(str)
    namefile1 = namedir + '/catalog_2007_2009.pkl'
    namefile2 = namedir + '/catalog_2002_2011.pkl'
    df1 = pickle.load(open(namefile1, 'rb'))
    df2 = pickle.load(open(namefile2, 'rb'))

    # Filter LFE
    df1 = df1.loc[df1['cc'] * df1['nchannel'] >= 1.1]
    df2 = df2.loc[df2['cc'] * df2['nchannel'] >= 1.5]
    # Get time series
    X1 = np.zeros(nw, dtype=int)
    # Loop on LFEs
    for j in range(0, len(df1)):
        myYear = df1['year'].iloc[j]
        myMonth = df1['month'].iloc[j]
        myDay = df1['day'].iloc[j]
        myHour = df1['hour'].iloc[j]
        myMinute = df1['minute'].iloc[j]
        mySecond = int(floor(df1['second'].iloc[j]))
        myMicrosecond = int(1000000.0 * (df1['second'].iloc[j] - mySecond))
        t = datetime(myYear, myMonth, myDay, myHour, myMinute, mySecond, \
            myMicrosecond)
        # Add LFE to appropriate time window
        if ((tbegin <= t) and (t < tbegin + timedelta(seconds=nw * window))):
            dt = t - tbegin
            duration = dt.days * 86400.0 + dt.seconds + dt.microseconds * \
                0.000001
            index = int(duration / window)
            X1[index] = X1[index] + 1
    print(np.sum(X1[X1 == 1]))

    X2 = np.zeros(nw, dtype=int)
    # Loop on LFEs
    for j in range(0, len(df2)):
        myYear = df2['year'].iloc[j]
        myMonth = df2['month'].iloc[j]
        myDay = df2['day'].iloc[j]
        myHour = df2['hour'].iloc[j]
        myMinute = df2['minute'].iloc[j]
        mySecond = int(floor(df2['second'].iloc[j]))
        myMicrosecond = int(1000000.0 * (df2['second'].iloc[j] - mySecond))
        t = datetime(myYear, myMonth, myDay, myHour, myMinute, mySecond, \
            myMicrosecond)
        # Add LFE to appropriate time window
        if ((tbegin <= t) and (t < tbegin + timedelta(seconds=nw * window))):
            dt = t - tbegin
            duration = dt.days * 86400.0 + dt.seconds + dt.microseconds * \
                0.000001
            index = int(duration / window)
            X2[index] = X2[index] + 1
    print(np.sum(X2[X2 == 1]))

    # Plot figure
    plt.figure(1, figsize=(16, 16))
    plt.subplot(211)
    plt.stem(np.arange(0, len(X1)), X1, 'k-', markerfmt=' ', basefmt=' ')
    plt.xlim([-0.5, len(X1) - 0.5])
    plt.xlabel('Time (days) since 2007/07/01', fontsize=24)
    plt.ylabel('Number of LFEs', fontsize=24)
    plt.title('Family {} ({:d} LFEs) - FAME stations'.format(templates[i][0].astype(str), np.sum(X1)), \
        fontsize=24)
    plt.subplot(212)
    plt.stem(np.arange(0, len(X2)), X2, 'k-', markerfmt=' ', basefmt=' ')
    plt.xlim([-0.5, len(X2) - 0.5])
    plt.xlabel('Time (days) since 2007/07/01', fontsize=24)
    plt.ylabel('Number of LFEs', fontsize=24)
    plt.title('Family {} ({:d} LFEs) - Permanent networks'.format(templates[i][0].astype(str), np.sum(X2)), \
        fontsize=24)
    plt.savefig('comparison/' + templates[i][0].astype(str) + '_daily_LFEs.eps', format='eps')
    plt.tight_layout()
    plt.close(1)
