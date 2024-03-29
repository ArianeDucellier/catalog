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
tbegin = datetime(2004, 1, 1, 0, 0, 0)
tend = datetime(2012, 1, 1, 0, 0, 0)

# Threshold for filtering the catalog
threshold = pd.read_csv('threshold_cc.txt', sep=r'\s{1,}', header=None, engine='python')
threshold.columns = ['family', 'threshold_FAME', 'threshold_perm']

# We construct the time series by counting the number of LFEs
# per one-day-long time window
window = 86400.0

# Length of the time series
dt = tend - tbegin
duration = dt.days * 86400.0 + dt.seconds + dt.microseconds * 0.000001
nw = int(duration / window)

# Indices of LFE families
#indices = [[2, 5, 6], \
#           [12, 13, 25, 20, 23], \
#           [40, 29, 33, 36, 41]]
indices = [[40, 29, 33, 36, 41]]

# Names of LFE families
#names = [['B1', 'B2', 'B3'], \
#         ['C1', 'C2', 'C3', 'C4', 'C5'], \
#         ['G1', 'G2', 'G3', 'G4', 'G5']]
names = [['G1', 'G2', 'G3', 'G4', 'G5']]

#plt.figure(1, figsize=(24, 15))
plt.figure(1, figsize=(8, 15))

# Loop on templates
for column in range(0, len(indices)):
    for (count, i) in enumerate(indices[column]):
        longitude = templates[i][3]

        # Open LFE catalog
        namedir = 'catalogs/' + templates[i][0].astype(str)
        namefile = namedir + '/catalog_2004_2011.pkl'
        df = pickle.load(open(namefile, 'rb'))

        # Filter LFE
        df = df.loc[df['cc'] * df['nchannel'] >= threshold['threshold_perm'].iloc[i]]

        # Get time series
        X = np.zeros(nw, dtype=int)
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
            # Add LFE to appropriate time window
            if ((tbegin <= t) and (t < tbegin + timedelta(seconds=nw * window))):
                dt = t - tbegin
                duration = dt.days * 86400.0 + dt.seconds + dt.microseconds * \
                    0.000001
                index = int(duration / window)
                X[index] = X[index] + 1

        # Plot figure
        plt.subplot2grid((5, len(indices)), (4 - count, column))
        plt.plot((1, 1))
        plt.stem(2004 + np.arange(0, len(X)) / 365.25, X, 'k-', markerfmt=' ', basefmt=' ')
        plt.xlim([2004, 2012])
        if count == 0:
            plt.xlabel('Time (years)', fontsize=24)
        if column == 0:
            plt.ylabel('Number of LFEs', fontsize=24)
        plt.title('Family {} at {}'.format(names[column][count], longitude), \
            fontsize=24)

plt.tight_layout()
#plt.savefig('LFEdistribution_perm/set_daily_LFEs.eps', format='eps')
plt.savefig('LFEdistribution_perm/set_daily_LFEs_seriesG.eps', format='eps')
plt.close(1)
