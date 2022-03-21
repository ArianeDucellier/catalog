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
#xmin = 2004.0
#xmax = 2012.0

xmin = (2008, 4, 27, 6, 28, 1)
xmax = (2008, 5, 27, 6, 28, 1)
xmin = ymdhms2day(xmin[0], xmin[1], xmin[2], xmin[3], xmin[4], xmin[5])
xmax = ymdhms2day(xmax[0], xmax[1], xmax[2], xmax[3], xmax[4], xmax[5])

# Space boundaries
latmin = 39.4
latmax = 41.8

# Read tremor file
tremor = pd.read_csv('../tremor/tremor.xyddtt', sep=' ', header=None)
tremor.columns = ['longitude', 'latitude', 'downdip', 'alongstrike', 'year', 'epoch']
tremor = tremor.loc[(tremor['year'] >= xmin) & (tremor['year'] <= xmax)]
year_tremor = tremor['year']
lat_tremor = tremor['latitude']

# Start figure and plot tremor
params = {'legend.fontsize': 24, \
          'xtick.labelsize':24, \
          'ytick.labelsize':24}
pylab.rcParams.update(params)
plt.figure(1, figsize=(20, 10))
plt.scatter(year_tremor, lat_tremor, c='r')

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

all_LFEs = 0
kept_LFEs = 0

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

        all_LFEs = all_LFEs + len(df)

        # Filter LFEs
        df = df.loc[df['cc'] * df['nchannel'] >= threshold['threshold_perm'].iloc[i]]

        kept_LFEs = kept_LFEs + len(df)

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
        time = time[nbLFEs > 3]
        nbLFEs = nbLFEs[nbLFEs > 3]
        plt.scatter(time, np.repeat(latitude, np.shape(time)[0]), s=1 + nbLFEs * 5, c='k', label='Numbers of LFEs')
#        plt.scatter(time, np.repeat(latitude, np.shape(time)[0]), c=nbLFEs, cmap='autumn')

# Plot time line
#tbegin = ymdhms2day(2007, 9, 25, 0, 0, 0)
#tend = ymdhms2day(2009, 5, 14, 0, 0, 0)
#plt.arrow(tbegin, 39.2, tend - tbegin, 0, head_width=0.05, \
#    head_length=0.03, linewidth=4, color='grey', length_includes_head=True)
#plt.arrow(tend, 39.2, tbegin - tend, 0, head_width=0.05, \
#    head_length=0.03, linewidth=4, color='grey', length_includes_head=True)
#plt.annotate('FAME operating', (2008.35, 39.25), fontsize=24, color='grey')

# Families
plt.annotate('A', (xmax + 0.06 * (xmax - xmin), 40.09000), fontsize=24, color='grey')
plt.annotate('B', (xmax + 0.06 * (xmax - xmin), 40.38000), fontsize=24, color='grey')
plt.annotate('C', (xmax + 0.06 * (xmax - xmin), 40.65000), fontsize=24, color='grey')
plt.annotate('D', (xmax + 0.06 * (xmax - xmin), 40.92000), fontsize=24, color='grey')
plt.annotate('E', (xmax + 0.08 * (xmax - xmin), 40.97000), fontsize=24, color='grey')
plt.annotate('F', (xmax + 0.06 * (xmax - xmin), 41.02000), fontsize=24, color='grey')
plt.annotate('G', (xmax + 0.08 * (xmax - xmin), 41.10000), fontsize=24, color='grey')

# End figure
plt.axvline(0.5 * (xmin + xmax), linewidth=2, color='grey')
plt.xlim([xmin, xmax + 0.2 * (xmax - xmin)])
plt.ylim([39.1, latmax])
plt.xlabel('Time (years)', fontsize=24)
plt.ylabel('Latitude', fontsize=24)
plt.title('Tremor and LFEs', fontsize=24)

msizes = [16, 51, 101, 251, 501]
mlabels = ['3', '10', '20', '50', '100']
markers = []
for (size, label) in zip(msizes, mlabels):
   markers.append(plt.scatter([],[], s=size, c='k', label=label))
plt.legend(handles=markers, frameon=False, title='Number of LFEs', loc=4, title_fontsize=24)

plt.tight_layout()
plt.savefig('LFEdistribution_perm/zoom_teq5.eps', format='eps')
plt.close(1)
