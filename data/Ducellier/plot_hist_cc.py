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

# Loop on templates
for i in range(0, np.shape(templates)[0]):

    # Open LFE catalog
    namedir = 'catalogs/' + templates[i][0].astype(str)
    namefile = namedir + '/catalog_2002_2011.pkl'
    df = pickle.load(open(namefile, 'rb'))
    time = pd.DataFrame({'year':df['year'], 'month':df['month'], 'day':df['day']})
    time = pd.to_datetime(time)
    df['time'] = time

    print(templates[i][0].astype(str), np.quantile(df['cc'] * df['nchannel'], 0.8))

    # Plot figure
    plt.figure(1, figsize=(15, 5))
    ax1 = plt.subplot(131)
    plt.hist(df['cc'])
    plt.xlabel('Cross-correlation value', fontsize=24)
    plt.ylabel('Number of LFEs', fontsize=24)
    plt.title('Cross-correlation', fontsize=30)
    ax2 = plt.subplot(132)
    plt.hist(df['cc'] * df['nchannel'])
    plt.xlabel('CC x nb channels', fontsize=24)
    plt.ylabel('Number of LFEs', fontsize=24)
    plt.title('CC x nb channels', fontsize=30)
    ax3 = plt.subplot(133)
    plt.plot(df['time'], df['nchannel'], 'ko')
    plt.xlabel('Time', fontsize=24)
    plt.ylabel('Number of channels', fontsize=24)
    plt.title('Number of channels', fontsize=30)
    plt.savefig('histograms/' + templates[i][0].astype(str) + '_cc.eps', format='eps')
    ax1.clear()
    ax2.clear()
    ax3.clear()
    plt.close(1)

