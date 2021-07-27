"""
Script to plot cumulative number of LFEs
"""

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

from datetime import datetime, timedelta

# List of LFE families
templates = np.loadtxt('../Plourde_2015/templates_list.txt', \
    dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
    'eZ', 'nb'), \
         'formats': ('S13', 'S3', np.float, np.float, np.float, \
    np.float, np.float, np.int)}, \
    skiprows=1)

tbegin = datetime(2007, 7, 1)
nb_days = 731

# Threshold for filtering the catalog
threshold = pd.read_csv('threshold_cc.txt', sep=r'\s{1,}', header=None, engine='python')
threshold.columns = ['family', 'threshold_FAME', 'threshold_perm']

# Loop on templates
for i in range(0, np.shape(templates)[0]):

    # Open LFE catalog (FAME data)
    namedir = 'catalogs/' + templates[i][0].astype(str)
    namefile = namedir + '/catalog_2007_2009.pkl'
    df_FAME = pickle.load(open(namefile, 'rb'))

    # Open LFE catalog (permanent networks)
    namefile = namedir + '/catalog_2004_2011.pkl'
    df_perm = pickle.load(open(namefile, 'rb'))

    # Filter LFEs
    df_FAME = df_FAME.loc[df_FAME['cc'] * df_FAME['nchannel'] >= threshold['threshold_FAME'].iloc[i]]
    df_perm = df_perm.loc[df_perm['cc'] * df_perm['nchannel'] >= threshold['threshold_perm'].iloc[i]]

    # Add date column
    df = pd.DataFrame({'year': df_FAME['year'], 'month': df_FAME['month'], 'day': df_FAME['day']})
    df = pd.to_datetime(df)
    df_FAME['date'] = df

    df = pd.DataFrame({'year': df_perm['year'], 'month': df_perm['month'], 'day': df_perm['day']})
    df = pd.to_datetime(df)
    df_perm['date'] = df

    # Cumulative number of LFEs
    cum_FAME = np.zeros(nb_days)
    cum_perm = np.zeros(nb_days)
    for j in range(0, nb_days):
        df_sub = df_FAME[(df_FAME['date'] >= tbegin) & (df_FAME['date'] < tbegin + timedelta(j + 1))]
        cum_FAME[j] = len(df_sub)
        df_sub = df_perm[(df_perm['date'] >= tbegin) & (df_perm['date'] < tbegin + timedelta(j + 1))]
        cum_perm[j] = len(df_sub)

    # Normalize
    nb_FAME = int(np.max(cum_FAME))
    nb_perm = int(np.max(cum_perm))
    cum_FAME = cum_FAME / nb_FAME
    cum_perm = cum_perm / nb_perm

    # Plot
    params = {'legend.fontsize': 24, \
              'xtick.labelsize':24, \
              'ytick.labelsize':24}
    pylab.rcParams.update(params)
    plt.figure(1, figsize=(10, 10))
    plt.plot(np.arange(0, nb_days), cum_FAME, 'r-', label='FAME ({:d} LFEs)'.format(nb_FAME))
    plt.plot(np.arange(0, nb_days), cum_perm, 'b-', label='networks ({:d} LFEs)'.format(nb_perm))
    plt.yticks([], [])
    plt.xlabel('Time (days) since {:02d}/{:02d}/{:04d}'. \
        format(tbegin.month, tbegin.day, tbegin.year), fontsize=24)
    plt.ylabel('Normalized number of LFEs', fontsize=24)
    plt.title('Cumulative number of LFEs', fontsize=24)
    plt.legend(loc=4, fontsize=20)
    plt.savefig('cumulative/' + templates[i][0].astype(str) + '.eps', format='eps')
    plt.close(1)
