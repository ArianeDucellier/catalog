import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

from math import floor
from scipy.io import loadmat

from date import ymdhms2matlab, matlab2ymdhms

params = {'legend.fontsize': 24, \
          'xtick.labelsize':24, \
          'ytick.labelsize':24}
pylab.rcParams.update(params)

# Number of boxes for plotting
nbox = 20

# List of LFE families
templates = np.loadtxt('../data/Plourde_2015/templates_list.txt', \
    dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
    'eZ', 'nb'), \
         'formats': ('S13', 'S3', np.float, np.float, np.float, \
    np.float, np.float, np.int)}, \
    skiprows=1)

# Threshold for filtering the catalog
threshold = pd.read_csv('../data/Ducellier/threshold_cc.txt', sep=r'\s{1,}', \
    header=None, engine='python')
threshold.columns = ['family', 'threshold_FAME', 'threshold_perm']

# Loop on templates
for i in range(0, np.shape(templates)[0]):

    # Look only at good data
    if threshold['threshold_perm'].iloc[i] > 0.0:

        # Load data
        filename = '../data/tidal_stress/' + templates[i][0].astype(str) + '.pkl'
        tidal_stress = pickle.load(open(filename, 'rb'))
        time = tidal_stress[0]
        norm = tidal_stress[1]
        shear = tidal_stress[2]
        vol = tidal_stress[3]

        # Open LFE catalog
        namedir = '../data/Ducellier/catalogs/' + templates[i][0].astype(str)
        filename = namedir + '/catalog_2004_2011.pkl'
        df = pickle.load(open(filename, 'rb'))
        df = df[['year', 'month', 'day', 'hour', 'minute', 'second', \
           'cc', 'nchannel']]
        df = df.astype({'year': int, 'month': int, 'day': int, 'hour': int, \
            'minute': int, 'second': float, 'cc': float, 'nchannel': int})
        df = df.loc[df['cc'] * df['nchannel'] >= \
             threshold['threshold_perm'].iloc[i]]
        df.reset_index(inplace=True)

        # Convert LFE catalog to date in matlab format
        date = np.zeros(len(df))
        for k in range(len(df)):
            date[k] = ymdhms2matlab(df.year.loc[k], df.month.loc[k], \
                df.day.loc[k], df.hour.loc[k], df.minute.loc[k], \
                int(floor(df.second.loc[k])))
        df['date'] = date

        # Look for value of tidal stress at the time of each LFE 
        shear_LFE = np.repeat(np.nan, len(df))
        norm_LFE = np.repeat(np.nan, len(df))
        vol_LFE = np.repeat(np.nan, len(df))
        for k in range(len(df)):
            indices = np.where(24 * np.abs(time - df.date.loc[k]) <= 1)[1]
            if len(indices) > 0:
                shear_LFE[k] = np.interp(df.date.loc[k], time[0, indices], shear[indices])
                norm_LFE[k] = np.interp(df.date.loc[k], time[0, indices], norm[indices])
                vol_LFE[k] = np.interp(df.date.loc[k], time[0, indices], vol[indices])
        df['shear'] = shear_LFE
        df['norm'] = norm_LFE
        df['vol'] = vol_LFE

        # Plot
        plt.figure(1, figsize=(20, 18))

        # Shear
        min_stress = np.min(shear)
        max_stress = np.max(shear)
        boxes = np.linspace(min_stress, max_stress, nbox + 1)
        size_box = (max_stress - min_stress) / (4 * nbox)

        count = np.zeros(nbox)
        count_shear = np.zeros(nbox)
        for k in range(0, nbox):
            count[k] = len(df.loc[(df.shear >= boxes[k]) & (df.shear <= boxes[k + 1])])
            count_shear[k] = np.sum((shear >= boxes[k]) & (shear <= boxes[k + 1]))
        normalize = np.sum(count) / np.sum(count / count_shear)

        ax1 = plt.subplot(321)
        expected = len(df) * count_shear / np.sum(count_shear)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Shear stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        ax2 = plt.subplot(322)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), normalize * expected / count_shear, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), normalize * count / count_shear, size_box)
        plt.xlabel('Shear stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        # Normal
        min_stress = np.min(norm)
        max_stress = np.max(norm)
        boxes = np.linspace(min_stress, max_stress, nbox + 1)
        size_box = (max_stress - min_stress) / (4 * nbox)

        count = np.zeros(nbox)
        count_norm = np.zeros(nbox)
        for k in range(0, nbox):
            count[k] = len(df.loc[(df.norm >= boxes[k]) & (df.norm <= boxes[k + 1])])
            count_norm[k] = np.sum((norm >= boxes[k]) & (norm <= boxes[k + 1]))
        normalize = np.sum(count) / np.sum(count / count_norm)

        ax3 = plt.subplot(323)
        expected = len(df) * count_norm / np.sum(count_norm)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Normal stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        ax4 = plt.subplot(324)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), normalize * expected / count_norm, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), normalize * count / count_norm, size_box)
        plt.xlabel('Normal stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        # Volume
        min_stress = np.min(vol)
        max_stress = np.max(vol)
        boxes = np.linspace(min_stress, max_stress, nbox + 1)
        size_box = (max_stress - min_stress) / (4 * nbox)

        count = np.zeros(nbox)
        count_vol = np.zeros(nbox)
        for k in range(0, nbox):
            count[k] = len(df.loc[(df.vol >= boxes[k]) & (df.vol <= boxes[k + 1])])
            count_vol[k] = np.sum((vol >= boxes[k]) & (vol <= boxes[k + 1]))
        normalize = np.sum(count) / np.sum(count / count_vol)

        ax5 = plt.subplot(325)
        expected = len(df) * count_vol / np.sum(count_vol)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Volume stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        ax6 = plt.subplot(326)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), normalize * expected / count_vol, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), normalize * count / count_vol, size_box)
        plt.xlabel('Volume stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        plt.suptitle('Influence of tidal stress on LFE occurrence', fontsize=30)
        plt.savefig('tidal_stress/' + templates[i][0].astype(str) + '.eps', format='eps')
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()
        ax6.clear()
        plt.close(1)
