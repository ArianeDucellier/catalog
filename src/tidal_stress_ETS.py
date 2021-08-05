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

# Friction
mu = 0.1

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

# Arrival time of 2008 ETS event
begin_time = pd.read_csv('../data/For_Ken_California_LFE/begin_time.txt', sep=' ', header=None)
begin_time.columns = ['latitude', 'longitude', 'begin_time']

# Loop on templates
for i in range(0, np.shape(templates)[0] - 5):

    # Look only at good data
    if threshold['threshold_perm'].iloc[i] > 0.0:

        # Look for beginning time of ETS
        latitude = templates[i][2]
        longitude = templates[i][3]
        subset = begin_time.loc[(begin_time['latitude'] == latitude) & \
            (begin_time['longitude'] == longitude)]
        time_ETS = subset.iloc[0]['begin_time']

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

        # Divide catalog between first and last part of ETS event
        t2 = ymdhms2matlab(2008, 4, 1, 0, 0, 0) + time_ETS
        t1 = t2 - 1.5
        t3 = t2 + 10
        df1 = df.loc[(df['date'] >= t1) & (df['date'] <= t2)].reset_index()
        df2 = df.loc[(df['date'] >= t2) & (df['date'] <= t3)].reset_index()
        
        # Look for value of tidal stress at the time of each LFE
        if len(df1) > 0:
            shear_LFE1 = np.repeat(np.nan, len(df1))
            norm_LFE1 = np.repeat(np.nan, len(df1))
            vol_LFE1 = np.repeat(np.nan, len(df1))
            for k in range(len(df1)):
                indices = np.where(24 * np.abs(time - df1.date.loc[k]) <= 1)[1]
                if len(indices) > 0:
                    shear_LFE1[k] = np.interp(df1.date.loc[k], time[0, indices], shear[indices])
                    norm_LFE1[k] = np.interp(df1.date.loc[k], time[0, indices], norm[indices])
                    vol_LFE1[k] = np.interp(df1.date.loc[k], time[0, indices], vol[indices])
            df1['shear'] = shear_LFE1
            df1['norm'] = norm_LFE1
            df1['vol'] = vol_LFE1
            df1['coulomb'] = shear_LFE1 + mu * (norm_LFE1 + vol_LFE1)

        # Look for value of tidal stress at the time of each LFE
        if len(df2) > 0:
            shear_LFE2 = np.repeat(np.nan, len(df2))
            norm_LFE2 = np.repeat(np.nan, len(df2))
            vol_LFE2 = np.repeat(np.nan, len(df2))
            for k in range(len(df2)):
                indices = np.where(24 * np.abs(time - df2.date.loc[k]) <= 1)[1]
                if len(indices) > 0:
                    shear_LFE2[k] = np.interp(df2.date.loc[k], time[0, indices], shear[indices])
                    norm_LFE2[k] = np.interp(df2.date.loc[k], time[0, indices], norm[indices])
                    vol_LFE2[k] = np.interp(df2.date.loc[k], time[0, indices], vol[indices])
            df2['shear'] = shear_LFE2
            df2['norm'] = norm_LFE2
            df2['vol'] = vol_LFE2
            df2['coulomb'] = shear_LFE2 + mu * (norm_LFE2 + vol_LFE2)

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
            if len(df1) > 0:
                count[k] = len(df1.loc[(df1.shear >= boxes[k]) & (df1.shear <= boxes[k + 1])])
            count_shear[k] = np.sum((shear >= boxes[k]) & (shear <= boxes[k + 1]))

        ax1 = plt.subplot(321)
        expected = len(df1) * count_shear / np.sum(count_shear)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Shear stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)
        plt.title('Start of ETS event', fontsize=24)

        count = np.zeros(nbox)
        count_shear = np.zeros(nbox)
        for k in range(0, nbox):
            if len(df2) > 0:
                count[k] = len(df2.loc[(df2.shear >= boxes[k]) & (df2.shear <= boxes[k + 1])])
            count_shear[k] = np.sum((shear >= boxes[k]) & (shear <= boxes[k + 1]))

        ax2 = plt.subplot(322)
        expected = len(df2) * count_shear / np.sum(count_shear)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Shear stress (kPa)', fontsize=24)
        plt.title('End of ETS event', fontsize=24)

        # Normal
        min_stress = np.min(norm)
        max_stress = np.max(norm)
        boxes = np.linspace(min_stress, max_stress, nbox + 1)
        size_box = (max_stress - min_stress) / (4 * nbox)

        count = np.zeros(nbox)
        count_norm = np.zeros(nbox)
        for k in range(0, nbox):
            if len(df1) > 0:
                count[k] = len(df1.loc[(df1.norm >= boxes[k]) & (df1.norm <= boxes[k + 1])])
            count_norm[k] = np.sum((norm >= boxes[k]) & (norm <= boxes[k + 1]))

        ax3 = plt.subplot(323)
        expected = len(df1) * count_norm / np.sum(count_norm)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Normal stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        count = np.zeros(nbox)
        count_norm = np.zeros(nbox)
        for k in range(0, nbox):
            if len(df2) > 0:
                count[k] = len(df2.loc[(df2.norm >= boxes[k]) & (df2.norm <= boxes[k + 1])])
            count_norm[k] = np.sum((norm >= boxes[k]) & (norm <= boxes[k + 1]))

        ax4 = plt.subplot(324)
        expected = len(df2) * count_norm / np.sum(count_norm)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Normal stress (kPa)', fontsize=24)

        # Volume
        min_stress = np.min(vol)
        max_stress = np.max(vol)
        boxes = np.linspace(min_stress, max_stress, nbox + 1)
        size_box = (max_stress - min_stress) / (4 * nbox)

        count = np.zeros(nbox)
        count_vol = np.zeros(nbox)
        for k in range(0, nbox):
            if len(df1) > 0:
                count[k] = len(df1.loc[(df1.vol >= boxes[k]) & (df1.vol <= boxes[k + 1])])
            count_vol[k] = np.sum((vol >= boxes[k]) & (vol <= boxes[k + 1]))

        ax5 = plt.subplot(325)
        expected = len(df1) * count_vol / np.sum(count_vol)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Volume stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        count = np.zeros(nbox)
        count_vol = np.zeros(nbox)
        for k in range(0, nbox):
            if len(df2) > 0:
                count[k] = len(df2.loc[(df2.vol >= boxes[k]) & (df2.vol <= boxes[k + 1])])
            count_vol[k] = np.sum((vol >= boxes[k]) & (vol <= boxes[k + 1]))

        ax6 = plt.subplot(326)
        expected = len(df2) * count_vol / np.sum(count_vol)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Volume stress (kPa)', fontsize=24)

        plt.suptitle('Influence of tidal stress on LFE occurrence', fontsize=30)
        plt.savefig('tidal_stress_ETS/' + templates[i][0].astype(str) + '.eps', format='eps')
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()
        ax6.clear()
        plt.close(1)

        # Plot Coulomb stress
        plt.figure(1, figsize=(20, 6))

        coulomb = shear + mu * (norm + vol)
        min_stress = np.min(coulomb)
        max_stress = np.max(coulomb)
        boxes = np.linspace(min_stress, max_stress, nbox + 1)
        size_box = (max_stress - min_stress) / (4 * nbox)

        count = np.zeros(nbox)
        count_coulomb = np.zeros(nbox)
        for k in range(0, nbox):
            if len(df1) > 0:
                count[k] = len(df1.loc[(df1.coulomb >= boxes[k]) & (df1.coulomb <= boxes[k + 1])])
            count_coulomb[k] = np.sum((coulomb >= boxes[k]) & (coulomb <= boxes[k + 1]))

        ax1 = plt.subplot(121)
        expected = len(df1) * count_coulomb / np.sum(count_coulomb)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Coulomb stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)
        plt.title('Start of ETS event', fontsize=24)

        count = np.zeros(nbox)
        count_coulomb = np.zeros(nbox)
        for k in range(0, nbox):
            if len(df2) > 0:
                count[k] = len(df2.loc[(df2.coulomb >= boxes[k]) & (df2.coulomb <= boxes[k + 1])])
            count_coulomb[k] = np.sum((coulomb >= boxes[k]) & (coulomb <= boxes[k + 1]))

        ax2 = plt.subplot(122)
        expected = len(df2) * count_coulomb / np.sum(count_coulomb)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
        plt.xlabel('Coulomb stress (kPa)', fontsize=24)
        plt.title('End of ETS event', fontsize=24)

        plt.suptitle('Influence of tidal stress on LFE occurrence', fontsize=30)
        plt.tight_layout()
        plt.savefig('tidal_stress_ETS/' + templates[i][0].astype(str) + '_coulomb.eps', format='eps')
        ax1.clear()
        ax2.clear()
        plt.close(1)
