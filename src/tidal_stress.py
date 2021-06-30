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

# Tidal stress file
data = loadmat('../data/For_Ken_California_LFE/LFEs2.mat')

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

        # Look for corresponding location of tidal stress
        latitude = templates[i][2]
        longitude = templates[i][3]
        for k in range(0, len(data['LFE']['TempLabel'][0][0])):
            if ((data['LFE']['lat'][0][0][k] == latitude) & \
                (data['LFE']['lon'][0][0][k] == longitude)):
                namefile = data['LFE']['TempLabel'][0][0][k][0][0]
                j = k
        print(templates[i][0].astype(str), namefile, i, j)

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
        shear = np.repeat(np.nan, len(df))
        norm = np.repeat(np.nan, len(df))
        vol = np.repeat(np.nan, len(df))
        for k in range(len(df)):
            indices = np.where(24 * 4 * np.abs(data['LFE']['tvec'][0][0][0, j] - df.date.loc[k]) <= 1)[1]
            if len(indices) > 0:
                shear[k] = np.mean(data['LFE']['shear'][0][0][0, j][0, indices])
                norm[k] = np.mean(data['LFE']['norm'][0][0][0, j][0, indices])
                vol[k] = np.mean(data['LFE']['vol'][0][0][0, j][0, indices])
        df['shear'] = shear
        df['norm'] = norm
        df['vol'] = vol

        # Plot
        plt.figure(1, figsize=(10, 18))

        # Shear
        min_stress = np.min(data['LFE']['shear'][0][0][0, j])
        max_stress = np.max(data['LFE']['shear'][0][0][0, j])
        boxes = np.linspace(min_stress, max_stress, nbox + 1)

        count = np.zeros(nbox)
        count_shear = np.zeros(nbox)
        for k in range(0, nbox):
            count[k] = len(df.loc[(df.shear >= boxes[k]) & (df.shear <= boxes[k + 1])])
            count_shear[k] = np.sum((data['LFE']['shear'][0][0][0, j] >= boxes[k]) & (data['LFE']['shear'][0][0][0, j] <= boxes[k + 1]))
        print(count)
        print(count_shear)
        norm = np.sum(count) / np.sum(count / count_shear)

        ax1 = plt.subplot(311)
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), norm * count / count_shear, 0.05)
        plt.xlabel('Shear stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        # Normal
        min_stress = np.min(data['LFE']['norm'][0][0][0, j])
        max_stress = np.max(data['LFE']['norm'][0][0][0, j])
        boxes = np.linspace(min_stress, max_stress, nbox + 1)

        count = np.zeros(nbox)
        count_norm = np.zeros(nbox)
        for k in range(0, nbox):
            count[k] = len(df.loc[(df.norm >= boxes[k]) & (df.norm <= boxes[k + 1])])
            count_norm[k] = np.sum((data['LFE']['norm'][0][0][0, j] >= boxes[k]) & (data['LFE']['norm'][0][0][0, j] <= boxes[k + 1]))
        norm = np.sum(count) / np.sum(count / count_norm)

        ax2 = plt.subplot(312)
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), norm * count / count_norm, 0.05)
        plt.xlabel('Normal stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        # Volume
        min_stress = np.min(data['LFE']['vol'][0][0][0, j])
        max_stress = np.max(data['LFE']['vol'][0][0][0, j])
        boxes = np.linspace(min_stress, max_stress, nbox + 1)

        count = np.zeros(nbox)
        count_vol = np.zeros(nbox)
        for k in range(0, nbox):
            count[k] = len(df.loc[(df.vol >= boxes[k]) & (df.vol <= boxes[k + 1])])
            count_vol[k] = np.sum((data['LFE']['vol'][0][0][0, j] >= boxes[k]) & (data['LFE']['vol'][0][0][0, j] <= boxes[k + 1]))
        norm = np.sum(count) / np.sum(count / count_vol)

        ax3 = plt.subplot(313)
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), norm * count / count_vol, 0.05)
        plt.xlabel('Volume stress (kPa)', fontsize=24)
        plt.ylabel('Number of LFEs', fontsize=24)

        plt.suptitle('Influence of tidal stress on LFE occurrence', fontsize=30)
        plt.savefig('tidal_stress/' + templates[i][0].astype(str) + '.eps', format='eps')
        ax1.clear()
        ax2.clear()
        ax3.clear()
        plt.close(1)
