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
B = 0.5

# LFE families
#families = [0, 0, 5, 14, 20, 28, 40]
#names = ['A (all LFEs)', 'A', 'B2', 'D', 'C4', 'F', 'G1']
#families = [0, 5, 14, 20, 28, 40]
#names = ['A', 'B2', 'D', 'C4', 'F', 'G1']
families = [0, 0]
names = ['A (2008 ETS event)', 'A (all LFEs)']

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

#plt.figure(1, figsize=(28, 10))
#plt.figure(1, figsize=(14, 15))
plt.figure(1, figsize=(14, 5))

# Loop on templates
for i in range(0, len(families)):

    # Look for beginning time of ETS
    latitude = templates[families[i]][2]
    longitude = templates[families[i]][3]
    subset = begin_time.loc[(begin_time['latitude'] == latitude) & \
        (begin_time['longitude'] == longitude)]
    time_ETS = subset.iloc[0]['begin_time']

    # Load data
    filename = '../data/tidal_stress/' + templates[families[i]][0].astype(str) + '.pkl'
    tidal_stress = pickle.load(open(filename, 'rb'))
    time = tidal_stress[0]
    norm = tidal_stress[1]
    shear = tidal_stress[2]
    vol = tidal_stress[3]

    # Open LFE catalog
    namedir = '../data/Ducellier/catalogs/' + templates[families[i]][0].astype(str)
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
#    if i == 0:
    if i == 1:
        df2 = df.copy(deep=True)
    else:
        t2 = ymdhms2matlab(2008, 4, 1, 0, 0, 0) + time_ETS
        t3 = t2 + 10
        df2 = df.loc[(df['date'] >= t2) & (df['date'] <= t3)].reset_index()
        
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
        df2['coulomb'] = shear_LFE2 + mu * (norm_LFE2 - B * vol_LFE2)

        # Plot Coulomb stress
        coulomb = shear + mu * (norm - B * vol)
        min_stress = np.min(coulomb)
        max_stress = np.max(coulomb)
        boxes = np.linspace(min_stress, max_stress, nbox + 1)
        size_box = (max_stress - min_stress) / (4 * nbox)

        count = np.zeros(nbox)
        count_coulomb = np.zeros(nbox)
        for k in range(0, nbox):
            if len(df2) > 0:
                count[k] = len(df2.loc[(df2.coulomb >= boxes[k]) & (df2.coulomb <= boxes[k + 1])])
            count_coulomb[k] = np.sum((coulomb >= boxes[k]) & (coulomb <= boxes[k + 1]))

#        if i == 0:
#            plt.subplot2grid((2, 4), (0, 3))
#        elif i == 1:
#            plt.subplot2grid((2, 4), (0, 0))
#        elif i == 2:
#            plt.subplot2grid((2, 4), (0, 1))
#        elif i == 3:
#            plt.subplot2grid((2, 4), (0, 2))
#        elif i == 4:
#            plt.subplot2grid((2, 4), (1, 0))
#        elif i == 5:
#            plt.subplot2grid((2, 4), (1, 1))
#        else:
#            plt.subplot2grid((2, 4), (1, 2))

#        if i == 0:
#            plt.subplot2grid((3, 2), (0, 0))
#        elif i == 1:
#            plt.subplot2grid((3, 2), (0, 1))
#        elif i == 2:
#            plt.subplot2grid((3, 2), (1, 0))
#        elif i == 3:
#            plt.subplot2grid((3, 2), (1, 1))
#        elif i == 4:
#            plt.subplot2grid((3, 2), (2, 0))
#        else:
#            plt.subplot2grid((3, 2), (2, 1))

        if i == 0:
            plt.subplot2grid((1, 2), (0, 0))
        else:
            plt.subplot2grid((1, 2), (0, 1))

        expected = len(df2) * count_coulomb / np.sum(count_coulomb)
        plt.plot(0.5 * (boxes[1:] + boxes[:-1]), expected, color='black')
        plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, size_box)
#        if (i == 0) or (i >= 4):
#            plt.xlabel('Coulomb stress (kPa)', fontsize=24)
#        if (i == 1) or (i == 4):
#            plt.ylabel('Number of LFEs', fontsize=24)
#        if (i == 4) or (i == 5):
        plt.xlabel('Coulomb stress (kPa)', fontsize=24)
        if (i == 0) or (i == 2) or (i == 4):
            plt.ylabel('Number of LFEs', fontsize=24)
        plt.title('Family ' + names[i], fontsize=24)

plt.tight_layout()
plt.savefig('tidal_stress_ETS/coulomb_A.eps', format='eps')
plt.close(1)
