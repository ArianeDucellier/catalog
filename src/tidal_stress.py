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

i = 0

nbox = 20

data = loadmat('../data/tidal_stress/LFEs2.mat')

namefile = data['LFE']['TempLabel'][0][0][i][0][0]

threshold = pd.read_csv('../data/Ducellier/threshold_cc.txt', sep=r'\s{1,}', header=None, engine='python')
threshold.columns = ['family', 'threshold_FAME', 'threshold_perm']

namedir = '../data/Ducellier/catalogs/' + data['LFE']['TempLabel'][0][0][0][0][0]
namefile = namedir + '/catalog_2004_2011.pkl'
df = pickle.load(open(namefile, 'rb'))
df = df[['year', 'month', 'day', 'hour', 'minute', 'second', 'cc', 'nchannel']]
df = df.astype({'year': int, 'month': int, 'day': int, 'hour': int, 'minute': int, 'second': float, 'cc': float, 'nchannel': int})
df = df.loc[df['cc'] * df['nchannel'] >= threshold['threshold_perm'].iloc[i]]
df.reset_index(inplace=True)

date = np.zeros(len(df))
for j in range(len(df)):
    date[j] = ymdhms2matlab(df.year.loc[j], df.month.loc[j], df.day.loc[j], df.hour.loc[j], df.minute.loc[j], int(floor(df.second.loc[j])))
df['date'] = date

shear = np.zeros(len(df))
norm = np.zeros(len(df))
vol = np.zeros(len(df))
for j in range(len(df)):
    indices = np.where(24 * 4 * np.abs(data['LFE']['tvec'][0][0][0, i] - df.date.loc[j]) <= 1)[1]
    if len(indices) > 0:
        shear[j] = np.mean(data['LFE']['shear'][0][0][0, i][0, indices])
        norm[j] = np.mean(data['LFE']['norm'][0][0][0, i][0, indices])
        vol[j] = np.mean(data['LFE']['vol'][0][0][0, i][0, indices])
df['shear'] = shear
df['norm'] = norm
df['vol'] = vol

# Shear
min_stress = np.min(data['LFE']['shear'][0][0][0, i])
max_stress = np.max(data['LFE']['shear'][0][0][0, i])
boxes = np.linspace(min_stress, max_stress, nbox + 1)

count = np.zeros(nbox)
for j in range(0, nbox):
    count[j] = len(df.loc[(df.shear >= boxes[j]) & (df.shear <= boxes[j + 1])])

plt.figure(1, figsize=(10, 6))
plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, 0.05)
plt.xlabel('Shear stress (kPa)', fontsize=24)
plt.ylabel('Number of LFEs', fontsize=24)
plt.title('Influence of tidal stress on LFE occurrence', fontsize=24)
plt.tight_layout()
plt.savefig('tidal_stress/shear.eps', format='eps')
plt.close(1)

# Normal
min_stress = np.min(data['LFE']['norm'][0][0][0, i])
max_stress = np.max(data['LFE']['norm'][0][0][0, i])
boxes = np.linspace(min_stress, max_stress, nbox + 1)

count = np.zeros(nbox)
for j in range(0, nbox):
    count[j] = len(df.loc[(df.norm >= boxes[j]) & (df.norm <= boxes[j + 1])])

plt.figure(1, figsize=(10, 6))
plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, 0.05)
plt.xlabel('Normal stress (kPa)', fontsize=24)
plt.ylabel('Number of LFEs', fontsize=24)
plt.title('Influence of tidal stress on LFE occurrence', fontsize=24)
plt.tight_layout()
plt.savefig('tidal_stress/normal.eps', format='eps')
plt.close(1)

# Volume
min_stress = np.min(data['LFE']['vol'][0][0][0, i])
max_stress = np.max(data['LFE']['vol'][0][0][0, i])
boxes = np.linspace(min_stress, max_stress, nbox + 1)

count = np.zeros(nbox)
for j in range(0, nbox):
    count[j] = len(df.loc[(df.vol >= boxes[j]) & (df.vol <= boxes[j + 1])])

plt.figure(1, figsize=(10, 6))
plt.bar(0.5 * (boxes[1:] + boxes[:-1]), count, 0.05)
plt.xlabel('Volume stress (kPa)', fontsize=24)
plt.ylabel('Number of LFEs', fontsize=24)
plt.title('Influence of tidal stress on LFE occurrence', fontsize=24)
plt.tight_layout()
plt.savefig('tidal_stress/vol.eps', format='eps')
plt.close(1)