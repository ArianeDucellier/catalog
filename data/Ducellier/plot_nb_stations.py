"""
Script to read the number of stations used to detect LFEs
and plot it as a function of time
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# List of LFE families
templates = np.loadtxt('../Plourde_2015/templates_list.txt', \
    dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
    'eZ', 'nb'), \
         'formats': ('S13', 'S3', np.float, np.float, np.float, \
    np.float, np.float, np.int)}, \
    skiprows=1)

for i in range(0, np.shape(templates)[0]):
    filename = 'nstations/' + templates[i][0].astype(str) + '.txt'
    df = pd.read_csv(filename, sep=r'\s{1,}', header=None, \
        engine='python')
    df.columns = ['year', 'month', 'day', 'hour', 'nchannel']
    date = pd.to_datetime(df.drop(columns='nchannel'))
    df['date'] = date
    plt.figure(1, figsize=(20, 10))
    plt.plot(df['date'], df['nchannel'], 'bo')
    plt.xlabel('Time', fontsize=24)
    plt.ylabel('Number of cgannels', fontsize=24)
    plt.title('Family {}'.format(templates[i][0].astype(str)), \
        fontsize=24)
    plt.savefig('nstations/' + templates[i][0].astype(str) + '.eps', format='eps')
    plt.close(1)
