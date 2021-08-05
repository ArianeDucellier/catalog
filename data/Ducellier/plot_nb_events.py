"""
This script plots the number of events as a function
of distance from the up-dip limit of tremor
"""

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from itertools import cycle

params = {'legend.fontsize': 24, \
          'xtick.labelsize':24, \
          'ytick.labelsize':24}
pylab.rcParams.update(params)

events = pd.read_csv('LFEdistribution_perm/nb_LFE_events.txt', sep=r'\s{1,}', header=None, engine='python')
events.columns = ['group', 'name', 'lon', 'lat', 'x', 'y', 'nb']

groups = events['group'].unique()

plt.figure(1, figsize=(8, 8))
cycol = cycle('bgrcmk')

for group in groups:
    events_sub = events.loc[events['group'] == group]
    color = next(cycol)
    if group == 'A':
        plt.arrow(events_sub['x'].iloc[0], 50, 0, 2, width=0.3, color=color)
        plt.annotate(events_sub['name'].iloc[0] + ' (' + str(int(events_sub['nb'].iloc[0])) + ')', \
            (events_sub['x'].iloc[0] + 1, 50), fontsize=24, color=color)
    else:
        plt.scatter(events_sub['x'], events_sub['nb'], c=color, s=100)
        for i in range(0, len(events_sub)):
            if events_sub['name'].iloc[i] == 'G2':
                plt.annotate(events_sub['name'].iloc[i], \
                    (events_sub['x'].iloc[i] - 5, events_sub['nb'].iloc[i]), \
                    fontsize=24, color=color)
            else:
                plt.annotate(events_sub['name'].iloc[i], \
                    (events_sub['x'].iloc[i] + 1, events_sub['nb'].iloc[i]), \
                    fontsize=24, color=color)

plt.xlabel('Distance in the east direction (km)', fontsize=24)
plt.ylabel('Number of LFE events', fontsize=24)
plt.xlim([2, 52])
plt.ylim([5, 54])
plt.savefig('LFEdistribution_perm/nb_LFE_events.eps', format='eps')
plt.close(1)