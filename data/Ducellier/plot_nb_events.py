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
events.columns = ['group', 'name', 'lon', 'lat', 'x', 'y', 'nb', 'rt']

event_B3 = events.loc[events['name'] == 'B3']
events = events.loc[events['name'] != 'B3']

groups = events['group'].unique()

plt.figure(1, figsize=(8, 8))
cycol = cycle('bgrcmk')

for group in groups:
    events_sub = events.loc[events['group'] == group]
    color = next(cycol)
    if group == 'B':
        color_B = color
    plt.scatter(events_sub['x'], events_sub['rt'], c=color, s=100)
    for i in range(0, len(events_sub)):
        if events_sub['name'].iloc[i] in ['C2', 'C4']:
            plt.annotate(events_sub['name'].iloc[i], \
                (events_sub['x'].iloc[i] - 5, events_sub['rt'].iloc[i]), \
                fontsize=24, color=color)
        elif events_sub['name'].iloc[i] == 'G2':
            plt.annotate(events_sub['name'].iloc[i], \
                (events_sub['x'].iloc[i] + 1, events_sub['rt'].iloc[i] + 10), \
                fontsize=24, color=color)
        elif events_sub['name'].iloc[i] == 'G3':
            plt.annotate(events_sub['name'].iloc[i], \
                (events_sub['x'].iloc[i] + 1, events_sub['rt'].iloc[i] - 10), \
                fontsize=24, color=color)
        else:
            plt.annotate(events_sub['name'].iloc[i], \
                (events_sub['x'].iloc[i] + 1, events_sub['rt'].iloc[i]), \
                fontsize=24, color=color)

plt.arrow(event_B3['x'].iloc[0], 275, 0, 20, width=0.3, color=color_B)
plt.annotate(event_B3['name'].iloc[0] + ' (' + str(int(event_B3['rt'].iloc[0])) + ')', \
    (event_B3['x'].iloc[0] + 1, 275), fontsize=24, color=color_B)

plt.xlabel('Distance in the east direction (km)', fontsize=24)
plt.ylabel('Recurrence time (days)', fontsize=24)
plt.xlim([2, 52])
plt.ylim([0, 300])
plt.tight_layout()
plt.savefig('LFEdistribution_perm/nb_LFE_events.eps', format='eps')
plt.close(1)