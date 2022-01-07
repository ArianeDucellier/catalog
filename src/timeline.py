import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from datetime import datetime, timedelta

station_file = '../data/Ducellier/stations_temporary.txt'

tmin = datetime(2007, 7, 1)
tmax = datetime(2009, 7, 1)
ticks_time = [datetime(2007, 7, 1), datetime(2008, 1, 1), \
    datetime(2008, 7, 1), datetime(2009, 1, 1), datetime(2009, 7, 1)]
ticks_labels = ['2007-07-01', '2008-01-01', '2008-07-01',\
    '2009-01-01', '2009-07-01']

#tmin = datetime(2004, 1, 1)
#tmax = datetime(2012, 1, 1)
#ticks_time = [datetime(2004, 1, 1), datetime(2006, 1, 1), \
#    datetime(2008, 1, 1), datetime(2010, 1, 1), datetime(2012, 1, 1)]
#ticks_labels = ['2004-01-01', '2006-01-01', '2008-01-01',\
#    '2010-01-01', '2012-01-01']

staloc = pd.read_csv(station_file, sep=r'\s{1,}', header=None, \
    engine='python')
staloc.columns = ['station', 'network', 'channels', 'location', 'server', \
    'latitude', 'longitude', 'time_on', 'time_off']

staloc['time_on']= pd.to_datetime(staloc['time_on'])
staloc['time_off']= pd.to_datetime(staloc['time_off'])

staloc['channels'] = staloc.apply(lambda x: len(x['channels'].split(',')), axis=1)
staloc['channels'] = staloc['channels'].astype(int)

staloc = staloc.groupby(['station', 'network', 'channels', 'server', \
    'latitude', 'longitude']).agg({ \
    'time_on':lambda x: min(x), \
    'time_off':lambda x: max(x)}).reset_index()

t = np.arange(tmin, tmax, timedelta(days=1)).astype(datetime)
nstations = np.zeros(len(t))

for i in range(0, len(t)):
    df = staloc.loc[(staloc['time_on'] <= t[i]) & (staloc['time_off'] >= t[i])]
    nstations[i] = int(df['channels'].sum())

# Figure
params = {'legend.fontsize': 20, \
          'xtick.labelsize':20, \
          'ytick.labelsize':20}
pylab.rcParams.update(params)
plt.figure(1, figsize=(10, 6))
ax = plt.axes()
plt.plot(t, nstations, 'b-')
plt.xlabel('Time', fontsize=24)
plt.ylabel('Number of stations', fontsize=24)
plt.ylim([0, np.max(nstations) + 1])
ax.set_xticks(ticks_time)
ax.set_xticklabels(ticks_labels)
plt.title('Evolution of the number of channels', fontsize=24)
plt.tight_layout()
plt.savefig('timeline_FAME.eps', format='eps')
plt.close(1)
