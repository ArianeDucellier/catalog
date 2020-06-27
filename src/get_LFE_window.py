import numpy as np
import pandas as pd
import pickle

from math import ceil, cos, pi, sin, sqrt
from sklearn import linear_model

def pick_arrivals(time, data, t0, threshold):
    """
    """
    N = int(ceil((time[-1] - time[0]) / t0))
    maxima = []
    tmaxima = []
    RMS = np.sqrt(np.mean(np.square(data)))
    for i in range(0, N):
        subtime = time[(time >= i * t0) & (time < (i + 1) * t0)]
        subdata = data[(time >= i * t0) & (time < (i + 1) * t0)]
        imax = np.argmax(subdata)
        tmax = subtime[imax]
        if np.max(subdata) > threshold * RMS:
            maxima.append(np.max(subdata))
            tmaxima.append(tmax)
    if len(maxima) > 0:
        imax = np.argmax(np.array(maxima))
        tS = tmaxima[imax]
        maxima.remove(np.max(np.array(maxima)))
        tmaxima.remove(tS)
        if (len(maxima) > 0):
            imax = np.argmax(np.array(maxima))
            tP = tmaxima[imax]
            if tP > tS:
                (tP, tS) = (tS, tP)
            return (tP, tS)
        else:
            return (np.nan, tS)
    else:
        return (np.nan, np.nan)
    
def get_time_arrival(directory, family_file, station_file, t0, threshold):
    """
    """
    families = pd.read_csv(family_file, sep=r'\s{1,}', header=None, engine='python')
    families.columns = ['family', 'stations']

    staloc = pd.read_csv(station_file, sep=r'\s{1,}', header=None, engine='python')
    staloc.columns = ['station', 'network', 'channels', 'location', \
        'server', 'latitude', 'longitude', 'time_on', 'time_off']

    df = pd.DataFrame(columns=['family', 'station', 'channel', 'tP', 'tS'])

    for i in range(0, len(families)):
        family = families['family'].iloc[i]
        stations = families['stations'].iloc[i].split(',')
        namedir = directory + '/' + family
        for station in stations:
            for ir in range(0, len(staloc)):
                if (station == staloc['station'][ir] and staloc['location'][ir] == '--'):
                    channels = staloc['channels'][ir]
            channels = channels.split(',')
            for channel in channels:
                filename = namedir + '/' + station + '_' + channel + '.pkl'
                data = pickle.load(open(filename, 'rb'))
                template = data[0]
                dt = dt = template.stats.delta
                nt = template.stats.npts
                time = dt * np.arange(0, nt)
                (tP, tS) = pick_arrivals(time, template, t0, threshold)
                if (tS != np.nan):
                    i0 = len(df.index)
                    df.loc[i0] = [family, station, channel, tP, tS]
    return df

def get_origin_time(family_file, station_file, df):
    """
    """
    LFEloc = np.loadtxt(family_file, \
        dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
        'eZ', 'nb'), \
             'formats': ('S13', 'S3', np.float, np.float, np.float, \
        np.float, np.float, np.int)}, \
        skiprows=1)
    staloc = pd.read_csv(station_file, sep=r'\s{1,}', header=None, engine='python')
    staloc.columns = ['station', 'network', 'channels', 'location', \
        'server', 'latitude', 'longitude', 'time_on', 'time_off']
    df['distance'] = np.zeros(len(df))
    a = 6378.136
    e = 0.006694470
    # Construct matrix for S-wave arrival time
    df.dropna(subset=['tS'], inplace=True)
    df.reset_index(drop=True, inplace=True)
#    X = np.zeros((len(df), len(LFEloc) + 1))
    X = np.zeros((len(df), len(LFEloc)))
    Y = df['tS'].to_numpy()
    for i in range(0, len(df)):
        family = df['family'].iloc[i]
        station = df['station'].iloc[i]
        for j in range(0, len(LFEloc)):
            if family == LFEloc[j][0].decode('utf-8'):
                icol = j
                lat0 = LFEloc[j][2]
                lon0 = LFEloc[j][3]
        subdf = staloc.loc[staloc['station'] == station]
        lat = subdf['latitude'].iloc[0]
        lon = subdf['longitude'].iloc[0]
        dx = (pi / 180.0) * a * cos(lat0 * pi / 180.0) / sqrt(1.0 - e * e * \
            sin(lat0 * pi / 180.0) * sin(lat0 * pi / 180.0))
        dy = (3.6 * pi / 648.0) * a * (1.0 - e * e) / ((1.0 - e * e * sin(lat0 * \
            pi / 180.0) * sin(lat0 * pi / 180.0)) ** 1.5)
        x = dx * (lon - lon0)
        y = dy * (lat - lat0)
        distance = sqrt(x ** 2.0 + y ** 2.0)
        df.at[i, 'distance'] = distance
        X[i, icol] = 1
#        X[i, len(LFEloc)] = distance
        Y[i] = Y[i] - distance * 0.25
    # Solve linear system
    regr = linear_model.LinearRegression(fit_intercept=False)
    regr.fit(X, Y)
    tori = regr.coef_[0 : len(LFEloc)]
#    sS = regr.coef_[len(LFEloc)]
    # Construct matrix for P-wave arrival time
#    df.dropna(subset=['tP'], inplace=True)
#    df.drop(df[df.tS - df.tP < 0.125 * df.distance].index, inplace=True)
#    df.reset_index(drop=True, inplace=True)
#    X = np.zeros((len(df), 1))
#    Y = np.zeros(len(df))
#    for i in range(0, len(df)):
#        family = df['family'].iloc[i]
#        for j in range(0, len(LFEloc)):
#            if family == LFEloc[j][0].decode('utf-8'):
#                icol = j
#        X[i] = df['distance'].iloc[i]
#        Y[i] = df['tP'].iloc[i] - tori[icol]
    # Solve linear system
#    regr = linear_model.LinearRegression(fit_intercept=False)
#    regr.fit(X, Y)
#    sP = regr.coef_[0]
#    return (tori, sS, sP)
    return tori

if __name__ == '__main__':

    directory = 'templates_2007_2009'
    family_file = '../data/Ducellier/families_permanent.txt'
    station_file = '../data/Ducellier/stations_permanent.txt'
    t0 = 1.0
    threshold = 3.0
    df = get_time_arrival(directory, family_file, station_file, t0, threshold)

    family_file = '../data/Plourde_2015/templates_list.txt'
#    (tori, sS, sP) = get_origin_time(family_file, station_file, df)
    tori = get_origin_time(family_file, station_file, df)
#    pickle.dump((tori, sS, sP), open('tori.pkl', 'wb'))
    pickle.dump(tori, open('tori.pkl', 'wb'))
    