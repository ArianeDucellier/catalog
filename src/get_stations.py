import pandas as pd

from datetime import datetime

def f(x):
    """
    Concatenate channels
    """
    result = '%s' % ','.join(x)
    result = list(set(result.split(',')))
    result = '%s' % ','.join(result)
    return result

def filter_stations(network):
    """
    Keep only relevant stations from the network
    """
    df = pd.read_csv('../data/networks/' + network + '.txt', sep='\t')

    # Filter channels
    channels = ['BHE', 'BHN', 'BHZ', 'BH1', 'BH2', \
                'EHE', 'EHN', 'EHZ', 'EH1', 'EH2', \
                'HHE', 'HHN', 'HHZ', 'HH1', 'HH2', \
                'SHE', 'SHN', 'SHZ', 'SH1', 'SH2']
    df = df.loc[df['Cha'].isin(channels)]

    # Filter recording dates
    startdate = pd.to_datetime(df['Start time'], format='%Y/%m/%d,%H:%M:%S')
    enddate = df['End time'].str.replace('3000', '2025')
    enddate = pd.to_datetime(enddate, format='%Y/%m/%d,%H:%M:%S')
    df['Start time'] = startdate
    df['End time'] = enddate
    mask = (df['Start time'] <= datetime(2009, 7, 1)) & (df['End time'] >= datetime(2007, 7, 1))
    df = df.loc[mask]

    # Group channels
    df = df.groupby(['Stat', 'Net', 'Lo', 'Lat', 'Lon', 'Elev', 'Depth']).agg({ \
        'Cha':f, \
        'Start time':lambda x: min(x), \
        'End time':lambda x: max(x)}).reset_index()

    return df
