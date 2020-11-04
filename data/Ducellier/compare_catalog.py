import pandas as pd
import pickle

from datetime import datetime, timedelta

def compare_catalog(family, tbegin, tend, dt, thresh1, thresh2):
    """
    """
    # Read FAME catalog
    namefile = 'catalogs/' + family + '/catalog_2007_2009.pkl'
    df1 = pickle.load(open(namefile, 'rb'))
    df1 = df1[['year', 'month', 'day', 'hour', 'minute', 'second', \
        'cc', 'nchannel']]
    df1 = df1.astype({'year': int, 'month': int, 'day': int, \
        'hour': int, 'minute': int, 'second': float, \
        'cc': float, 'nchannel': int})
    date = pd.to_datetime(df1.drop(columns=['cc', 'nchannel']))
    df1['date'] = date
    df1 = df1[(df1['date'] >= tbegin) & (df1['date'] <= tend)]
    df1_filter = df1.loc[df1['cc'] * df1['nchannel'] >= thresh1]

    # Read network catalog
    namefile = 'catalogs/' + family + '/catalog_2004_2011.pkl'
    df2 = pickle.load(open(namefile, 'rb'))
    df2 = df2[['year', 'month', 'day', 'hour', 'minute', 'second', \
        'cc', 'nchannel']]
    df2 = df2.astype({'year': int, 'month': int, 'day': int, \
        'hour': int, 'minute': int, 'second': float, \
        'cc': float, 'nchannel': int})
    date = pd.to_datetime(df2.drop(columns=['cc', 'nchannel']))
    df2['date'] = date
    df2['date'] = df2['date'] - timedelta(seconds=dt)
    df2 = df2[(df2['date'] >= tbegin) & (df2['date'] <= tend)]
    df2_filter = df2.loc[df2['cc'] * df2['nchannel'] >= thresh2]

    # LFEs in filtered FAME catalog but not in (unfiltered) network catalog
    df_all = pd.concat([df2, df1_filter], ignore_index=True)
    df_merge = df_all.merge(df2.drop_duplicates(), \
        on=['date'], \
        how='left', indicator=True)
    df_added_FAME = df_merge[df_merge['_merge'] == 'left_only']

    # LFEs in filtered network catalog but not in (unfiltered) FAME catalog
    df_all = pd.concat([df1, df2_filter], ignore_index=True)
    df_merge = df_all.merge(df1.drop_duplicates(), \
        on=['date'], \
        how='left', indicator=True)
    df_added_network = df_merge[df_merge['_merge'] == 'left_only']
    
    print(family)
    print('FAME: false detections = {}, total = {}, ratio = {}'.format( \
         len(df_added_FAME), len(df1_filter),len(df_added_FAME) / len(df1_filter)))
    print('Network: false detections = {}, total = {}, ratio = {}'.format( \
         len(df_added_network), len(df2_filter),len(df_added_network) / len(df2_filter)))

if __name__ == '__main__':

    # Get the families, stations, and duration of the template
    families = pd.read_csv('families_permanent_timeinterval.txt', sep=r'\s{1,}', header=None, \
        engine='python')
    families.columns = ['family', 'stations', 'dt', 'unused']

    # Beginning and end of the period we are looking at
    tbegin = datetime(2007, 9, 25, 0, 0, 0)
    tend = datetime(2009, 5, 14, 0, 0, 0)

    # Thresholds
    thresh1 = 2.6
    thresh2 = 1.2

    for i in range(21, 22): #len(families)):
        family = families['family'].iloc[i]
        dt = families['dt'].iloc[i]
        compare_catalog(family, tbegin, tend, dt, thresh1, thresh2)
