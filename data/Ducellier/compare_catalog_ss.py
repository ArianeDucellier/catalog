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

    # Shift FAME catalog to the right
    df1r = df1.copy()
    df1rr = df1.copy()
    df1rrr = df1.copy()
    df1r['second'] = df1r['second'] + 0.05
    df1rr['second'] = df1rr['second'] + 0.10
    df1rrr['second'] = df1rrr['second'] + 0.15

    # Shift FAME catalog to the left
    df1l = df1.copy()
    df1ll = df1.copy()
    df1lll = df1.copy()
    df1l['second'] = df1l['second'] - 0.05
    df1ll['second'] = df1ll['second'] - 0.10
    df1lll['second'] = df1lll['second'] - 0.15

    # Shift network catalog to the right
    df2r = df2.copy()
    df2rr = df2.copy()
    df2rrr = df2.copy()
    df2r['second'] = df2r['second'] + 0.05
    df2rr['second'] = df2rr['second'] + 0.10
    df2rrr['second'] = df2rrr['second'] + 0.15

    # Shift network catalog to the left
    df2l = df2.copy()
    df2ll = df2.copy()
    df2lll = df2.copy()
    df2l['second'] = df2l['second'] - 0.05
    df2ll['second'] = df2ll['second'] - 0.10
    df2lll['second'] = df2lll['second'] - 0.15

    # LFEs in filtered FAME catalog but not in (unfiltered) network catalog
    df2_all = pd.concat([df2, df2r, df2rr, df2rrr, df2l, df2ll, df2lll], ignore_index=True)
    df_all = pd.concat([df2_all, df1_filter], ignore_index=True)
    df_merge = df_all.merge(df2_all.drop_duplicates(), \
        on=['date'], \
        how='left', indicator=True)
    df_added_FAME = df_merge[df_merge['_merge'] == 'left_only']

    # LFEs in filtered network catalog but not in (unfiltered) FAME catalog
    df1_all = pd.concat([df1, df1r, df1rr, df1rrr, df1l, df1ll, df1lll], ignore_index=True)
    df_all = pd.concat([df1_all, df2_filter], ignore_index=True)
    df_merge = df_all.merge(df1_all.drop_duplicates(), \
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
    thresh1 = 4.4
    thresh2 = 1.9

    for i in range(63, 64): #len(families)):
        family = families['family'].iloc[i]
        dt = families['dt'].iloc[i]
        compare_catalog(family, tbegin, tend, dt, thresh1, thresh2)
