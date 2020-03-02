"""
This module contains functions to get instrument response
from the IRIS DMC or the NCEDC website
"""

import obspy
import obspy.clients.fdsn.client as fdsn
import urllib.request

def get_from_IRIS(station, network):
    """
    """
    fdsn_client = fdsn.Client('IRIS')
    inventory = fdsn_client.get_stations(network=network, \
        station=station, level='response')
    inventory.write('../data/response/' + network + '_' + station + '.xml', \
        format='STATIONXML')

def get_from_NCEDC(station, network):
    """
    """
    url = 'http://service.ncedc.org/fdsnws/station/1/query?net=' + \
        network + '&sta=' + station + \
        '&level=response&format=xml&includeavailability=true'
    s = urllib.request.urlopen(url)
    contents = s.read()
    file = open('../data/response/' + network + '_' + station + '.xml', 'wb')
    file.write(contents)
    file.close()

if __name__ == '__main__':

    station_file = '../data/Ducellier/stations_temporary.txt'

    # Get the network, channels, and location of the stations
    staloc = pd.read_csv(station_file, sep=r'\s{1,}', header=None, engine='python')
    staloc.columns = ['station', 'network', 'channels', 'location', \
        'server', 'latitude', 'longitude', 'time_on', 'time_off']

    # Download the response of the stations
    for ir in range(0, len(staloc)):
        station = staloc['station'][ir]
        network = staloc['network'][ir]
        server = staloc['server'][ir]
        # First case: we can get the data from IRIS
        if (server == 'IRIS'):
            get_from_IRIS(station, network)
        # Second case: we get the data from NCEDC
        elif (server == 'NCEDC'):
            get_from_NCEDC(station, network)
        else:
            raise ValueError('You can only download data from IRIS and NCEDC')
