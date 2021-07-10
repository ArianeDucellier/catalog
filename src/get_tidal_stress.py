"""
Script to read tidal stress file from H. Houston and
interpolate the values at the LFE locations
"""

import numpy as np
import pickle

from math import cos, pi, sin, sqrt
from scipy.io import loadmat

# Earth's radius and ellipticity
a = 6378.136
e = 0.006694470

# List of LFE families
templates = np.loadtxt('../data/Plourde_2015/templates_list.txt', \
    dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
    'eZ', 'nb'), \
         'formats': ('S13', 'S3', np.float, np.float, np.float, \
    np.float, np.float, np.int)}, \
    skiprows=1)

# Read tidal stress file
data = loadmat('../data/For_Ken_California_LFE/stressonfp_4Ken_extended.mat')
xmnr = data['xmnr']
ymnr = data['ymnr']

# Loop on templates
for i in range(0, np.shape(templates)[0]):

    lat0 = templates[i][2]
    lon0 = templates[i][3]
    dx = (pi / 180.0) * a * cos(lat0 * pi / 180.0) / sqrt(1.0 - e * e * \
        sin(lat0 * pi / 180.0) * sin(lat0 * pi / 180.0))
    dy = (3.6 * pi / 648.0) * a * (1.0 - e * e) / ((1.0 - e * e * sin(lat0 * \
        pi / 180.0) * sin(lat0 * pi / 180.0)) ** 1.5)
    x = (xmnr - lon0) * dx
    y = (ymnr - lat0 ) * dy
    dist = np.sqrt(np.power(x, 2.0) + np.power(y, 2.0))
    closest_dists = np.sort(dist, axis=None)[0:4]
    indices = np.argsort(dist, axis=None)[0:4]
    indices = np.unravel_index(indices, dist.shape)
    ind1 = indices[0]
    ind2 = indices[1]
    weights = np.reshape(closest_dists / np.sum(closest_dists), (4, 1))
    norm = np.sum(weights * data['gridnorm'][ind1, ind2, :], axis=0)
    shear = np.sum(weights * data['gridshear'][ind1, ind2, :], axis=0)
    vol = np.sum(weights * data['gridvol'][ind1, ind2, :], axis=0)
    time = data['tvec']
    filename = '../data/tidal_stress/' + templates[i][0].astype(str) + '.pkl'
    pickle.dump([time, norm, shear, vol], open(filename, 'wb'))
