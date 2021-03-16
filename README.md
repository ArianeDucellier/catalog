# Catalog

This repository contains Python scripts to read the initial LFEs times from the catalog from Plourde et al. (2015), compute the corresponding templates, and find new LFEs for the period 2004-2011.

I have developed two catalogs:
- The first catalog was obtained with the data from the FAME experiment (2007-2009).
- The second catalog was obtained with the data from the permanent seismic networks (2004-2011).

The catalogs are in data/Ducellier/catalogs. There is one directory for each LFE family. Each directory contains two files corresponding to the two catalogs. Each file contains a pandas data frame with the following columns: year, month, day, hour, minute, second (of the event), cc (average cross-correlation between waveform and template), nchannel (number of seismic channels recording at that time).
