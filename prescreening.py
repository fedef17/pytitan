#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import pickle

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#import planetaryimage as plim
#from scipy.misc import imsave    # requires pillow as well
import scipy.io as io

############################################################

cart_in = '/home/federico/TITANO/VIMS_data/NEW_COLL_HCN-CH4-C2H2_sza110_limb10/'

cart_out = '/home/federico/TITANO/prescreening_2013-19/'
if not os.path.exists(cart_out):
    os.mkdir(cart_out)

pixs = io.readsav(cart_in + 'PIXs_HCN-CH4-C2H2_sza110_limb10.sav')

datasav = pixs['comppix']
allattrs = datasav.dtype.names # questa è la lista di tutti gli attributi del dataset
print('all attributes: ', allattrs)

# allcubi = np.unique(np.concatenate(datasav.cubo))
# print('tutti i cubi disponibili: ', allcubi)
# for cub in allcubi:
#     print(cub, np.sum(datasav.cubo == cub))

# Definisco delle scatole
years = np.arange(2013, 2018)

lats = np.arange(-90, 91, 10)
latgrid = np.mean([lats[:-1], lats[1:]], axis = 0)

# szas = np.array(0, 40, 60, 70, 80, 90)
szas = np.arange(0, 111, 5)
szagrid = np.mean([szas[:-1], szas[1:]], axis = 0)

alts = np.arange(350, 1101, 50)
altgrid = np.mean([alts[:-1], alts[1:]], axis = 0)

# Range di wl da considerare:
wl0 = 2700.
wl1 = 3700.

allwl_med = np.median(np.stack(datasav.wl), axis = 0)
okwl = (allwl_med >= wl0) & (allwl_med <= wl1)
wldim = np.sum(okwl)

##
stat_latsza = np.zeros((len(years), len(latgrid), len(szagrid)))
coeffs_latsza = dict()
#np.empty((len(years), len(latgrid), len(szagrid)))
wl_med = np.zeros((len(years), wldim))

ok_szas = dict()

stat_all = np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid)))
median_all = np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))
mean_all = np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))
std_all = np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))

# Carico risultati analisi

wl_med, stat_latsza, ok_szas, coeffs_latsza, stat_all, median_all, mean_all, std_all = pickle.load(open(cart_out + 'stat_pixs_2013-2019_szastep5.p', 'rb'))

fillo = open(cart_out + 'stat_lat.dat', 'w')
forma = '{:6d}'*18 + '\n'
statlat = stat_latsza.sum(axis = -1)
fillo.write('Number of spectra at all szas available for each year/lat \n')
fillo.write('anno  |  ' + forma.format(*map(int, latgrid)))
for iye, ye in enumerate(range(2013, 2018)):
    fillo.write(str(ye) +'  |  '+ forma.format(*map(int, statlat[iye, :])))
fillo.close()


bestsza_year = dict()
nbest_year = dict()
for iye, ye in enumerate(np.arange(2013, 2018)):
    bestsza = []
    nbest = []
    for ilat, lat in enumerate(np.arange(-85.0, 86.0, 10.)):
        if (ye, lat) in ok_szas.keys():
            szaok = ok_szas[(ye, lat)]
            if len(szaok) == 0:
                nbest.append(0)
                bestsza.append(np.nan)
                continue
            szaok = np.min(szaok)
            isza = np.argmin(np.abs(szagrid-szaok))
            print(isza, szaok, szagrid)
            print(ye, lat, szaok, stat_latsza[(iye, ilat, isza)])
            nbest.append(stat_latsza[(iye, ilat, isza)])
            bestsza.append(szaok)
        else:
            nbest.append(0)
            bestsza.append(np.nan)

    bestsza_year[ye] = bestsza
    nbest_year[ye] = nbest

fillo = open(cart_out + 'stat_lat_szamin.dat', 'w')
forma2 = '{:6.0f}'*18 + '\n'
fillo.write('Best sza available for each box \n')
fillo.write('anno  |  ' + forma.format(*map(int, latgrid)))
for iye, ye in enumerate(range(2013, 2018)):
    fillo.write(str(ye) +'  |  '+ forma2.format(*bestsza_year[ye]))

fillo.write('\n Number of spectra at best sza available for each box \n')

for iye, ye in enumerate(range(2013, 2018)):
    fillo.write(str(ye) +'  |  '+ forma.format(*map(int, nbest_year[ye])))

fillo.close()
