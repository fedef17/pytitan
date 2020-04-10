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

# Carico risultati analisi

wl_med, stat_latsza, ok_szas, coeffs_latsza, cub_latsza, pix_latsza, bestsza, pix_bestsza, stat_all, median_all, mean_all, std_all = pickle.load(open(cart_out + 'stat_pixs_2013-2017_v2.p', 'rb'))

bestsza_data = pickle.load(open(cart_out + 'bestsza_data_2013-2017_v2.p', 'rb'))
allsza_data = pickle.load(open(cart_out + 'allsza_data_2013-2017_v2.p', 'rb'))


fillo = open(cart_out + 'stat_lat_szamin_bestcube.dat', 'w')
forma2 = '{:6.0f}'*18 + '\n'
fillo.write('Best sza available for each box \n')
fillo.write('anno  |  ' + forma.format(*map(int, latgrid)))
for iye, ye in enumerate(range(2013, 2018)):
    cose = []
    for latg in latgrid:
        if (ye, latg) in bestsza.keys():
            cose.append(bestsza[(ye, latg)])
        else:
            cose.append(np.nan)
    fillo.write(str(ye) +'  |  '+ forma2.format(*cose))

fillo.write('\n Max number of spectra in the same cube at best sza available for each box \n')

for iye, ye in enumerate(range(2013, 2018)):
    cose = []
    for latg in latgrid:
        if (ye, latg) in bestsza_data.keys():
            cose.append(len(bestsza_data[(ye, latg)]))
        else:
            cose.append(0)
    fillo.write(str(ye) +'  |  '+ forma.format(*cose))

fillo.close()
