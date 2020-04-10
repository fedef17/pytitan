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

#cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/titano_leonardo/'
cart_in = '/home/federico/TITANO/VIMS_data/NEW_COLL_HCN-CH4-C2H2_sza110_limb10/'

cart_out = '/home/federico/TITANO/prescreening_2013-19/'
if not os.path.exists(cart_out):
    os.mkdir(cart_out)

# nomecub = 'C1809689567_1_ir.cub'
# cubo = plim.CubeFile.open(cart_in + nomecub)

#pixs = io.readsav(cart_in + 'PIXs_HCN-CH4-C2H2_season_sza90_all_nu2.sav')
pixs = io.readsav(cart_in + 'PIXs_HCN-CH4-C2H2_sza110_limb10.sav')
#pixs = io.readsav(cart_in + 'PIXs_2013-2017_RC19.sav')

datasav = pixs['comppix']
allattrs = datasav.dtype.names # questa è la lista di tutti gli attributi del dataset
print('all attributes: ', allattrs)

allcubi = np.unique(np.concatenate(datasav.cubo))
print('tutti i cubi disponibili: ', allcubi)
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

for iye, ye in enumerate(range(2013, 2018)):
    cond_time = (datasav.year >= ye) & (datasav.year < ye+1)
    print(ye, np.sum(cond_time))

    allwl = np.stack([wlo[okwl] for wlo in datasav.wl[cond_time]])
    wl_med[iye, :] = np.median(allwl, axis = 0)

    for ilat, (lat1, lat2, latg) in enumerate(zip(lats[:-1], lats[1:], latgrid)):
        cond_lat = (datasav.lat >= lat1) & (datasav.lat < lat2)
        totlat = np.sum(np.all([cond_time, cond_lat], axis = 0))
        print('measures in latitude box {} - {}: {}'.format(lat1, lat2, totlat))

        if totlat == 0: continue

        ok_szas[(ye, latg)] = []

        for isza, (sza1, sza2, szag) in enumerate(zip(szas[:-1], szas[1:], szagrid)):
            cond_sza = (datasav.sza >= sza1) & (datasav.sza < sza2)
            totsza = np.sum(np.all([cond_time, cond_lat, cond_sza], axis = 0))
            print('measures in sza box {} - {}: {}'.format(sza1, sza2, totsza))

            if totsza == 0: continue

            cond_tot = np.all([cond_time, cond_lat, cond_sza], axis = 0)
            stat_latsza[iye, ilat, isza] = np.sum(cond_tot)

            tot_alts = np.concatenate(datasav.alt[cond_tot])
            if np.min(tot_alts) < 500. and np.max(tot_alts) > 1000 and totsza > 10:
                ok_szas[(ye, latg)].append(np.mean(np.concatenate(datasav.sza[cond_tot])))

            allspet = np.stack([spe[okwl] for spe in datasav.spet[cond_tot]])
            allalts = np.concatenate(datasav.alt[cond_tot])

            coeffs_latsza[(ye, latg, szag)] = []
            for iwl, wla in enumerate(wl_med[iye, :]):
                fitco = np.polyfit(allalts, allspet[:, iwl], deg = 3, cov = False)
                coeffs_latsza[(ye, latg, szag)].append(fitco)

            for ialt, (alt1, alt2) in enumerate(zip(alts[:-1], alts[1:])):
                cond_alt = (datasav.alt >= alt1) & (datasav.alt < alt2)

                cond_tot = np.all([cond_time, cond_lat, cond_sza, cond_alt], axis = 0)
                nutot = np.sum(cond_tot)
                #print('measures in alt box {} - {}: {}'.format(alt1, alt2, nutot))

                if nutot == 0: continue

                allspet = np.stack([spe[okwl] for spe in datasav.spet[cond_tot]])

                stat_all[iye, ilat, isza, ialt] = nutot
                median_all[iye, ilat, isza, ialt, :] = np.median(allspet, axis = 0)
                mean_all[iye, ilat, isza, ialt, :] = np.mean(allspet, axis = 0)
                std_all[iye, ilat, isza, ialt, :] = np.std(allspet, axis = 0)


pickle.dump([wl_med, stat_latsza, ok_szas, coeffs_latsza, stat_all, median_all, mean_all, std_all], open(cart_out + 'stat_pixs_2013-2019_szastep5.p', 'wb'))


szas = np.array([0, 40, 50, 60, 70, 75, 80, 85, 90, 95, 100, 105, 110])
szagrid = np.mean([szas[:-1], szas[1:]], axis = 0)

stat_latsza = np.zeros((len(years), len(latgrid), len(szagrid)))
coeffs_latsza = dict()
#np.empty((len(years), len(latgrid), len(szagrid)))
wl_med = np.empty((len(years), wldim))

stat_all = np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid)))
median_all = np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))
mean_all = np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))
std_all = np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))

for iye, ye in enumerate(range(2013, 2018)):
    cond_time = (datasav.year >= ye) & (datasav.year < ye+1)
    print(ye, np.sum(cond_time))

    allwl = np.stack([wlo[okwl] for wlo in datasav.wl[cond_time]])
    wl_med[iye, :] = np.median(allwl, axis = 0)

    for ilat, (lat1, lat2, latg) in enumerate(zip(lats[:-1], lats[1:], latgrid)):
        cond_lat = (datasav.lat >= lat1) & (datasav.lat < lat2)
        totlat = np.sum(np.all([cond_time, cond_lat], axis = 0))
        print('measures in latitude box {} - {}: {}'.format(lat1, lat2, totlat))

        if totlat == 0: continue

        for isza, (sza1, sza2, szag) in enumerate(zip(szas[:-1], szas[1:], szagrid)):
            cond_sza = (datasav.sza >= sza1) & (datasav.sza < sza2)
            totsza = np.sum(np.all([cond_time, cond_lat, cond_sza], axis = 0))
            print('measures in sza box {} - {}: {}'.format(sza1, sza2, totsza))

            if totsza == 0: continue

            cond_tot = np.all([cond_time, cond_lat, cond_sza], axis = 0)
            stat_latsza[iye, ilat, isza] = np.sum(cond_tot)

            allspet = np.stack([spe[okwl] for spe in datasav.spet[cond_tot]])
            allalts = np.concatenate(datasav.alt[cond_tot])

            coeffs_latsza[(ye, latg, szag)] = []
            for iwl, wla in enumerate(wl_med[iye, :]):
                fitco = np.polyfit(allalts, allspet[:, iwl], deg = 3, cov = False)
                coeffs_latsza[(ye, latg, szag)].append(fitco)

            for ialt, (alt1, alt2) in enumerate(zip(alts[:-1], alts[1:])):
                cond_alt = (datasav.alt >= alt1) & (datasav.alt < alt2)

                cond_tot = np.all([cond_time, cond_lat, cond_sza, cond_alt], axis = 0)
                nutot = np.sum(cond_tot)
                #print('measures in alt box {} - {}: {}'.format(alt1, alt2, nutot))

                if nutot == 0: continue

                allspet = np.stack([spe[okwl] for spe in datasav.spet[cond_tot]])

                stat_all[iye, ilat, isza, ialt] = nutot
                median_all[iye, ilat, isza, ialt, :] = np.median(allspet, axis = 0)
                mean_all[iye, ilat, isza, ialt, :] = np.mean(allspet, axis = 0)
                std_all[iye, ilat, isza, ialt, :] = np.std(allspet, axis = 0)


pickle.dump([wl_med, stat_latsza, coeffs_latsza, stat_all, median_all, mean_all, std_all], open(cart_out + 'stat_pixs_2013-2019_lesssza.p', 'wb'))

# Convert to dictionary


# Create a DataArray with all cube data: alt, lat, lon, sza, pha
# Then a Dataset
