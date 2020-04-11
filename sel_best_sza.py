#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import pickle

import numpy as np
import matplotlib.pyplot as plt

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

# print('tutti i cubi disponibili: ', allcubi)
# allcubi = np.unique(np.concatenate(datasav.cubo))
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
stat_latsza = dict()
coeffs_latsza = dict()
cub_latsza = dict()
pix_latsza = dict()

ok_szas = dict()

bestsza = dict()
pix_bestsza = dict()

#np.empty((len(years), len(latgrid), len(szagrid)))
wl_med = dict()#np.zeros((len(years), wldim))

stat_all = dict()
median_all = dict() #np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))
mean_all = dict() #np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))
std_all = dict()#np.zeros((len(years), len(latgrid), len(szagrid), len(altgrid), wldim))

#for iye, ye in enumerate(range(2013, 2018)):
for iye, ye in enumerate([2017]):
    cond_time = (datasav.year >= ye) & (datasav.year < ye+1)
    print(ye, np.sum(cond_time))

    allwl = np.stack([wlo[okwl] for wlo in datasav.wl[cond_time]])
    #wl_med[iye, :] = np.median(allwl, axis = 0)
    wl_med[ye] = np.median(allwl, axis = 0)

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

            if totsza == 0:
                stat_latsza[(ye, latg, szag)] = 0
                continue

            cond_tot = np.all([cond_time, cond_lat, cond_sza], axis = 0)

            # Check sui cubi
            allcubi = np.concatenate(datasav.cubo[cond_tot])
            tot_alts = np.concatenate(datasav.alt[cond_tot])
            tot_szas = np.concatenate(datasav.sza[cond_tot])

            unicub = np.unique(allcubi)
            tutticubs = []
            for cub in unicub:
                okcub = allcubi == cub
                if np.min(tot_alts[okcub]) < 550. and np.max(tot_alts[okcub]) > 950 and np.sum(okcub) > 8:
                    tutticubs.append((cub, np.sum(okcub)))
                    print(cub, np.sum(okcub))

            if len(tutticubs) == 0:
                stat_latsza[(ye, latg, szag)] = 0
                continue
            else:
                stat_latsza[(ye, latg, szag)] = np.sum(cond_tot)
                cub_latsza[(ye, latg, szag)] = unicub
                pix_latsza[(ye, latg, szag)] = np.where(cond_tot)[0]
                ok_szas[(ye, latg)].append(np.mean(np.concatenate(datasav.sza[cond_tot])))

            # metto in bestsza il cubo più abbondante in questo sza
            n_cubs = [cos[1] for cos in tutticubs]
            best_cub = tutticubs[np.argmax(n_cubs)][0]

            cond_cub = datasav.cubo == best_cub # adding this condition
            cond_tot_cub = np.all([cond_time, cond_lat, cond_sza, cond_cub], axis = 0)
            print(tutticubs)
            print('best', best_cub, np.sum(cond_tot_cub))

            if len(ok_szas[(ye, latg)]) == 1:
                # this is the best sza
                pix_bestsza[(ye, latg)] = np.where(cond_tot_cub)[0]
                bestsza[(ye, latg)] = np.mean(np.concatenate(datasav.sza[cond_tot_cub]))

                if latg == 5.0: sys.exit()

            allspet = np.stack([spe[okwl] for spe in datasav.spet[cond_tot]])
            allalts = np.concatenate(datasav.alt[cond_tot])

            coeffs_latsza[(ye, latg, szag)] = []
            for iwl, wla in enumerate(wl_med[ye]):
                fitco = np.polyfit(allalts, allspet[:, iwl], deg = 3, cov = False)
                coeffs_latsza[(ye, latg, szag)].append((wla, fitco))

            for ialt, (alt1, alt2, altg) in enumerate(zip(alts[:-1], alts[1:], altgrid)):
                cond_alt = (datasav.alt >= alt1) & (datasav.alt < alt2)

                cond_tot_alt = np.all([cond_tot, cond_alt], axis = 0)
                nutot = np.sum(cond_tot_alt)
                #print('measures in alt box {} - {}: {}'.format(alt1, alt2, nutot))

                stat_all[(ye, latg, szag, altg)] = nutot
                if nutot == 0: continue

                allspet = np.stack([spe[okwl] for spe in datasav.spet[cond_tot_alt]])

                median_all[(ye, latg, szag, altg)] = np.median(allspet, axis = 0)
                mean_all[(ye, latg, szag, altg)] = np.mean(allspet, axis = 0)
                std_all[(ye, latg, szag, altg)] = np.std(allspet, axis = 0)

pickle.dump([wl_med, stat_latsza, ok_szas, coeffs_latsza, cub_latsza, pix_latsza, bestsza, pix_bestsza, stat_all, median_all, mean_all, std_all], open(cart_out + 'stat_pixs_2013-2017_v2.p', 'wb'))

# Salvo un datasav ridotto ai soli cosi che ci interessano
bestdata = dict()

for ke in pix_bestsza:
    gigi = datasav[pix_bestsza[ke]]
    for att in allattrs:
        if len(gigi.__getattribute__(att)[0]) == 1:
            gigi.__setattr__(att, np.concatenate(gigi.__getattribute__(att)))
    bestdata[ke] = gigi

pickle.dump(bestdata, open(cart_out + 'bestsza_data_2013-2017_v2.p', 'wb'))


# Salvo un datasav ridotto al miglior cubo di ogni sza
bestdata = dict()

for ke in pix_latsza:
    gigi = datasav[pix_latsza[ke]]
    for att in allattrs:
        if len(gigi.__getattribute__(att)[0]) == 1:
            gigi.__setattr__(att, np.concatenate(gigi.__getattribute__(att)))
    bestdata[ke] = gigi

pickle.dump(bestdata, open(cart_out + 'allsza_data_2013-2017_v2.p', 'wb'))
