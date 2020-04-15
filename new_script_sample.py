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

import lib_vims_select as lvs
############################################################

#result = lvs.read_obs('osservato.dat')

cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/titano_leonardo/'

cart_out = cart_in + 'figures/'
if not os.path.exists(cart_out):
    os.mkdir(cart_out)

wl_med, stat_latsza, ok_szas, coeffs_latsza, cub_latsza, pix_latsza, bestsza, pix_bestsza, stat_all, median_all, mean_all, std_all = pickle.load(open(cart_in + 'stat_pixs_2013-2017_v2.p', 'rb'))

bestsza_data = pickle.load(open(cart_in + 'bestsza_data_2013-2017_v2.p', 'rb'))

# Definisco delle scatole
years = np.arange(2013, 2018)

lats = np.arange(-90, 91, 10)
latgrid = np.mean([lats[:-1], lats[1:]], axis = 0)

szas = np.arange(0, 111, 5)
szagrid = np.mean([szas[:-1], szas[1:]], axis = 0)

alts = np.arange(350, 1101, 50)
altgrid = np.mean([alts[:-1], alts[1:]], axis = 0)

# bestsza_data contiene tutti i dati che ci servono (spero)
# for ye in years:
#     wlgrid = wl_med[ye] # le wl per questo anno
#     allfigs = []
#     filename = cart_out + 'fig_bestspectra_{}.pdf'.format(ye)
#     for latg in latgrid:
#         if (ye, latg) in bestsza_data.keys(): # ci sono dati?
#             fig = plt.figure()
#             plt.title('year {} - lat {}: sza {:6.1f}'.format(ye, latg, bestsza[(ye, latg)]))
#             for pixel in bestsza_data[(ye, latg)]:
#                 plt.plot(wlgrid, pixel.spet)#, label = pixel.alt)
#             #plt.legend()
#             #fig.savefig(cart_fig + 'spet_{}_{}.pdf'.format(ye, int(latg)))
#             allfigs.append(fig)
#
#     lvs.plot_pdfpages(filename, allfigs)
#
#     wlgrid = wl_med[ye] # le wl per questo anno
#     allfigs = []
#     filename = cart_out + 'medianspectra_{}.pdf'.format(ye)
#     for latg in latgrid:
#         if (ye, latg) in bestsza_data.keys(): # ci sono dati?
#             fig = plt.figure()
#             plt.title('year {} - lat {}: sza {:6.1f}'.format(ye, latg, bestsza[(ye, latg)]))
#             szag = szagrid[np.argmin(np.abs(szagrid - bestsza[(ye, latg)]))]
#
#             for altg in altgrid:
#                 if (ye, latg, szag, altg) in median_all.keys():
#                     okmedian = median_all[(ye, latg, szag, altg)]
#                     plt.plot(wlgrid, okmedian)#, label = pixel.alt)
#             #plt.legend()
#             #fig.savefig(cart_fig + 'spet_{}_{}.pdf'.format(ye, int(latg)))
#             allfigs.append(fig)
#
#     lvs.plot_pdfpages(filename, allfigs)


# prova out observ.dat
lvs.write_obs_frompixels(cart_out + 'test_observ.dat', bestsza_data[(2015, 15.0)])

cart_out = cart_in
# stampo cubo per ogni lat
forma = '{:12d}'*18 + '\n'

fillo = open(cart_out + 'stat_annolat_bestcube.dat', 'w')
forma2 = '{:12s}'*18 + '\n'
fillo.write('Best cube available for each box \n')
fillo.write('anno  |  ' + forma.format(*map(int, latgrid)))
for iye, ye in enumerate(range(2013, 2018)):
    cose = []
    for latg in latgrid:
        if (ye, latg) in bestsza.keys():
            cose.append(bestsza_data[(ye, latg)].cubo[0])
        else:
            cose.append('none')
    print(cose)
    print(forma2.format(*cose))
    fillo.write(str(ye) +'  |  '+ forma2.format(*cose))

fillo.close()

fillo = open(cart_out + 'sequenze_best_2013-2017.dat', 'w')
forma2 = '{:12s}'*18 + '\n'


uadists = [9.9, 9.9, 10.0, 10.0, 10.1]
# _7892_50NW      2008 fliby  7892  50.11   66.52  305130.0  9.305
fillo.write('Best VIMS sequences available in years 2013-2017\n')
stringa = '{:9s}'.format('nome') + '  {:4s}  {:4s}  {:11s}  {:6s}  {:6s}  {:8s}  {}\n'.format('anno', 'flyby', 'cubo', 'latm', 'szam', 'titdist(km)', 'sundist(ua)')
fillo.write(stringa)
for ye, uadist in zip(range(2013, 2018), uadists):
    cose = []
    for latg in latgrid:
        if (ye, latg) in bestsza.keys():
            cubo = bestsza_data[(ye, latg)].cubo[0]
            latok = np.mean(bestsza_data[(ye, latg)].lat)
            szaok = np.mean(bestsza_data[(ye, latg)].sza)
            dist = np.mean(bestsza_data[(ye, latg)].dist)
            if latg < 0:
                nome = '{:4s}_{:02d}S'.format(cubo[-4:], int(abs(latg)))
            else:
                nome = '{:4s}_{:02d}N'.format(cubo[-4:], int(abs(latg)))
            flyby = ''
            stringa = nome + '  {:4d}  {:4s}  {:11s}  {:6.2f}  {:6.2f}  {:8.3e}  {:6.1f}\n'.format(ye, flyby, cubo, latok, szaok, dist, uadist)
            fillo.write(stringa)
fillo.close()
