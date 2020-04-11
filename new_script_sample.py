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
for ye in years:
    wlgrid = wl_med[ye] # le wl per questo anno
    allfigs = []
    filename = cart_out + 'fig_bestspectra_{}.pdf'.format(ye)
    for latg in latgrid:
        if (ye, latg) in bestsza_data.keys(): # ci sono dati?
            fig = plt.figure()
            plt.title('year {} - lat {}: sza {:6.1f}'.format(ye, latg, bestsza[(ye, latg)]))
            for pixel in bestsza_data[(ye, latg)]:
                plt.plot(wlgrid, pixel.spet, label = pixel.alt)
            plt.legend()
            #fig.savefig(cart_fig + 'spet_{}_{}.pdf'.format(ye, int(latg)))
            allfigs.append(fig)

    lvs.plot_pdfpages(filename, allfigs)


# prova out observ.dat
lvs.write_obs_frompixels(cart_out + 'test_observ.dat', bestsza_data[(2015, 15.0)])
