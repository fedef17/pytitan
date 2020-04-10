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

cart_in = '/home/federico/TITANO/VIMS_data/NEW_COLL_HCN-CH4-C2H2_sza110_limb10/'

cart_out = '/home/federico/TITANO/prescreening_2013-19/'
if not os.path.exists(cart_out):
    os.mkdir(cart_out)

pixs = io.readsav(cart_in + 'PIXs_HCN-CH4-C2H2_sza110_limb10.sav')

datasav = pixs['comppix']
allattrs = datasav.dtype.names # questa è la lista di tutti gli attributi del dataset
print('all attributes: ', allattrs)
