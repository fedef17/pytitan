#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import pickle

import numpy as np
import matplotlib.pyplot as plt

import scipy.io as io
from datetime import datetime

############################################################

def trova_spip(ifile, hasha = '#', read_past = False):
    """
    Trova il '#' nei file .dat
    """
    gigi = 'a'
    while gigi != hasha :
        linea = ifile.readline()
        gigi = linea[0]
    else:
        if read_past:
            return linea[1:]
        else:
            return


def read_obs(filename):
    """
    Reads files of VIMS observations. (gbb_2015 format)
    :param filename:
    :return:
    """
    infile = open(filename, 'r')
    trova_spip(infile)
    line = infile.readline()
    cosi = line.split()
    n_freq = int(cosi[0])
    n_limb = int(cosi[1])
    trova_spip(infile)
    dists = []
    while len(dists) < n_limb:
        line = infile.readline()
        dists += list(map(float, line.split()))
    dists = np.array(dists)
    trova_spip(infile)
    alts = []
    while len(alts) < n_limb:
        line = infile.readline()
        alts += list(map(float, line.split()))
    alts = np.array(alts)
    trova_spip(infile)
    data = [line.split() for line in infile]
    data_arr = np.array(data)
    freq = np.array([float(r) for r in data_arr[:, 0]])
    obs = data_arr[:, 1:2*n_limb+2:2]
    obs = obs.astype(float)
    flags = data_arr[:, 2:2*n_limb+2:2]
    flags = flags.astype(int)

    infile.close()
    return n_freq, n_limb, dists, alts, freq, obs, flags


def writevec(ifile, vec, n_per_line, format_str):
    """
    Writes a vector in formatted output to a file.
    :param file: File to write to.
    :param vec: Vector to be written
    :param n_per_line: Number of elements of vector per line
    :param format_str: String format of each number written
    :return: nada
    """
    n = len(vec)
    com = n/n_per_line
    for i in range(com):
        i1 = i*n_per_line
        i2 = i1 + n_per_line
        strin = n_per_line*format_str+'\n'
        ifile.write(strin.format(*vec[i1:i2]))
    nres = n - com * n_per_line
    i1 = com * n_per_line
    if(nres > 0):
        strin = nres*format_str+'\n'
        ifile.write(strin.format(*vec[i1:n]))

    return


def write_obs_frompixels(filename, pixels, wl_range = None):
    """
    Wrap for pixels as in bestsza_data.
    """
    if type(pixels) == list:
        pixels = np.array(pixels)

    n_freq = len(pixels[0].wl)
    n_limb = len(pixels)

    ord_pix = np.argsort([pix.alt for pix in pixels])[::-1]
    dists = [pix.dist for pix in pixels[ord_pix]]
    alts = [pix.alt for pix in pixels[ord_pix]]
    freq = pixels[0].wl

    obs = np.stack([pix.spet for pix in pixels[ord_pix]]).T # Creates a matrix with frequencies in different lines and observations in different columns
    flags = np.stack([pix.spet for pix in pixels[ord_pix]]).T

    write_obs(n_freq, n_limb, dists, alts, freq, obs, flags, filename, wl_range = wl_range)

    return


def write_obs(n_freq, n_limb, dists, alts, freq, obs, flags, filename, old_file = 'None', wl_range = None):
    """
    Writes files of VIMS observations. (gbb_2015 format)
    :param filename:
    :return:
    """

    if wl_range is not None:
        print('Selecting wl range: {}'.format(wl_range))
        okwl = (freq >= wl_range[0]) & (freq <= wl_range[1])
        freq = freq[okwl]
        obs = obs[okwl, :]
        flags = flags[okwl, :]
        n_freq = len(freq)

    infile = open(filename, 'w')
    data = datetime.now()
    infile.write('Modified on: {}\n'.format(data))
    infile.write('Original file: {}\n'.format(old_file))
    infile.write('\n')
    infile.write('Number of spectral points, number of tangent altitudes:\n')
    infile.write('{:1s}\n'.format('#'))
    infile.write('{:12d}{:12d}\n'.format(n_freq,n_limb))
    infile.write('\n')
    infile.write('Altitudes of satellite (km):\n')
    infile.write('{:1s}\n'.format('#'))
    writevec(infile,dists,8,'{:15.4e}')
    #infile.write((8*'{:15.4e}'+'\n').format(dists))
    infile.write('\n')
    infile.write('Tangent altitudes (km): \n')
    infile.write('{:1s}\n'.format('#'))
    writevec(infile,alts,8,'{:10.2f}')
    #infile.write((8*'{:10.2f}'+'\n').format(alts))
    infile.write('\n')
    infile.write('Wavelength (nm), spectral data (W m^-2 nm^-1 sr^-1):\n')
    infile.write('{:1s}\n'.format('#'))

    for fr, ob, fl in zip(freq, obs, flags):
        strin = '{:10.4f}'.format(fr)
        for oo, ff in zip(ob,fl):
            strin = strin + '{:15.4e}{:3d}'.format(oo, int(ff))
        strin = strin + '\n'
        infile.write(strin)

    infile.close()
    return


def read_input_prof_gbb(filename, ptype, n_alt_max =None, n_alt = 151, alt_step = 10.0, n_gas = 86, n_lat = 4, read_bad_names = False):
    """
    Reads input profiles from gbb standard formatted files (in_temp.dat, in_pres.dat, in_vmr_prof.dat).
    Profile order is from surface to TOA.
    type = 'vmr', 'temp', 'pres'
    :return: profiles
    read_bad_names is to read profiles of stuff that is not in the HITRAN list.
    """
    alts = np.linspace(0,(n_alt-1)*alt_step,n_alt)

    infile = open(filename, 'r')

    if(ptype == 'vmr'):
        print(ptype)
        trova_spip(infile)
        trova_spip(infile)
        proftot = []
        mol_names = []
        for i in range(n_gas):
            try:
                lin = infile.readline()
                print(lin)
                num = int(lin.split()[0])
                nome = lin.split()[1]
            except:
                break

            try:
                nome = find_molec_metadata(num, 1)['mol_name']
            except:
                if not read_bad_names:
                    break

            prof = []
            while len(prof) < n_alt:
                line = infile.readline()
                prof += list(map(float, line.split()))
            prof = np.array(prof[::-1])*1.e-6 # SETTING UNITY TO ABSOLUTE FRACTION, NOT PPM
            if n_alt_max is not None and n_alt_max < n_alt:
                prof = prof[:n_alt_max]
            proftot.append(prof)
            try:
                mol_names.append(find_molec_metadata(i+1, 1)['mol_name'])
            except:
                mol_names.append(nome)

            for j in range(n_lat-1): # to skip other latitudes
                prof = []
                while len(prof) < n_alt:
                    line = infile.readline()
                    prof += list(map(float, line.split()))

        proftot = dict(zip(mol_names,proftot))

    if(ptype == 'temp' or ptype == 'pres'):
        print(ptype)
        trova_spip(infile)
        trova_spip(infile)
        prof = []
        while len(prof) < n_alt:
            line = infile.readline()
            prof += list(map(float, line.split()))
        proftot = np.array(prof[::-1])
        if n_alt_max is not None and n_alt_max < n_alt:
            proftot = proftot[:n_alt_max]

    return proftot


def write_input_prof_gbb(filename, prof, ptype, n_alt = 151, alt_step = 10.0, nlat = 4, descr = '', script=__file__):
    """
    Writes input profiles in gbb standard formatted files (in_temp.dat, in_pres.dat, in_vmr_prof.dat)
    Works both with normal vectors and with sbm.AtmProfile objects.
    :return:
    """
    from datetime import datetime

    alts = np.linspace(0,(n_alt-1)*alt_step,n_alt)

    infile = open(filename, 'w')
    data = datetime.now()
    infile.write(descr+'\n')
    infile.write('\n')
    infile.write('Processed through -{}- on: {}\n'.format(script,data))
    infile.write('{:1s}\n'.format('#'))

    n_per_line = 8

    if(ptype == 'vmr'):
        strin = '{:10.3e}'
        infile.write('VMR of molecules (ppmV)\n')
    elif(ptype == 'temp'):
        strin = '{:10.5f}'
        infile.write('Temperature (K)\n')
    elif(ptype == 'pres'):
        strin = '{:11.4e}'
        infile.write('Pressure (hPa)\n')
    else:
        raise ValueError('Type not recognized. Should be one among: {}'.format(['vmr','temp','pres']))

    infile.write('{:1s}\n'.format('#'))
    for i in range(nlat):
        writevec(infile,prof[::-1],n_per_line,strin)

    return


def read_sim_gbb(filename):
    """
    Read sim_*.dat or spet_*.dat files in gbb format.
    :return:
    """
    infile = open(filename, 'r')
    line = infile.readline()
    line = infile.readline()
    alt = line.split()[0]
    trova_spip(infile)

    data = [line.split() for line in infile]
    data_arr = np.array(data)
    freq = np.array([float(r) for r in data_arr[:, 1]])
    obs = np.array([float(r) for r in data_arr[:, 2]])
    sim = np.array([float(r) for r in data_arr[:, 3]])
    err = np.array([float(r) for r in data_arr[:, 4]])
    #flags = data_arr[:, 2:2*n_limb+2:2]
    #flags = flags.astype(np.int)
    infile.close()

    return alt, freq, obs, sim, err


def color_set(n, cmap = 'nipy_spectral', bright_thres = None, full_cb_range = False, only_darker_colors = False):
    """
    Gives a set of n well chosen (hopefully) colors, darker than bright_thres. bright_thres ranges from 0 (darker) to 1 (brighter).

    < full_cb_range > : if True, takes all cb values. If false takes the portion 0.05/0.95.
    """

    if bright_thres is None:
        if only_darker_colors:
            bright_thres = 0.6
        else:
            bright_thres = 1.0

    cmappa = cm.get_cmap(cmap)
    colors = []

    if full_cb_range:
        valori = np.linspace(0.0,1.0,n)
    else:
        valori = np.linspace(0.05,0.95,n)

    for cos in valori:
        colors.append(cmappa(cos))

    for i, (col,val) in enumerate(zip(colors, valori)):
        if color_brightness(col) > bright_thres:
            # Looking for darker color
            col2 = cmappa(val+1.0/(3*n))
            col3 = cmappa(val-1.0/(3*n))
            colori = [col, col2, col3]
            brighti = np.array([color_brightness(co) for co in colori]).argmin()
            colors[i] = colori[brighti]

    return colors


def plotta_sim_VIMS(nomefile, freq, obs, sim, all_sims = None, names = None, err=1.5e-8, title=None, auto = True, xscale=[-1,-1],yscale=[-1,-1],yscale_res=[-1,-1]):
    """
    Plots observed/simulated with residuals and single contributions.
    :param obs: Observed
    :param sim: Simulated total
    :param n_sims: Number of simulated mol. contributions
    :param all_sims: matrix with one sim per row
    :param names: names of each sim
    :return:
    """
    from matplotlib.font_manager import FontProperties
    import matplotlib.gridspec as gridspec

    fontP = FontProperties()
    fontP.set_size('small')

    if all_sims is None:
        n_sims = 2
    else:
        n_sims = all_sims.shape[0]+2

    fig = plt.figure(figsize=(12, 8))

    gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])
    ax1 = plt.subplot(gs[0])
    plt.ylabel('Radiance (W m$^{-2}$ nm$^{-1}$ sr$^{-1}$)')
    plt.title(title)
    ax2 = plt.subplot(gs[1])
    if not auto:
        ax1.set_xlim(xscale)
        ax1.set_ylim(yscale)
        ax2.set_xlim(xscale)
        ax2.set_ylim(yscale_res)

    plt.xlabel('Wavelength (nm)')

#    plt.subplot(211)
    colors = color_set(n_sims)
    ax1.plot(freq,obs,color=colors[0],label='Data',linewidth=1.0)
    ax1.scatter(freq,obs,color=colors[0],linewidth=1.0)
    ax1.errorbar(freq,obs,color=colors[0],yerr=err, linewidth=1.0)

    ax1.plot(freq,sim,color=colors[1],linewidth=3.0,label='Sim')
    i=1
    if all_sims is not None:
        for namu, simu, colu in zip(names, all_sims, colors[2:]):
            ax1.plot(freq, simu, color=colu, linestyle=li, label=namu, linewidth=2.0)
            i +=1
    ax1.grid()
    ax1.legend(loc=1,bbox_to_anchor=(1.05,1.1),fontsize='small',fancybox=1,shadow=1)

    ax2.grid()
    ax2.plot(freq,obs-sim,color='red',linewidth=3.0)
#    ax2.fill_between(freq,err*np.ones(len(freq)),-err*np.ones(len(freq)), facecolor=findcol(12,8)[0], alpha=0.1)
    ax2.plot(freq,err*np.ones(len(freq)),color='black',linestyle='--',linewidth=2.0)
    ax2.plot(freq,-err*np.ones(len(freq)),color='black',linestyle='--',linewidth=2.0)

    fig.savefig(nomefile)
    plt.close(fig)

    return


def plot_pdfpages(filename, figs):
    """
    Saves a list of figures to a pdf file.
    """
    from matplotlib.backends.backend_pdf import PdfPages

    pdf = PdfPages(filename)
    for fig in figs:
        pdf.savefig(fig)
    pdf.close()

    return
