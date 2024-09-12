import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.stats import binned_statistic
import os

def broken_pl(r, q1, rbreak, q2):
        pl = np.zeros_like(r)
        pl[r <= rbreak] = r[r<=rbreak]**-q1
        pl[r > rbreak] = rbreak**(q2-q1) * r[r>rbreak]**-q2
        return pl

def get_lines(files):
    enshift = []
    disc_r = []
    weight = []
    for i in range(len(files)):
        with pyfits.open(files[i]) as fits_file: 
            try: 
                enshift.append(np.array(fits_file['ENSHIFT'].data))    # energy shift from each pixel to observer
                disc_r.append(np.array(fits_file['RADIUS'].data))     # radius on the disc seen in each pixel
                #fits_file.info()
                weight.append(np.array(fits_file['WEIGHT'].data))
            except: 
                print('Error with ', files[i])
                fits_file.info()
                print(np.array(fits_file['ENSHIFT'].data))
                print(np.array(fits_file['RADIUS'].data))
                print(np.array(fits_file['WEIGHT'].data))
                 
                return


        enshift[i][np.isnan(enshift[i])] = 0

    disc_emis = []
    for i in range(len(files)):
        disc_emis.append(broken_pl(disc_r[i], 3, 5, 3))

    disc_flux = []

    for i in range(len(files)):
        disc_flux.append(disc_emis[i] * enshift[i]**3 * weight[i])
        disc_flux[i][np.isnan(disc_flux[i])] = 0

    # rest frame line energy
    line_en = 6.4

    # edges of the spectral bins
    bin_edges = np.arange(1,10,0.1)

    # sum up the flux in each bin
    lines = []
    en_edges = []
    energies = []
    for i in range(len(files)):
        line, en_edge, _ = binned_statistic(line_en*enshift[i].flatten(), disc_flux[i].flatten()
                                             , statistic='sum', bins=bin_edges)
        lines.append(line)
        en_edges.append(en_edge)
        energies.append(0.5*(en_edges[i][1:] + en_edges[i][:-1]))
        
    return energies, lines
