#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import os
import mgutils as mg, mgutils.constants as co


def model_galaxy_dist(fname=None, redo=False, scale_height=120, scale_radius=3000, n=20_000_000, rng=None, thin=1, scale=1):
    """ Following the description in Pretorius 2007, MNRAS 374 1495.
    """
    if rng is None:
        rng = np.random.default_rng()

    if fname is None:
        fname = f'model_galaxy_{scale_height}.dat'
    
    if not redo and os.path.isfile(fname):
        dist,gal_n = np.loadtxt(fname, unpack=True)
    
    else:
        radius = rng.exponential(scale_radius, n)
        # radius = rng.uniform(0,14000, n)
        height = rng.exponential(scale_height, n)
        azimut = rng.uniform(0,np.pi/8.,n)          # To truly follow Pretorius this should be up to pi, but for distances <2000pc pi/8 is fine.
        # azimut = np.zeros(n)

        r_earth = 7620
        dist = np.sqrt(radius**2 + r_earth**2 - 2*radius*r_earth*np.cos(azimut) + height**2)
        dist = np.sort(dist)

        # gal_n,_ = np.histogram(dist, bins=bins, weights=np.ones(len(dist))/len(dist))
        gal_n = np.array(range(len(dist))) / len(dist)

        gal_n *= scale / gal_n[mg.find_nearest(dist, 100, True)]

        np.savetxt(fname, np.column_stack([dist,gal_n]))

    return dist[::thin], gal_n[::thin]


if __name__ == '__main__':

    if "-h" in argv:
        print ("distances.py usage")
        raise SystemExit


    distance_limit = 800
    parallax_sig = 3


    table = Table.read('amcvn_catalogue_secret.fits')
    oi = table['Confirmed'] 
    ok = (table['Gaia_parallax'] / table['Gaia_parallax_err'] > parallax_sig) & \
                    (table['Distance_mean'] < distance_limit)

    distance = table['Distance_mean']
    
    # bins = np.logspace(np.log10(distance.min()), np.log10(distance_limit), 11)
    bins = np.sort(distance[oi&ok])
    bin_log_centres = 10**((np.log10(bins[:-1]) + np.log10(bins[1:])) / 2)


    ax = mg.formatGraph(1, xlabel='Distance [pc]', ylabel='Number', grid=False)

    n2,_,_ = ax.hist(distance[ok], bins, histtype='step', cumulative=True, color='r', lw=2, ls='--', label='Candidate')
    n1,_,_ = ax.hist(distance[oi&ok], bins, histtype='step', cumulative=True, color='k', lw=2, label='Confirmed')


    ### Uniform galactic distribution
    # gal_bins = np.logspace(np.log10(distance.min()), np.log10(distance_limit), 10)
    # gal_bin_log_centres = 10**((np.log10(gal_bins[:-1]) + np.log10(gal_bins[1:])) / 2)
    for scale_height in [300]:
        gal_dist, gal_n = model_galaxy_dist(scale_height=scale_height, redo=False)
        gal_ok = (gal_dist > 76) & (gal_dist < distance_limit)
        # scale_factor = 1 / gal_n[0]
        plt.plot(gal_dist[gal_ok], gal_n[gal_ok], 'k--', label=f"Galactic distribution\n(H={scale_height}pc)")


    ax.set_xscale('log')
    ax.set_xticks([70,100,200,300,400,500,700])
    ax.set_xticklabels([70,100,200,300,400,500,700])
    ax.set_xlim(60, 750)

    ax.set_yscale('log')
    ax.set_yticks([1,3,5,10,30,50])
    ax.set_yticklabels([1,3,5,10,30,50])
    ax.set_ylim(0.7, 10**(np.log10(n2.max())*1.15))

    ax.legend(loc='upper left', framealpha=1)

    plt.savefig('distances.png')
    plt.savefig('distances.pdf')
    plt.show()