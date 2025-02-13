#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import copy
import mgutils as mg, mgutils.constants as co





if __name__ == '__main__':

    if "-h" in argv:
        print ("sky_location.py usage")
        raise SystemExit


    table = Table.read('amcvn_catalogue.fits')

    ra = table['RA']*u.deg
    dec = table['Dec']*u.deg
    close = table['Distance'] < 300
    confirmed = table['Confirmed']
    ra[ra > 180*u.deg] -= 360*u.deg
    # ra,dec = wrap_coords(ra, dec)

    plt.rcParams['text.usetex'] = True
    fig = plt.figure(1, figsize=(8,4.5))
    ax_radec_coords = fig.add_subplot(111, projection="mollweide")
    ax_radec_coords.grid(True)


    # ax_radec_coords.scatter(ra[~close&~confirmed].to(u.rad), dec[~close&~confirmed].to(u.rad), color='C3', marker='o')
    # ax_radec_coords.scatter(ra[close].to(u.rad), dec[close].to(u.rad), color='C1', marker='D')
    # ax_radec_coords.scatter(ra[~close&confirmed].to(u.rad), dec[~close&confirmed].to(u.rad), color='C0', marker='s')
    ax_radec_coords.scatter(ra[~confirmed].to(u.rad), dec[~confirmed].to(u.rad), color='C3', marker='o', label='Candidates')
    ax_radec_coords.scatter(ra[confirmed].to(u.rad), dec[confirmed].to(u.rad), color='C0', marker='s', label='Confirmed')

    print(np.sum(close&~confirmed))


    ### Add Galactic plane
    l = np.linspace(-180,180,300)*u.deg
    b = np.zeros(len(l))*u.deg
    plane = coordinates.SkyCoord(l=l, b=b, frame='galactic')
    plane_ra = copy.deepcopy(plane.transform_to('icrs').ra)
    plane_ra[plane_ra > 180*u.deg] = plane_ra[plane_ra > 180*u.deg] - 360*u.deg
    plane_ra = plane_ra.to(u.rad)
    plane_dec = plane.icrs.dec.to(u.rad)
    plane_dec = plane_dec[np.argsort(plane_ra)]
    plane_ra = plane_ra[np.argsort(plane_ra)]
    
    ax_radec_coords.plot(plane_ra[plane_ra < 180*u.deg], plane_dec[plane_ra < 180*u.deg], 'k-')
    ax_radec_coords.plot(plane_ra[plane_ra > 180*u.deg] - 360*u.deg, plane_dec[plane_ra > 180*u.deg], 'k-')


    print("Northern hemisphere:", np.sum(close&(dec>0)))
    print("Southern hemisphere:", np.sum(close&(dec<0)))

    plt.legend(loc='lower right', framealpha=1)

    plt.savefig('sky_location.pdf')
    plt.savefig('sky_location.png')
    plt.show()
    



