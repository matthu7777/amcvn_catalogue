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
    ok = table['Has_eROSITA']
    ra[ra > 180*u.deg] -= 360*u.deg
    # ra,dec = wrap_coords(ra, dec)

    coords = coordinates.SkyCoord(ra=ra, dec=dec, frame='icrs')
    l = np.array(copy.deepcopy(coords.transform_to('galactic').l))
    b = copy.deepcopy(coords.transform_to('galactic').b)
    l[l > 180] = l[l > 180] - 360


    fig1 = plt.figure(1, figsize=(8,4.5))
    ax_radec_coords = fig1.add_subplot(111, projection="mollweide")
    ax_radec_coords.grid(True)

    fig2 = plt.figure(2, figsize=(8,4.5))
    ax_gal_coords = fig2.add_subplot(111, projection="mollweide")
    ax_gal_coords.grid(True)


    ax_radec_coords.scatter(ra[~ok].to(u.rad), dec[~ok].to(u.rad), color='C3', marker='o')
    ax_radec_coords.scatter(ra[ok].to(u.rad), dec[ok].to(u.rad), color='C0', marker='s')

    ax_gal_coords.scatter(l[~ok]*np.pi/180, b[~ok].to(u.rad), color='C3', marker='o')
    ax_gal_coords.scatter(l[ok]*np.pi/180, b[ok].to(u.rad), color='C0', marker='s')



    # ### Add Galactic plane
    # l = np.linspace(-180,180,300)*u.deg
    # b = np.zeros(len(l))*u.deg
    # plane = coordinates.SkyCoord(l=l, b=b, frame='galactic')
    # plane_ra = copy.deepcopy(plane.transform_to('icrs').ra)
    # plane_ra[plane_ra > 180*u.deg] = plane_ra[plane_ra > 180*u.deg] - 360*u.deg
    # plane_ra = plane_ra.to(u.rad)
    # plane_dec = plane.icrs.dec.to(u.rad)
    # plane_dec = plane_dec[np.argsort(plane_ra)]
    # plane_ra = plane_ra[np.argsort(plane_ra)]
    
    # ax_radec_coords.plot(plane_ra[plane_ra < 180*u.deg], plane_dec[plane_ra < 180*u.deg], 'k-')
    # ax_radec_coords.plot(plane_ra[plane_ra > 180*u.deg] - 360*u.deg, plane_dec[plane_ra > 180*u.deg], 'k-')


    for row, l1, b1 in zip(table, l, b):
        if l1 < 0 and ~row['Has_eROSITA']:
            print(row['Name'])




    plt.savefig('sky_location.pdf')
    plt.savefig('sky_location.png')
    plt.show()
    



