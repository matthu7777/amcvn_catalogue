#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import mgutils as mg, mgutils.constants as co


if __name__ == '__main__':

    if "-h" in argv:
        print ("colour_colour.py usage")
        raise SystemExit



    table = Table.read('amcvn_catalogue.fits')
    hybrid = table['Has_optical_hydrogen']
    confirmed = table['Confirmed']
    oi = confirmed & ~hybrid
    ok = (table['Gaia_parallax'] / table['Gaia_parallax_err'] > parallax_sig) 

    period = table['Period']
    ok_period = oi & ok & np.isfinite(period)


    
