#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from scipy.special import erf
import mgutils as mg, mgutils.constants as co


def effective_volume(sigma_min):
    hz = 300

    return 2 * np.pi * hz**3 / sigma_min / hz * np.exp(- 1 / 2 / sigma_min**2 / hz**2) + \
            np.pi**(3/2) * np.sqrt(2) * hz**3 * (1 / sigma_min**2 / hz**2 - 1) * erf(1 / np.sqrt(2) / sigma_min / hz)



if __name__ == '__main__':

    if "-h" in argv:
        print ("temp_volumes.py usage")
        raise SystemExit


    ax = mg.formatGraph(xlabel="Maximum Distance, $d_{\\rm max}$ [pc]", ylabel="Effective Volume, $V_{\\rm eff}$ [pc$^3$]", grid=False)


    dmax = np.logspace(1,5,100)
    sigma_min = 1 / dmax

    veff_vals = effective_volume(sigma_min)

    ax.plot(dmax, veff_vals, 'C0', lw=3)


    ax.plot(dmax, 2362*dmax**2, 'k--', lw=2)
    ax.plot(dmax, 4/3*np.pi*dmax**3, 'k:', lw=2)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(10,1e5)

    plt.show()