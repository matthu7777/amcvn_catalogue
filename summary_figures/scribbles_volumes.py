#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from scipy.special import erf
from scipy.integrate import quad
import mgutils as mg, mgutils.constants as co


def effective_volume(sigma_min, hz=300):
    """ This is the solution to Eq 15 of Rix2021, solved for me by Anton Biryukov
    """
    return 2 * np.pi * hz**3 / sigma_min / hz * np.exp(- 1 / 2 / sigma_min**2 / hz**2) + \
            np.pi**(3/2) * np.sqrt(2) * hz**3 * (1 / sigma_min**2 / hz**2 - 1) * erf(1 / np.sqrt(2) / sigma_min / hz)


def find_veff_period(periods, dlim=300):
    dmax = np.minimum(dlim, find_dmax(periods))
    return effective_volume(1/dmax)



def find_absmag(p):
    """ This is the fit that is done by the script magnitudes.py. 
    """
    print("WARNING!! UPDATE ABS MAGNITUDE FIT")
    return p*0.06821 + 8.3467



def find_dmax(p):
    maglim = 19.5
    return 10**((maglim - find_absmag(p))/5 + 1)



def period_prob_volume(porb, alpha, norm=1):
    return porb**alpha * find_veff_period(porb) / norm








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



    ax2 = mg.formatGraph(2, xlabel="Period [days]", ylabel="Veff(dmax) / Veff(300pc)", grid=False)
    periods = np.linspace(1e-6,70,100)
    veffs = find_veff_period(periods, 300)




    ax2.plot(periods, veffs/effective_volume(1/300))

    ax2.set_yscale('log')

    print("WARNING!! UPDATE ABS MAGNITUDE FIT")




    quad(period_prob_volume, 0, 70, )
    




    
    plt.show()