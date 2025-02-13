#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import mgutils as mg, mgutils.constants as co



def convert_v_g(v, bprp):
    g = v - 0.01760 - 0.006860 * bprp - 0.1732 * bprp**2 
    return g



if __name__ == '__main__':

    if "-h" in argv:
        print ("magnitudes.py usage")
        raise SystemExit

    ax = mg.formatGraph(xlabel="Orbital period [days]", ylabel="Absolute $G$ magnitude", grid=False, figsize=(7,6))


    #### Read table

    table = Table.read('amcvn_catalogue.fits')

    porb = table['Period']
    distance = table['Distance']
    absG = table['Gaia_Gmag'] - 5 * np.log10(distance/10)
    absG_err = table['Gaia_Gmag_err']
    dist_cont = 5 * (np.log10(table['Distance_84']) - np.log10(table['Distance_16'])) / 2
    # dist_cont = 5 * np.log10(np.array([mg.functApp(np.divide, 1000,0, p,pe)[1] for p,pe in zip(table['Gaia_parallax'],table['Gaia_parallax_err'])])/distance + 1)
    absG_err = np.sqrt(table['Gaia_Gmag_err']**2 + dist_cont**2)

    ok = np.isfinite(porb) & (absG > 0) & (absG_err > 0) & table['Confirmed'] & (table['Period_method'] != 'Predicted') & \
                (table['Gaia_parallax']/table['Gaia_parallax_err'] > 5)

    amcvn = ok & table['AM_CVn']
    hecv = ok & table['He_CV']
    sdb = ok & table['sdB_donor']
    bd = ok & table['BD_donor']
    other = ok & ( ~amcvn & ~hecv & ~sdb & ~bd )

    direct = ok & (table['Disk_state'] == 'Direct')
    high = ok & (table['Disk_state'] == 'High')
    outburst = ok & (table['Disk_state'] == 'Outburst')
    low = ok & (table['Disk_state'] == 'Low')
    hydrogen = ok & table['Has_optical_hydrogen'] & ~(table['Disk_state'] == 'Direct')



    plt.errorbar(porb[direct], absG[direct], absG_err[direct], fmt='*', markersize=12, color='C5', label="Direct impact")
    plt.errorbar(porb[high], absG[high], absG_err[high], fmt='*', markersize=12, color='C1', label="High state")
    plt.errorbar(porb[amcvn&(low|outburst)], absG[amcvn&(low|outburst)], absG_err[amcvn&(low|outburst)], fmt='o', markersize=6, color='C0', label="Low state/outbursting")
    plt.errorbar(porb[hecv], absG[hecv], absG_err[hecv], fmt='D', markersize=6, color='C4', label="He CV")
    plt.errorbar(porb[sdb], absG[sdb], absG_err[sdb], fmt='^', markersize=6, color='C2', label="sdB donor")
    plt.errorbar(porb[bd], absG[bd], absG_err[bd], fmt='v', markersize=6, color='k', label="BD donor")
    plt.errorbar(porb[other], absG[other], absG_err[other], fmt='s', markersize=6, color='C3', label="Other UCB")


    ## Model

    ok = table['AM_CVn'] & table['Confirmed'] & (table['Gaia_parallax'] / table['Gaia_parallax_err'] > 5)
    median_bprp = np.median(table['Gaia_BPRP'][ok])

    x,v = np.loadtxt("mags_disc.csv", unpack=True, delimiter=',', skiprows=1)
    v = v[np.argsort(x)]
    x = np.sort(x)
    g = convert_v_g(v, median_bprp)
    plt.plot(x, g, 'k:', alpha=0.5, lw=2)

    x,v = np.loadtxt("mags_wd_high_mdot.csv", unpack=True, delimiter=',', skiprows=1)
    v = v[np.argsort(x)]
    x = np.sort(x)
    g = convert_v_g(v, median_bprp)
    plt.plot(x, g, 'k-', alpha=0.5, lw=2)

    x,v = np.loadtxt("mags_wd_low_mdot.csv", unpack=True, delimiter=',', skiprows=1)
    v = v[np.argsort(x)]
    x = np.sort(x)
    g = convert_v_g(v, median_bprp)
    plt.plot(x, g, 'k-', alpha=0.5, lw=2)




    ## Fit to dependence


    oi = amcvn & (low|outburst) & ok

    err_inflated = np.sqrt(absG_err**2 + 0.1**2)

    in_mask = absG[oi] > 7 - 0.08*porb[oi]

    popt, pcov, mask = mg.iterativeFit(mg.slope, porb[oi], absG[oi], err_inflated[oi], 60,3, p0=[0.08,7], verbose=True, inmask=in_mask)

    # plt.errorbar(porb[oi], absG[oi], err_inflated[oi], fmt='k.')
    # plt.errorbar(porb[oi][~mask], absG[oi][~mask], err_inflated[oi][~mask], fmt='rD')



    xplot = np.linspace(25,70,100)
    yplot = mg.slope(xplot, *popt)

    plt.plot(xplot,yplot,'C0--')



    plt.xlim(3,75)
    plt.ylim(18, 4)

    plt.legend(loc='lower left', ncol=2)
    plt.savefig("magnitudes.pdf")
    plt.show()