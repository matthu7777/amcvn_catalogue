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
        print ("cmd.py usage")
        raise SystemExit


    parallax_sig = 5

    table = Table.read('amcvn_catalogue_secret.fits')
    oi = table['Confirmed'] 
    ok = (table['Gaia_parallax'] / table['Gaia_parallax_err'] > parallax_sig) 


    bprp = table['Gaia_BPmag'] - table['Gaia_RPmag']
    bprp_err = np.sqrt(table['Gaia_BPmag_err']**2 + table['Gaia_RPmag_err']**2)
    
    absG = table['Gaia_Gmag'] - 5 * np.log10(table['Distance_mean']/10)
    dist_cont = 5 * (np.log10(table['Distance_84percentile']) - np.log10(table['Distance_16percentile'])) / 2
    absG_err = np.sqrt(table['Gaia_BPmag_err']**2 + dist_cont**2)

    
    ax = mg.formatGraph(xlabel='$BP-RP$', ylabel='$M_G$', grid=False)

    ax.errorbar(bprp[~oi&ok], absG[~oi&ok], absG_err[~oi&ok], bprp_err[~oi&ok], fmt='C3o', label='Candidate')
    ax.errorbar(bprp[oi&ok], absG[oi&ok], absG_err[oi&ok], bprp_err[oi&ok], fmt='C0s', label='Confirmed')

    ax.invert_yaxis()


    ### Background

    with fits.open('/home/matthew/work/general_data/hrd/hrd_wd_b_gt_30deg.fits') as hdul:
        bg_wd = hdul[1].data
    wd_bprp = bg_wd['bp_rp']
    wd_absG = bg_wd['phot_g_mean_mag'] - 5 * np.log10(100/bg_wd['parallax'])
    wd_ok = np.isfinite(wd_bprp) & np.isfinite(wd_absG) & (bg_wd['parallax_over_error'] > 5)
    mg.kernel_plot(wd_bprp[wd_ok], wd_absG[wd_ok], '.background_wds.fits', cmap='Greys')

    with fits.open('/home/matthew/work/general_data/hrd/hrd_background_edr3_100pc.fits') as hdul:
        bg_ms = hdul[1].data
    thin = 1
    ms_bprp = (bg_ms['bp_rp'])[::thin]
    ms_absG = (bg_ms['phot_g_mean_mag'] - 5 * np.log10(100/bg_ms['parallax']))[::thin]
    ms_ok = np.isfinite(ms_bprp) & np.isfinite(ms_absG) & (bg_ms['parallax'] / bg_ms['parallax_error'] > 10)[::thin] & (ms_absG < 5 + 5*ms_bprp) & \
                (ms_bprp < 2.5)
    mg.kernel_plot(ms_bprp[ms_ok], ms_absG[ms_ok], '.background_mss.fits', cmap='Greys')

    ax.set_xlim(-0.8,2)
    ax.set_ylim(16,3.9)

    ax.legend(loc='lower right', framealpha=1)

    plt.savefig('cmd.png')
    plt.savefig('cmd.pdf')

    plt.show()