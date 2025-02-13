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
        print ("colour-colour.py usage")
        raise SystemExit


    xlim = (-0.56, 0.18)
    ylim = (0.8, -0.8)




    ax = mg.formatGraph(xlabel='\\textit{g\'-r\'}', ylabel='\\textit{u\'-g\'}', grid=False)

    ### Background

    with fits.open('gentile_fusillo_gedr3_wd_catalogue_clean.fits') as hdul:
        back = hdul[1].data
    ok = np.isfinite(back['umag']) & np.isfinite(back['gmag']) & np.isfinite(back['rmag'])

    # mg.kernel_plot(back['gmag'][ok]-back['rmag'][ok], back['umag'][ok]-back['gmag'][ok], \
    #                             '.wd_ugr.fits', ax=ax, cmap='Greys', thin=3)

    plt.hist2d(back['gmag'][ok]-back['rmag'][ok], back['umag'][ok]-back['gmag'][ok], bins=50, \
                    range=(xlim, ylim[::-1]), cmap='Greys')




    ### UCBs


    table = Table.read('x_sdss.fits')

    ok = table['Confirmed'] & (table['gmag'] < 20.5)
    classical = (~table['Has_optical_hydrogen'] | (table['Name'] == 'HM Cnc')) & ok
    sdb = (table['Name'] == 'ZTF J2130+4420') | (table['Name'] == 'ZTF J2055+4651') & ok
    hecv = ~classical & ~sdb & ok


    ug = table['umag'] - table['gmag']
    gr = table['gmag'] - table['rmag']
    ri = table['rmag'] - table['imag']
    uge = np.sqrt(table['e_umag']**2 + table['e_gmag']**2)
    gre = np.sqrt(table['e_gmag']**2 + table['e_rmag']**2)
    rie = np.sqrt(table['e_rmag']**2 + table['e_imag']**2)

    plt.errorbar(gr[classical],ug[classical],uge[classical],gre[classical], 'C0o', label="AM CVn")
    plt.errorbar(gr[hecv],ug[hecv],uge[hecv],gre[hecv], 'mD', label="He CV")


    # plt.errorbar(0.1, -0.6, np.nanmedian(np.sqrt(table['e_umag'][classical]**2 + table['e_gmag'][classical]**2)), \
    #                         np.nanmedian(np.sqrt(table['e_gmag'][classical]**2 + table['e_rmag'][classical]**2)), 'C0-')

    selected_roelofs = ug < np.minimum(gr*1.35 + 0.32, 0.14) - uge
    selected_roelofs &= (-0.42 + gre < gr) & (gr < 0.02 - gre)
    selected_roelofs &= (-0.33 + rie < ri) & (ri < 0.03 - rie)
    selected_roelofs &= table['gmag'] < 20.5

    selected_carter14 = ug < np.minimum(gr*2.83 + 1.05, 0.35) - uge
    selected_carter14 &= (-0.48 + gre < gr) & (gr < 0.02 - gre)
    selected_carter14 &= (-0.35 + rie < ri) & (ri < 0.03 - rie)
    selected_carter14 &= table['gmag'] < 20.5

    print(f"{np.sum(selected_roelofs & classical)} out of {np.sum(classical & (table['gmag'] < 20.5))} classical AM CVns with colours and within mag limit are within selection box")
    print(f"{np.sum(selected_carter14 & classical)} out of {np.sum(classical & (table['gmag'] < 20.5))} classical AM CVns with colours and within mag limit are within selection box")
    
    
    






    ### Roelofs/Carter selections

    gr_box = np.linspace(-0.42, 0.02, 100)
    ug_box = np.minimum(1.35*gr_box + 0.32, 0.14)

    gr_boxc = np.linspace(-0.48, 0.02, 100)
    ug_boxc = np.minimum(2.83*gr_boxc + 1.05, 0.35)


    plt.plot(gr_box, ug_box, 'k:', lw=2)
    plt.plot([-0.42,-0.42], [ug_box.min(), -2], 'k:', lw=2)
    plt.plot([0.02,0.02], [0.14, -2], 'k:', lw=2)

    plt.plot(gr_boxc, ug_boxc, 'k--', lw=2)
    plt.plot([-0.48,-0.48], [ug_boxc.min(), -2], 'k--', lw=2)
    plt.plot([0.02,0.02], [0.35, -2], 'k--', lw=2)

    plt.gca().invert_yaxis()

    plt.xlim(*xlim)
    plt.ylim(*ylim)
    

    plt.legend(loc='upper left', framealpha=1)

    # for row in table:
    #     plt.annotate(row['Name'], (row['gmag']-row['rmag'],row['umag']-row['gmag']))




    plt.savefig("colour-colour.pdf")
    plt.show()