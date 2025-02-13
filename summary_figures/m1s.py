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
        print ("m1s.py usage")
        raise SystemExit


    ax = mg.formatGraph(xlabel='Accretor mass [$M_\\odot$]', grid=False, figsize=(6.4,5.5))


    table = Table.read('amcvn_catalogue.fits')


    j = 0


    amcvn = table['AM_CVn']
    hecv = table['He_CV']
    sdb = table['sdB_donor']
    bd = table['BD_donor']
    other = ( ~amcvn & ~hecv & ~sdb & ~bd )

    for k,row in enumerate(table):
        if np.isfinite(row['M1']) and np.isfinite(row['M1_err']):

            fmt = 'C0o' if amcvn[k] else \
                  'C4D' if hecv[k] else \
                  'C2^' if sdb[k] else \
                  'kv' if bd[k] else \
                  'C3s'
            
            print(row['Q_method'])

            plt.errorbar(row['M1'], j+0.2, 0, row['M1_err'], fmt=fmt)

            plt.annotate(row['Name'], (0.02,j), color=fmt[:-1])

            j -= 1

    # plt.fill_betweenx(np.linspace(j-2, 1, 20), 0.83-0.17, 0.83+0.17, color='lightgrey')
    plt.fill_betweenx(np.linspace(j-2, 2, 20), 0.95, 1.05, color='lightgrey')

    # plt.plot([1,1], [j-2,1], 'k--')


    plt.xlim(0, 1.15)
    plt.ylim(j-0.6, 2.5)
    plt.yticks([])


    ### CV masses

    plt.annotate("CVs (Pala+ 2019)", (0.83, j-0.3), color='k', fontsize=12, horizontalalignment='center')
    plt.annotate("SN Ia double-detonation\nrange (approx)", (1.0, (j+1)/2), color='k', fontsize=12, horizontalalignment='center', verticalalignment='center', rotation=90)

    plt.arrow(0.83, j+0.5, -0.17, 0, length_includes_head=True, head_width=0.3, head_length=0.03, color='k', overhang=0.8)
    plt.arrow(0.83, j+0.5, 0.17, 0, length_includes_head=True, head_width=0.3, head_length=0.03, color='k', overhang=0.8)

    ### WD masses

    plt.annotate("Isolated WDs\n(O'Brien+ 2024)", (0.61, 0.5), color='k', fontsize=12, horizontalalignment='center')

    plt.arrow(0.61, 2, -0.05, 0, length_includes_head=True, head_width=0.3, head_length=0.03, color='k', overhang=0.8)
    plt.arrow(0.61, 2, 0.05, 0, length_includes_head=True, head_width=0.3, head_length=0.03, color='k', overhang=0.8)

    print(-j, "targets have measured M1")

    plt.savefig("m1s.pdf")

    plt.show()