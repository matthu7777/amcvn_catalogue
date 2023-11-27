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
        print ("periods.py usage")
        raise SystemExit


    table = Table.read('amcvn_catalogue_secret.fits')
    ok = table['Confirmed']


    direct = ok & (table['Disk_state'] == 'direct')
    high = ok & (table['Disk_state'] == 'high')
    outburst = ok & (table['Disk_state'] == 'outburst')
    low = ok & (table['Disk_state'] == 'low')
    eclipsing = ok & (table['Period_method'] == 'Eclipses')

    period = table['Period']
    period_subsets = [period[direct], period[high], period[low], period[outburst]]

    bins = np.logspace(np.log10(5), np.log10(70), 10)


    ### Plot

    ax = mg.formatGraph(1, xlabel='Orbital period [min]', ylabel='Number', grid=False)

    n1,_,_ = ax.hist(period[direct], bins, histtype='stepfilled', linewidth=2, facecolor='C0', hatch='//', edgecolor='k', label="Direct impact")
    n2,_,_ = ax.hist(period[high], bins, histtype='stepfilled', linewidth=2, facecolor='C1', hatch='\\\\', edgecolor='k', bottom=n1, label="High state")
    n3,_,_ = ax.hist(period[outburst], bins, histtype='stepfilled', linewidth=2, facecolor='C2', hatch='--', edgecolor='k', bottom=n1+n2, label="Outbursting")
    n4,_,_ = ax.hist(period[low], bins, histtype='stepfilled', linewidth=2, facecolor='C3', hatch='+', edgecolor='k', bottom=n1+n2+n3, label="Low state")


    # ax.hist(period_subsets, bins, density=False, histtype='stepfilled', stacked=True, color=['C0', 'C1', 'C3', 'C2'], label=labels, edgecolor='k', hatch=['//','\\','+','--'], fill=True, lw=2)
    # ax.hist(period[eclipsing], bins, density=False, histtype='step', stacked=False, color=['k'], lw=3, label="Eclipsing")
    ax.set_xscale('log')
    # ax.set_yscale('log')

    xticks = [5, 10, 20, 30, 50, 70]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)

    yticks = [0,5,10,15]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)

    ax.legend(loc='upper left')

    print(np.sum(outburst))
    print(np.sum(low))

    plt.savefig('periods.png')
    plt.savefig('periods.pdf')
    plt.show()