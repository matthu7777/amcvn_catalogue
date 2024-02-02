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


    table = Table.read('amcvn_catalogue.fits')
    ok = np.array(table['Confirmed']) & (~np.array(table['Has_optical_hydrogen']) | (table['Disk_state'] == 'direct'))


    direct = ok & (table['Disk_state'] == 'direct')
    high = ok & (table['Disk_state'] == 'high')
    outburst = ok & (table['Disk_state'] == 'outburst')
    low = ok & (table['Disk_state'] == 'low')
    eclipsing = ok & (table['Period_method'] == 'Eclipses')
    close = ok & (table['Distance'] < 300)

    period = table['Period']
    period_subsets = [period[direct], period[high], period[low], period[outburst]]

    bins = np.logspace(np.log10(5), np.log10(70), 10)


    ### Plot

    fig1, ax1 = mg.formatGraph(1, xlabel='Orbital period [min]', ylabel='Number', grid=False, returnFig=True)
    fig2, ax2 = mg.formatGraph(2, xlabel='Orbital period [min]', ylabel='Number', grid=False, returnFig=True)

    n1,_,_ = ax1.hist(period[direct], bins, histtype='stepfilled', linewidth=2, facecolor='C0', hatch='//', edgecolor='k', label="Direct impact")
    n2,_,_ = ax1.hist(period[high], bins, histtype='stepfilled', linewidth=2, facecolor='C1', hatch='\\\\', edgecolor='k', bottom=n1, label="High state")
    n3,_,_ = ax1.hist(period[outburst], bins, histtype='stepfilled', linewidth=2, facecolor='C2', hatch='--', edgecolor='k', bottom=n1+n2, label="Outbursting")
    n4,_,_ = ax1.hist(period[low], bins, histtype='stepfilled', linewidth=2, facecolor='C3', hatch='+', edgecolor='k', bottom=n1+n2+n3, label="Low state")
    ax1.set_xscale('log')


    ax2.hist(period[eclipsing], bins, density=False, histtype='step', stacked=False, color=['k'], lw=3, label="Eclipsing")
    ax2.hist(period[close], bins, density=False, histtype='step', stacked=False, color=['r'], lw=3, label="$D < 300$ pc")
    ax2.set_xscale('log')
    # ax.set_yscale('log')

    for ax in [ax1, ax2]:
        xticks = [5, 10, 20, 30, 50, 70]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)

        yticks = [0,5,10,15]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        ax.legend(loc='upper left')
    ax2.set_ylim(0,9)

    fig1.savefig('periods1.png')
    fig1.savefig('periods1.pdf')
    fig2.savefig('periods2.png')
    fig2.savefig('periods2.pdf')
    plt.show()