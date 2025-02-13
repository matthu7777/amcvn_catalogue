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
    ok = np.array(table['Confirmed']) & (table['Period_method'] != 'Predicted')

    print(len(table))
    print(np.sum(table['Confirmed']))


    direct = ok & (table['AM_CVn']) & (table['Disk_state'] == 'Direct')
    high = ok & (table['AM_CVn']) & (table['Disk_state'] == 'High')
    outburst = ok & (table['AM_CVn']) & (table['Disk_state'] == 'Outburst')
    low = ok & (table['AM_CVn']) & (table['Disk_state'] == 'Low')
    # eclipsing = ok & (table['Period_method'] == 'Eclipses')
    # close = ok & (table['Distance'] < 300)
    # erosita = ok & (table['Has_eROSITA'])
    hecv = ok & (table['He_CV'])
    other = ok & (~high & ~direct & ~outburst & ~low & ~hecv)

    period = table['Period']
    period_subsets = [period[direct], period[high], period[low], period[outburst]]

    bins = np.logspace(np.log10(5), np.log10(70), 10)



    ### Plot

    fig1, ax1 = mg.formatGraph(1, xlabel='Orbital period [min]', ylabel='Number', grid=False, returnFig=True)

    n1,_,_ = ax1.hist(period[direct], bins, histtype='stepfilled', linewidth=2, facecolor='C0', hatch='//', edgecolor='k', label="Direct impact")
    n2,_,_ = ax1.hist(period[high], bins, histtype='stepfilled', linewidth=2, facecolor='C1', hatch='\\\\', edgecolor='k', bottom=n1, label="High state")
    n3,_,_ = ax1.hist(period[outburst], bins, histtype='stepfilled', linewidth=2, facecolor='C2', hatch='--', edgecolor='k', bottom=n1+n2, label="Outbursting")
    n4,_,_ = ax1.hist(period[low], bins, histtype='stepfilled', linewidth=2, facecolor='C3', hatch='+', edgecolor='k', bottom=n1+n2+n3, label="Low state")
    n5,_,_ = ax1.hist(period[hecv], bins, histtype='stepfilled', linewidth=2, facecolor='C4', hatch='||', edgecolor='k', bottom=n1+n2+n3+n4, label="He CV")
    n6,_,_ = ax1.hist(period[other], bins, histtype='stepfilled', linewidth=2, facecolor='k', edgecolor='k', bottom=n1+n2+n3+n4+n5, label="Other")
    ax1.set_xscale('log')

    xticks = [5, 10, 20, 30, 50, 70]
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticks)

    yticks = [0,5,10,15,20,25]
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(yticks)

    ax1.legend(loc='upper left')
    ax1.set_ylim(0,29)

    fig1.savefig('periods1.png')
    fig1.savefig('periods1.pdf')
    plt.show()