#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import re
import mgutils as mg, mgutils.constants as co


if __name__ == '__main__':

    if "-h" in argv:
        print ("discovery_year.py usage")
        raise SystemExit

    separate = True


    table = Table.read('amcvn_catalogue.fits')

    if separate:
        ax1 = mg.formatGraph(1, xlabel='Year', grid=False, figsize=(6,4.5))
        ax2 = mg.formatGraph(2, xlabel='Year', grid=False, figsize=(6,4.5))
        axs = [ax1, ax2]
    else:
        fig, axs = mg.formatSubplots((2,1), sharex=True, xlabel='Year', grid=False, figsize=(6,8))


    ref_years = table['Discovery_year']

    ok = table['Confirmed']
    
    bins = np.arange(ref_years.min(), ref_years.max()+3)


    ### Simple plot of number known

    axs[0].set_ylabel('Number of known AM CVn binaries')
    axs[0].hist(ref_years, cumulative=True, bins=bins, log=False, histtype='step', color='r', lw=3, ls='--', label='Candidate')
    axs[0].hist(ref_years[ok], cumulative=True, bins=bins, log=False, histtype='step', color='k', lw=3, label='Confirmed')

    axs[0].set_ylim(0, len(table)*1.05)
    axs[0].legend(loc='upper left')


    ### Divided by disc state and other properties


    direct = table['Disk_state'] == 'direct'
    high = table['Disk_state'] == 'high'
    outburst = table['Disk_state'] == 'outburst'
    low = table['Disk_state'] == 'low'
    eclipsing = table['Period_method'] == 'Eclipses'

    axs[1].set_ylabel('Fraction of class known')
    axs[1].hist(ref_years[ok&direct], cumulative=True, bins=bins, histtype='step', weights=np.ones(np.sum(ok&direct))/np.sum(ok&direct), \
            label='Direct impact', lw=3, ls=(0, (3, 5, 1, 5, 1, 5)))
    axs[1].hist(ref_years[ok&high], cumulative=True, bins=bins, histtype='step', weights=np.ones(np.sum(ok&high))/np.sum(ok&high),\
             label='High state', lw=3, ls='-')
    axs[1].hist(ref_years[ok&outburst], cumulative=True, bins=bins, histtype='step', weights=np.ones(np.sum(ok&outburst))/np.sum(ok&outburst),\
             label='Outbursting', lw=3, ls='--')
    axs[1].hist(ref_years[ok&low], cumulative=True, bins=bins, histtype='step', weights=np.ones(np.sum(ok&low))/np.sum(ok&low),\
             label='Low state', lw=3, ls='-.')
    axs[1].hist(ref_years[ok&eclipsing], cumulative=True, bins=bins, histtype='step', weights=np.ones(np.sum(ok&eclipsing))/np.sum(ok&eclipsing),\
             label='Eclipsing', lw=3, ls=':')

    axs[1].legend(loc='upper left')
    axs[1].set_ylim(0,1)

    final_year = ref_years.max()+1
    axs[0].set_xlim(1965, final_year)
    axs[1].set_xlim(1965, final_year)
    print('X limit set to', axs[1].get_xlim()[1])
    
    ### Annotate some surveys

    axs[0].vlines([1998,2008,2013,2018], axs[0].get_ylim()[0], axs[0].get_ylim()[1], color='grey', alpha=0.5, ls='--')
    axs[1].vlines([1998,2008,2013,2018], axs[1].get_ylim()[0], axs[1].get_ylim()[1], color='grey', alpha=0.5, ls='--')

    axs[0].annotate('SDSS', (1998, axs[0].get_ylim()[1]*0.88), rotation=30, fontsize='small')
    axs[0].annotate('ASAS-SN', (2008, axs[0].get_ylim()[1]*0.85), rotation=30, fontsize='small')
    axs[0].annotate('PTF', (2013, axs[0].get_ylim()[1]*0.85), rotation=30, fontsize='small')
    axs[0].annotate('ZTF', (2018, axs[0].get_ylim()[1]*0.88), rotation=30, fontsize='small')

    if separate:
        plt.figure(1)
        plt.savefig('discovery_year_separate.png')
        plt.savefig('discovery_year_separate.pdf')
        plt.figure(2)
        plt.savefig('discovery_year_separate_by_class.png')
        plt.savefig('discovery_year_separate_by_class.pdf')
    else:
        plt.savefig('discovery_year.png')
        plt.savefig('discovery_year.pdf')

    plt.show()
    print('Done')