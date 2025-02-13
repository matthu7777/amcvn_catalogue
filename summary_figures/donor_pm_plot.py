#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import mgutils as mg, mgutils.constants as co

# formats = {1:'method',2:'class'}
formats = {1:'class',2:'class'}

def cvMassRadius(period, n=200, minmasslog=-2.5, maxmasslog=-0.1):
    """ Returns masses and radii corresponding to a CV or AM CVn of a certain period
    """
    mass = np.logspace(minmasslog, maxmasslog, num=n)
    rad = np.power((0.01**3 * (mass / 0.1) * np.power((period / 101.), 2)), 1./3)
    return mass, rad


def cv_period(m,r):
    return 101 * np.sqrt(r**3 / 0.01**3 / m * 0.1) / 60




if __name__ == '__main__':

    if "-h" in argv:
        print ("donor_pm_plot.py usage")
        raise SystemExit

    fig1, ax1 = mg.formatGraph(1, ylabel='Donor mass [$M_\\odot$]', xlabel='Period [days]', grid=False, figsize=(8,6.5), returnFig=True)
    fig2, ax2 = mg.formatGraph(2, ylabel='Donor mass [$M_\\odot$]', xlabel='Period [days]', grid=False, figsize=(8,6.5), returnFig=True)



    ## Wong
    for j,fname in enumerate([
                'wong_1e7.txt',
                'wong_he.txt',
                'wong_he_adiabatic.txt',
                'wong_2e7.txt',
                'wong_2e7_adiabatic.txt',
                'wong_3e7.txt',
                'wong_3e7_adiabatic.txt',
                'wong_5e6.txt',
                ]):
        logm,logr = np.loadtxt(fname, unpack=True, skiprows=1)
        label = "WD donor (Wong+ 2021)" if fname=='wong_1e7.txt' else \
                "He star donor (Wong+ 2021)" if fname=='wong_he.txt' else \
                "No cooling (Wong+ 2021)" if fname=='wong_2e7_adiabatic.txt' else \
                None
        fmt = 'k:' if 'adiabatic' in fname else 'k--' if fname=='wong_he.txt' else 'k-'
        m = 10**logm
        p = cv_period(10**logm, 10**logr)
        ax1.plot(p,m,fmt, label=label, lw=3)

    # ## Bauer
    # for j,fname in enumerate([
    #               'bauer37.txt',
    #               'bauer47.txt',
    #               ]):
    #     m,logr = np.loadtxt(fname, unpack=True, usecols=(1,5))
    #     p = cv_period(m, 10**logr)
    #     label = 'He star donor (Bauer+ 2021)' if j==0 else None
    #     ax1.plot(p,m,'k-.', lw=3, label=label)
        
    ## Rajamuthukumar
    for j,fname in enumerate([
                  'rajamuthukumar_nova.npy',
                  'rajamuthukumar_supernova.npy',
                  ]):
        m,r = np.load(fname).T
        p = cv_period(m, r)
        label = 'He star donor (Rajamuthukumar+ 2024)' if j==0 else None
        ax1.plot(p,m,'-.', lw=3, label=label, color='darkslategrey')
    
    ## Sarkar
    for j,fname in enumerate([
                  'sarkar15.txt',
                  'sarkar17.txt',
                  'sarkar19.txt',
                  'sarkar21.txt',
                  ]):
        m,p = np.loadtxt(fname, unpack=True)
        ok = m > 0.01
        label = 'Evolved CV, DD (Sarkar+ 2023)' if j==0 else None
        ax2.plot(p[ok],m[ok],'k-', lw=3, label=label)
    

    ## Belloni
    for j,fname in enumerate([
                  'belloni_1_1.dat',
                  'belloni_1_2.dat',
                  'belloni_1_3.dat',
                  'belloni_6_1.dat',
                  'belloni_6_2.dat',
                  'belloni_6_3.dat',
                  ]):
        m,logr = np.loadtxt(fname, unpack=True)
        p = cv_period(m, 10**logr)
        label = 'Evolved CV, CARB (Belloni+ 2023)' if j==0 else None
        ax2.plot(p,m,'k:', lw=3, label=label)
    


    #### Read table

    table = Table.read('amcvn_catalogue.fits')

    p = np.array(table['Period'])
    q = np.array(table['Q'])
    m1 = np.array(table['M1'])
    m2 = np.array(table['M2'])
    q_err = np.array(table['Q_err'])
    m1_err = np.array(table['M1_err'])
    m2_err = np.array(table['M2_err'])
    method = table['Q_method']


    eclipse = (method=='Eclipses')
    spec = (method=='Spectroscopy')
    superhumps = (method=='Superhumps') | (method=='Stage A Superhumps')

    amcvn = table['AM_CVn']
    hecv = table['He_CV']
    sdb = table['sdB_donor']
    # bd = table['BD_donor']
    # other = ( ~amcvn & ~hecv & ~sdb & ~bd )
    other = ( ~amcvn & ~hecv & ~sdb)

    assume_m1 = np.isnan(m1)
    m1[assume_m1] = 0.83
    m1_err[assume_m1] = 0.17
    m2[assume_m1] = (m1 * q)[assume_m1]
    m2_err[assume_m1] = (np.sqrt((m1_err/m1)**2 + (q_err/q)**2) * m2)[assume_m1]


    #### Plot stars

    for j, ax in enumerate([ax1, ax2]):
        if formats[j+1] == 'method':
            ax.errorbar(p[superhumps], m2[superhumps], m2_err[superhumps], color='grey', marker='o', ls='none', label='Superhumps', markersize=10)
            ax.errorbar(p[spec], m2[spec], m2_err[spec], color='g', marker='s', ls='none', label='Spectroscopic', markersize=10)
            ax.errorbar(p[eclipse], m2[eclipse], m2_err[eclipse], color='C0', marker='*', ls='none', label='Eclipsing', markersize=16)
        
        if formats[j+1] == 'class':
            ax.errorbar(p[amcvn], m2[amcvn], m2_err[amcvn], color='C0', marker='o', ls='none', label='Classical AM CVn', markersize=10)
            ax.errorbar(p[hecv], m2[hecv], m2_err[hecv], color='C4', marker='D', ls='none', label='He CV', markersize=10)
            ax.errorbar(p[sdb], m2[sdb], m2_err[sdb], color='C2', marker='^', ls='none', label='sdB donor', markersize=10)
            # ax.errorbar(p[bd], m2[bd], m2_err[bd], color='k', marker='v', ls='none', label='BD donor', markersize=10)
            ax.errorbar(p[other], m2[other], m2_err[other], color='C3', marker='s', ls='none', label='Other UCB', markersize=10)


    for ax in [ax1, ax2]:
        # ax.legend(loc='lower left', ncol=2)

        ax.set_xlim(0, 80)
        ax.set_ylim(0.006, 0.6)

        ax.set_yscale('log')
        # yticks = [0.003, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5]
        yticks = [0.01, 0.03, 0.05, 0.1, 0.3, 0.5]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)


    fig1.savefig('donor_pm_plot_1.pdf')
    fig2.savefig('donor_pm_plot_2.pdf')
    plt.show()

