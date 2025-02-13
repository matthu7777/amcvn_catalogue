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


def eggletonMR(m1,m2,period):
    """ Returns donor radii calculated using the Eggleton Roche lobe radius formula. Period in minutes
    """
    G = 6.67408e-11
    mSun = 1.98855e30
    rSun = 695.7e6
    a = (G*mSun*(m1+m2)*(period*60/(2.*np.pi))**2)**(1./3.)/rSun
    q = m2/m1
    return (a*0.49*q**(2./3)) / (0.6*q**(2./3)+np.log(1.+q**(1./3)))


def cvMassRadius(period, n=200, minmasslog=-2.5, maxmasslog=-0.1):
    """ Returns masses and radii corresponding to a CV or AM CVn of a certain period
    """
    mass = np.logspace(minmasslog, maxmasslog, num=n)
    rad = np.power((0.01**3 * (mass / 0.1) * np.power((period / 101.), 2)), 1./3)
    return mass, rad


def plot_stars(m1, m2, m2_err, p, colour='k', markersize=10, marker='o', label=None, ax=None):
    if ax is None:
        ax = plt.gca()
    ax.plot(np.log10([m2-m2_err,m2+m2_err]), np.log10(eggletonMR(m1,np.array([m2-m2_err,m2+m2_err]),p)), '-', linewidth=3, color=colour)  #Error bar line
    ax.plot(np.log10(m2), np.log10(eggletonMR(m1,m2,p)), color=colour, markersize=markersize, label=label, marker=marker, ls='')    # Marker
    return None


    
if __name__ == '__main__':

    if "-h" in argv:
        print ("donor_mr_plot.py usage")
        raise SystemExit

    fig1, ax1 = mg.formatGraph(1, xlabel='Donor mass [$M_\\odot$]', ylabel='Donor radius [$R_\\odot$]', grid=False, figsize=(8,6.5), returnFig=True)
    fig2, ax2 = mg.formatGraph(2, xlabel='Donor mass [$M_\\odot$]', ylabel='Donor radius [$R_\\odot$]', grid=False, figsize=(8,6.5), returnFig=True)




    #### Plot period diagonals
    for ax in [ax1, ax2]:
        for p in [10,20,30,40,50,60,70]:
            m, r = cvMassRadius(p*60.)
            ax.plot(np.log10(m),np.log10(r), color='grey', ls=':')
            ax.annotate(f"{p} min", (-0.2-0.1, np.log10(eggletonMR(1, 10**-0.2, p))-0.02), \
                        color='grey', rotation=27)


    #### Plot models

    # ## Deloye
    # for j,fname in enumerate([
    #               'deloye_b1.dat',
    #               'deloye_g1.dat',
    #               'deloye_r1.dat',
    #               'deloye_y1.dat',
    #               'deloye_b2.dat',
    #               'deloye_g2.dat',
    #               'deloye_r2.dat',
    #               'deloye_y2.dat',
    #               ]):
    #     x,y = np.loadtxt(fname, unpack=True)
    #     label = 'WD donor (Deloye+ 2007)' if j==0 else None
    #     plt.plot(x,y,'k-', lw=3, label=label)

    # ## Yungelson
    # for j,fname in enumerate([
    #               'yung_black.csv',
    #               'yung_blue.csv',
    #               'yung_cyan.csv',
    #               'yung_green.csv',
    #               'yung_red.csv',
    #               ]):
    #     x,y = np.loadtxt(fname, unpack=True, delimiter=',')
    #     label = 'He star donor (Yungelson+ 2008)' if j==0 else None
    #     ax1.plot(x,y,'k--', lw=3, label=label)

    # ## Goliasch
    # for j,fname in enumerate([
    #               'goliasch.txt',
    #               'Goliasch0_85.dat',
    #               ]):
    #     y,x = np.loadtxt(fname, unpack=True)
    #     label = 'Evolved CV, classic (Goliasch+ 2015)' if j==0 else None
    #     ax2.plot(x,y,'k--', lw=3, label=label)
    

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
        x,y = np.loadtxt(fname, unpack=True, skiprows=1)
        label = "WD donor (Wong+ 2021)" if fname=='wong_1e7.txt' else \
                "He star donor (Wong+ 2021)" if fname=='wong_he.txt' else \
                "No cooling (Wong+ 2021)" if fname=='wong_2e7_adiabatic.txt' else \
                None
        fmt = 'k:' if 'adiabatic' in fname else 'k--' if fname=='wong_he.txt' else 'k-'
        ax1.plot(x,y,fmt, label=label, lw=3)


    ## Belloni
    for j,fname in enumerate([
                  'belloni_1_1.dat',
                  'belloni_1_2.dat',
                  'belloni_1_3.dat',
                  'belloni_6_1.dat',
                  'belloni_6_2.dat',
                  'belloni_6_3.dat',
                  ]):
        x10,y = np.loadtxt(fname, unpack=True)
        x = np.log10(x10)
        label = 'Evolved CV, CARB (Belloni+ 2023)' if j==0 else None
        oi = 10**y < 0.17
        ax2.plot(x[oi],y[oi],'k:', lw=3, label=label)
    
    # ## Bauer
    # for j,fname in enumerate([
    #               'bauer37.txt',
    #               'bauer47.txt',
    #               ]):
    #     x10,y = np.loadtxt(fname, unpack=True, usecols=(1,5))
    #     x = np.log10(x10)
    #     label = 'He star donor (Bauer+ 2021)' if j==0 else None
    #     ax1.plot(x,y,'k-.', lw=3, label=label)

    ## Rajamuthukumar
    for j,fname in enumerate([
                  'rajamuthukumar_nova.npy',
                  'rajamuthukumar_supernova.npy',
                  ]):
        x10,y10 = np.load(fname).T
        x = np.log10(x10)
        y = np.log10(y10)
        # ok = m > 0.01
        label = 'He star donor (Rajamuthukumar+ 2024)' if j==0 else None
        ax1.plot(x,y,'-.', lw=3, label=label, color='darkslategrey')
    

    ## Sarkar
    for j,fname in enumerate([
                  'sarkar15.txt',
                  'sarkar17.txt',
                  'sarkar19.txt',
                  'sarkar21.txt',
                  ]):
        m,p = np.loadtxt(fname, unpack=True)
        ok = m > 0.01
        x = np.log10(m)[ok]
        y = np.log10(eggletonMR(0.8, m, p)[ok])
        label = 'Evolved CV, DD (Sarkar+ 2023)' if j==0 else None
        ax2.plot(x,y,'k-', lw=3, label=label)
    


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

    # classical = ~table['Has_optical_hydrogen'] | (table['Name'] == 'HM Cnc')
    # other = (table['Name'] == 'ZTF J2130+4420') | (table['Name'] == 'ZTF J2055+4651') | (table['Name'] == 'ATLAS J1138-5139') | (table['Name'] == 'SDSS J1507+5230')
    # hecv = ~classical & ~other

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
            plot_stars(m1[superhumps], m2[superhumps], m2_err[superhumps], p[superhumps], \
                                                colour='grey', marker='o', label='Superhumps', ax=ax)
            plot_stars(m1[spec], m2[spec], m2_err[spec], p[spec], \
                                                colour='g', marker='s', label='Spectroscopic', ax=ax)
            plot_stars(m1[eclipse], m2[eclipse], m2_err[eclipse], p[eclipse], \
                                                colour='C0', marker='*', label='Eclipsing', markersize=16, ax=ax)
        if formats[j+1] == 'class':
            plot_stars(m1[amcvn], m2[amcvn], m2_err[amcvn], p[amcvn], \
                                                colour='C0', marker='o', label='Classical AM CVn', ax=ax)
            plot_stars(m1[hecv], m2[hecv], m2_err[hecv], p[hecv], \
                                                colour='C4', marker='D', label='He CV', ax=ax)
            plot_stars(m1[sdb], m2[sdb], m2_err[sdb], p[sdb], \
                                                colour='C2', marker='^', label='sdB donor', ax=ax)
            # plot_stars(m1[bd], m2[bd], m2_err[bd], p[bd], \
            #                                     colour='k', marker='v', label='BD donor', ax=ax)
            plot_stars(m1[other], m2[other], m2_err[other], p[other], \
                                                colour='C3', marker='s', label='Other UCB', ax=ax)


    
    

    # Formatting
    for ax in [ax1, ax2]:
        ax.legend(loc='upper left', framealpha=0, frameon=False, ncols=1)
        ax.set_xlim(-2.2, -0.1)
        ax.set_ylim(-1.6, -0.58)

        xticks = [0.01,0.02,0.05,0.1,0.2,0.5]
        ax.set_xticks(np.log10(xticks))
        ax.set_xticklabels(xticks)

        yticks = [0.03,0.05,0.1,0.2]
        ax.set_yticks(np.log10(yticks))
        ax.set_yticklabels(yticks)

        ax.minorticks_off()


    fig1.savefig('donor_mr_plot_1.pdf')
    fig2.savefig('donor_mr_plot_2.pdf')
    plt.show()

    