#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import mgutils as mg, mgutils.constants as co



def eggletonMR(m1,m2,period):
    """ Returns donor radii calculated using the Eggleton Roche lobe radius formula. Period in seconds. Wrapper around eggletonRoche
    """
    G = 6.67408e-11
    mSun = 1.98855e30
    rSun = 695.7e6
    a = (G*mSun*(m1+m2)*(period*60/(2.*np.pi))**2)**(1./3.)/rSun
    q = m2/m1
    return (a*0.49*q**(2./3)) / (0.6*q**(2./3)+np.log(1.+q**(1./3)))


def cvMassRadius(period, n=200, minmasslog=-2.5, maxmasslog=-0.1):
    """ Returns masses and radii corresponding to a CV or AM CVn of a certain period. Wrapper around cvRad.
    """
    mass = np.logspace(minmasslog, maxmasslog, num=n)
    rad = np.power((0.01**3 * (mass / 0.1) * np.power((period / 101.), 2)), 1./3)
    return mass, rad


def plot_stars(m1, m2, m2_err, p, colour='k', markersize=12, marker='o', label=None):
    plt.plot(np.log10([m2-m2_err,m2+m2_err]), np.log10(eggletonMR(m1,np.array([m2-m2_err,m2+m2_err]),p)), '-', linewidth=3, color=colour)  #Error bar line
    plt.plot(np.log10(m2), np.log10(eggletonMR(m1,m2,p)), color=colour, markersize=markersize, label=label, marker=marker, ls='')    # Marker
    return None


    
if __name__ == '__main__':

    if "-h" in argv:
        print ("donor_mr_plot.py usage")
        raise SystemExit

    ax = mg.formatGraph(xlabel='Log donor mass [$M_\\odot$]', ylabel='Log donor radius [$R_\\odot$]', grid=False, figsize=(8,6))


    #### Plot period diagonals
    for p in [10,20,30,40,50,60,70]:
        m, r = cvMassRadius(p*60.)
        plt.plot(np.log10(m),np.log10(r),color='grey',ls=':')
    plt.annotate("10", (-1.23,-1.58),color='grey',)
    plt.annotate("20", (-1.62,-1.51),color='grey',)
    plt.annotate("30", (-1.85,-1.47),color='grey',)
    plt.annotate("40", (-2.01,-1.44),color='grey',)
    plt.annotate("50", (-2.16,-1.42),color='grey',)
    plt.annotate("60", (-2.25,-1.40),color='grey',)
    plt.annotate("70", (-2.35,-1.38),color='grey',)


    #### Plot models

    ## Deloye
    for j,fname in enumerate([
                #   'deloye_b1.dat',
                #   'deloye_g1.dat',
                #   'deloye_r1.dat',
                #   'deloye_y1.dat',
                  'deloye_b2.dat',
                  'deloye_g2.dat',
                  'deloye_r2.dat',
                  'deloye_y2.dat',
                  ]):
        x,y = np.loadtxt(fname, unpack=True)
        label = 'WD donor (Deloye+ 2007)' if j==0 else None
        plt.plot(x,y,'-', lw=3, label=label, color='grey')

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
    #     plt.plot(x,y,'k--', lw=3, label=label)

    # ## Goliasch
    # for j,fname in enumerate([
    #               'goliasch.txt',
    #               'Goliasch0_85.dat',
    #               ]):
    #     y,x = np.loadtxt(fname, unpack=True)
    #     label = 'Evolved CV (Goliasch \\& Nelson 2015)' if j==0 else None
    #     plt.plot(x,y,'k:', lw=3, label=label)
    






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

    classical = ~table['Has_optical_hydrogen'] | (table['Name'] == 'HM Cnc')

    eclipse = (method=='Eclipses')
    spec = (method=='Spectroscopy')
    superhumps = (method=='Superhumps') | (method=='Stage A Superhumps')

    sdb = (table['Name'] == 'ZTF J2130+4420') | (table['Name'] == 'ZTF J2055+4651') 

    hecv = ~classical & ~sdb




    assume_m1 = np.isnan(m1)
    m1[assume_m1] = 0.8
    m1_err[assume_m1] = 0.1
    m2[assume_m1] = (m1 * q)[assume_m1]
    m2_err[assume_m1] = (np.sqrt((m1_err/m1)**2 + (q_err/q)**2) * m2)[assume_m1]


    #### Plot stars

    plot_stars(m1[classical], m2[classical], m2_err[classical], p[classical], \
                                        colour='black', marker='o', label='Classical AM CVn')
    # plot_stars(m1[classical&spec], m2[classical&spec], m2_err[classical&spec], p[classical&spec], \
    #                                     colour='g', marker='s', label='Spectroscopic')
    # plot_stars(m1[classical&eclipse], m2[classical&eclipse], m2_err[classical&eclipse], p[classical&eclipse], \
    #                                     colour='b', marker='*', label='Eclipsing', markersize=18)
    plot_stars(m1[hecv], m2[hecv], m2_err[hecv], p[hecv], \
                                        colour='m', marker='D', label='He CV')
    plot_stars(m1[sdb], m2[sdb], m2_err[sdb], p[sdb], \
                                        colour='r', marker='s', label='sdB donor')
    
    plt.legend(loc='upper left', framealpha=0, frameon=False)
    plt.xlim(-2.4, -0.2)
    plt.ylim(-1.6, -0.65)

    # plt.savefig('donor_mr_plot.pdf')
    # plt.savefig('donor_mr_plot.png')
    plt.show()

    