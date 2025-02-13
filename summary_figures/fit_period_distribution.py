#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import emcee
import os
from scipy.special import erf
from scipy.integrate import quad
import mgutils as mg, mgutils.constants as co


def effective_volume(sigma_min, hz=300):
    """ This is the solution to Eq 15 of Rix2021, solved for me by Anton Biryukov
    """
    return 2 * np.pi * hz**3 / sigma_min / hz * np.exp(- 1 / 2 / sigma_min**2 / hz**2) + \
            np.pi**(3/2) * np.sqrt(2) * hz**3 * (1 / sigma_min**2 / hz**2 - 1) * erf(1 / np.sqrt(2) / sigma_min / hz)


def find_veff_period(periods, dlim=np.inf, maglim=19.5):
    dmax = np.minimum(dlim, find_dmax(periods, maglim))
    return effective_volume(1/dmax)



def find_absmag(p):
    """ This is the fit that is done by the script magnitudes.py. 
    """
    # print("WARNING!! UPDATE ABS MAGNITUDE FIT")
    return p*0.06821 + 8.3467



def find_dmax(p, maglim):
    return 10**((maglim - find_absmag(p))/5 + 1)


def period_prob_plain(porb, alpha, pmin=0, pmax=70):
    norm = 1/(alpha + 1) * (pmax**(alpha+1) - pmin**(alpha+1))
    return porb**alpha / norm


def period_prob_volume(porb, alpha, maglim=19.5):
    def period_prob_volume_unnormed(porb, alpha, maglim):
        return porb**alpha * find_veff_period(porb, dlim=dlim_vol, maglim=maglim)
    norm = quad(period_prob_volume_unnormed, pmin_vol, pmax, args=(alpha, maglim))[0]
    return period_prob_volume_unnormed(porb, alpha, maglim=maglim) / norm


def period_prob_eclipse(porb, alpha, xi=-0.2, maglim=19.5):
    def period_prob_eclipse_unnormed(porb, alpha, xi, maglim):
        beta = alpha + 2 / (9 * xi - 3)
        return (porb**alpha + porb**beta) * find_veff_period(porb, dlim=np.inf, maglim=maglim)
    norm = quad(period_prob_eclipse_unnormed, pmin_ecl, pmax, args=(alpha, xi, maglim))[0]
    prob = period_prob_eclipse_unnormed(porb, alpha, xi=xi, maglim=maglim) / norm
    prob[porb < pmin_ecl] = np.nan
    return prob
    


def ln_prior(alpha):
    if alpha > 0:
        return 0
    else:
        return -np.inf


def ln_prob_volume(alpha, porb_close):
    lp = ln_prior(alpha)
    ll = np.log(period_prob_volume(porb_close, alpha)).sum()
    if not np.isfinite(ll) or not np.isfinite(lp):
        return -np.inf
    return ll + lp


def ln_prob_eclipse(alpha, porb_eclipse):
    lp = ln_prior(alpha)
    ll = np.log(period_prob_eclipse(porb_eclipse, alpha)).sum()
    if not np.isfinite(ll) or not np.isfinite(lp):
        return -np.inf
    return ll + lp




def mcmc_fit(funct, *args, p0=1, redo=True, sname='.chains.npy', short=False):
    nwalkers, ndim = 16, 1
    nsteps = 2000 if short else 20000
    nburn, nthin = nsteps//2, 100
    pos = p0 + 1e-4*np.random.randn(nwalkers, ndim)

    if not os.path.isfile(sname) or redo:

        sampler = emcee.EnsembleSampler(nwalkers, ndim, funct, args=args)
        sampler.run_mcmc(pos, nsteps, progress=True)

        # samples = sampler.get_chain(discard=nburn, thin=nthin)
        flat_samples = sampler.get_chain(discard=nburn, thin=nthin, flat=True)

        np.save(sname, flat_samples)
    
    else:
        flat_samples = np.load(sname)

    best_fit = flat_samples.mean(axis=0)
    print("Best fit:", best_fit)
    errors = flat_samples.std(axis=0)
    print("Errors:", errors)

    return flat_samples, best_fit, errors





if __name__ == '__main__':

    if "-h" in argv:
        print ("periods.py usage")
        raise SystemExit



    table = Table.read('amcvn_catalogue.fits')
    ok = np.array(table['Confirmed']) & np.array(table['AM_CVn']) & np.isfinite(table['Period'])
    print(len(table))
    print(np.sum(table['Confirmed']))

    dlim_vol = 300
    maglim_vol = 19.5
    maglim_ecl = 19.5

    pmax = 70
    pmin_vol = 1e-5
    pmin_ecl = 25
    pmin_plot = 5
    nbins = 9
    nbins_ecl = 5

    direct = ok & (table['Disk_state'] == 'direct')
    high = ok & (table['Disk_state'] == 'high')
    outburst = ok & (table['Disk_state'] == 'outburst')
    low = ok & (table['Disk_state'] == 'low')
    eclipsing = ok & (table['Has_eclipses']) & (table['Gaia_Gmag'] < maglim_ecl) & (table['Period'] > pmin_ecl)
    close = ok & (table['Distance'] < 300) & (table['Gaia_Gmag'] < maglim_vol) & (table['Gaia_parallax'] / table['Gaia_parallax_err'] > 5)
    erosita = ok & (table['Has_eROSITA'])

    period = table['Period']
    period_subsets = [period[direct], period[high], period[low], period[outburst]]

    bins = np.logspace(np.log10(pmin_plot), np.log10(pmax), nbins+1)
    log_bin_width = (np.log10(pmax) - np.log10(pmin_plot)) / nbins

    bins_ecl = np.logspace(np.log10(pmin_ecl), np.log10(pmax), nbins_ecl+1)
    log_bin_width_ecl = (np.log10(pmax) - np.log10(pmin_ecl)) / nbins_ecl





    ### Plot histogram



    fig, axs = mg.formatSubplots((2,1), xlabel='Orbital period [min]', ylabel='Number', grid=False, sharex=True,)
    ax1, ax2 = axs

    n1,_,_ = ax1.hist(period[close], bins, density=False, histtype='step', stacked=False, color=['r'], lw=3, label=f"$D < {dlim_vol}$ pc, $G < {maglim_vol}$")
    n2,_,_ = ax2.hist(period[eclipsing], bins_ecl, density=False, histtype='step', stacked=False, color=['k'], lw=3, label=f"Eclipsing, $G < {maglim_vol}$")
    
    # Error bars
    ax1.errorbar(bins[:-1]*10**(log_bin_width/2), n1, np.sqrt(n1), fmt='none', color='r')
    ax2.errorbar(bins_ecl[:-1]*10**(log_bin_width_ecl/2), n2, np.sqrt(n2), fmt='none', color='k')



    print(np.sum(eclipsing), 'eclipsing systems')
    print(np.sum(close), 'within 300 pc')
    print(np.sum(eclipsing&close), 'overlap')

    # plt.show()




    ### Run MCMC volume
    flat_samples, alpha, alpha_err = mcmc_fit(ln_prob_volume, period[close], sname='.chains_volume.npy')
    alpha_volume, alpha_volume_err = alpha[0], alpha_err[0]
    
    
    xplot = np.linspace(5,70,100)
    iters = np.array([period_prob_volume(xplot, sample) for sample in flat_samples])
    lims = np.percentile(iters,[2.3,15.9,50,84.1,97.7], axis=0) * \
                xplot * (10**(log_bin_width/2) - 10**-(log_bin_width/2)) * np.sum(close)
    ax1.fill_between(xplot, lims[1], lims[-2], color='red', alpha=0.2, label=f'Best fit range ($\\alpha={alpha_volume:.1f} \\pm {alpha_volume_err:.1f}$)')

    
    
    ### Run MCMC eclipsers
    
    flat_samples, alpha, alpha_err = mcmc_fit(ln_prob_eclipse, period[eclipsing], sname='.chains_eclipse.npy')
    alpha_eclipse, alpha_eclipse_err = alpha[0], alpha_err[0]


    xplot = np.linspace(pmin_ecl,70,100)
    iters_eclipse = np.array([period_prob_eclipse(xplot, sample) for sample in flat_samples])
    lims_eclipse = np.percentile(iters_eclipse,[2.3,15.9,50,84.1,97.7], axis=0) * \
                xplot * (10**(log_bin_width_ecl/2) - 10**-(log_bin_width_ecl/2)) * np.sum(eclipsing)
    ax2.fill_between(xplot, lims_eclipse[1], lims_eclipse[-2], color='k', alpha=0.2, label=f'Best fit range ($\\alpha={alpha_eclipse:.1f} \\pm {alpha_eclipse_err:.1f}$)')


    ### Deloye model

    xplot = np.linspace(5,70,100)
    ax1.plot(xplot, period_prob_volume(xplot, 3.6)*xplot * (10**(log_bin_width/2) - 10**-(log_bin_width/2)) * np.sum(close), \
             'k--', label='Deloye+ (2007) prediction ($\\alpha = 3.6$)')
    ax2.plot(xplot, period_prob_eclipse(xplot, 3.6)*xplot * (10**(log_bin_width_ecl/2) - 10**-(log_bin_width_ecl/2)) * np.sum(eclipsing), \
             'k--', label='Deloye+ (2007) prediction ($\\alpha = 3.6$)')

    ### Formatting
    for ax in [ax1, ax2]:
        xticks = [15, 20, 30, 40, 50, 60, 70]
        ax.set_xlim(15,70)
        ax.set_xscale('log')
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)

        yticks = [0,5,10,15,20,25]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        ax.legend(loc='upper left')

        ax.set_ylim(0,11)

    plt.savefig("fit_period_distribution.pdf")



    ### Est fraction of high state systems
    f_high_mean = quad(period_prob_plain, 10, 20, args=alpha_volume)[0]
    f_high_up = quad(period_prob_plain, 10, 20, args=alpha_volume-alpha_volume_err)[0]
    f_high_down = quad(period_prob_plain, 10, 20, args=alpha_volume+alpha_volume_err)[0]

    print('Fraction of high state systems, volume fit:', f_high_down, f_high_mean, f_high_up)
    print(f"log(high) - log(all): {np.log10(f_high_mean):.2f} + {np.log10(f_high_up)-np.log10(f_high_mean):.2f} - {np.log10(f_high_mean)-np.log10(f_high_down):.2f}")



    f_high_mean = quad(period_prob_plain, 10, 20, args=alpha_eclipse)[0]
    f_high_up = quad(period_prob_plain, 10, 20, args=alpha_eclipse-alpha_eclipse_err)[0]
    f_high_down = quad(period_prob_plain, 10, 20, args=alpha_eclipse+alpha_eclipse_err)[0]

    print('Fraction of high state systems, eclipse fit:', f_high_down, f_high_mean, f_high_up)
    print(f"log(high) - log(all): {np.log10(f_high_mean):.2f} + {np.log10(f_high_up)-np.log10(f_high_mean):.2f} - {np.log10(f_high_mean)-np.log10(f_high_down):.2f}")




    plt.show()