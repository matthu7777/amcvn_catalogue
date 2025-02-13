#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import mgutils as mg, mgutils.constants as co



rng = np.random.default_rng(42)

n = 100_000

x_bright = rng.uniform(-1000,1000,n//10)
y_bright = rng.uniform(-1000,1000,n//10)
x_faint = rng.uniform(-1000,1000,n)
y_faint = rng.uniform(-1000,1000,n)

dist_bright = np.sqrt(x_bright**2 + y_bright**2)
dist_faint = np.sqrt(x_faint**2 + y_faint**2)


absmag_bright = rng.normal(8, 0.1, n//10)
absmag_faint = rng.normal(9, 0.1, n)

appmag_bright = absmag_bright + 5 * np.log10(dist_bright/10)
appmag_faint = absmag_faint + 5 * np.log10(dist_faint/10)

scale = np.sum(appmag_faint < 16) / np.sum(appmag_bright < 16)


plt.hist(appmag_bright, 20, cumulative=True, histtype='step', weights=np.ones(n//10)*scale)
plt.hist(appmag_faint, 20, cumulative=True, histtype='step')

plt.yscale('log')
plt.show()