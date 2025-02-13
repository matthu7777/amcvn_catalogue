#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
import mgutils as mg, mgutils.constants as co


def write_table(fname, oi):

    with open(fname, 'w') as f:
        f.write('')

        for row in data[oi]:
            rowstring = " & ".join([
                row['Name'],
                row['Coords'],
                f"{row['Period']:.2f}"+row['Period_note'] if row['Period']<10 else f"{row['Period']:.1f}"+row['Period_note'] if np.isfinite(row['Period']) else '--',
                row['Disk_state'] if row['Disk_state'] else '??',
                f"${row['Gaia_Gmag']:.2f} \\pm {row['Gaia_Gmag_err']:.2f}$" if np.isfinite(row['Gaia_Gmag']) else '--',
                f"${row['Q']:.2f} \\pm {row['Q_err']:.2f}$" if np.isfinite(row['Q']) else '--',
                "\\checkmark" if row['Has_spectrum'] else "$\\times$"
            ]) + " \\\\\n"
            f.write(rowstring)






if __name__ == '__main__':

    if "-h" in argv:
        print ("generate_table_latex.py usage")
        raise SystemExit


    with fits.open("amcvn_catalogue.fits") as hdul:
        data = hdul[1].data
    

    write_table("latex_table_amcvns.tex", data['AM_CVn'])

    write_table("latex_table_hecvs.tex", data['He_CV'])

    write_table("latex_table_others.tex", data['Confirmed']&(~data['AM_CVn'] & ~data['He_CV']))

    write_table("latex_table_candidates.tex", ~data['Confirmed'])
    


