#!/usr/bin/env python

"""
"""
    
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from astropy.io import fits
from astropy.table import Table
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import copy
import re
import os
from itertools import zip_longest
import mgutils as mg, mgutils.constants as co


def get_sheet_column(sheet, cols, col_name, length=None):
    vals = np.array(sheet.col_values(np.argwhere(cols==col_name)[0,0]+1)[1:], dtype=str)
    if length and len(vals) < length:
        vals = np.append(vals, ['']*(length-len(vals)))
        assert len(vals) == length
    return vals


def treat_reference(ref):
    return ref.replace('(ATel)','').strip()


def match_arrays(array, inds, length, dtype=str):
    fill_value = '' if dtype is str else np.nan
    new_array = [array[inds.index(j)] if j in inds else fill_value for j in range(length)]
    return np.array(new_array, dtype=dtype)




if __name__ == '__main__':

    if "-h" in argv:
        print ("generate_catalogue.py usage")
        raise SystemExit


    ### Read from Google sheet

    # Authorize the API
    scope = [
        'https://www.googleapis.com/auth/drive',
        'https://www.googleapis.com/auth/drive.file'
        ]
    file_name = 'client_key.json'
    creds = ServiceAccountCredentials.from_json_keyfile_name(file_name,scope)
    client = gspread.authorize(creds)


    # Fetch the sheet
    sheet = client.open('AM CVn population').sheet1

    # Fetch data from sheet
    names = np.array(sheet.col_values(1)[1:], dtype=str)
    cols = np.array(sheet.row_values(1), dtype=str)
    oi = np.array([
        j for j,(p,c) in enumerate(zip_longest(sheet.col_values(np.argwhere(cols=='Include in published table')[0,0]+1)[1:], \
                                       sheet.col_values(np.argwhere(cols=='Include in candidate table')[0,0]+1)[1:])) if (p=='y' or c=='y')
    ])

    # Initialise table
    table = Table()
    table['Name'] = names[oi]

    # Confirmed or not
    table['Confirmed'] = np.array(get_sheet_column(sheet, cols, 'Include in published table', length=len(names))[oi] == 'y', dtype=bool)

    # Other names
    table['Other_names'] = get_sheet_column(sheet, cols, 'Other names', length=len(names))[oi]

    # Coords
    coords_hex = get_sheet_column(sheet, cols, 'RA/Dec')[oi]
    table['RA'] = np.array([mg.convertRADec(radecstr)[0] for radecstr in coords_hex])
    table['Dec'] = np.array([mg.convertRADec(radecstr)[1] for radecstr in coords_hex])
    for c in coords_hex:
        assert len(c) == 23
        assert len(c.split())==2
        assert len(c.split(':'))==5
        assert '+' in c or '-' in c
    table['Coords_hex'] = coords_hex

    # Discovery paper
    table['Discovery_ref'] = [treat_reference(r) for r in get_sheet_column(sheet, cols, 'Discovery paper')[oi]]
    for d,c in zip(table['Discovery_ref'],table['Confirmed']):
        assert d and (d != '--' or ~c)
    
    # Comment
    table['Notes'] = get_sheet_column(sheet, cols, 'Notes', length=len(names))[oi]

    # Period
    period_strs = get_sheet_column(sheet, cols, 'Period (min)', length=len(names))[oi]
    table['Period'] = [float(re.findall(r'\d+[.]?\d*', p)[0]) if len(p)>0 else np.nan for p in period_strs]
    table['Period_note'] = [p.replace(re.findall(r'\d+[.]?\d*', p)[0],'').strip() if len(p)>0 else '' for p in period_strs]
    table['Period_method'] = get_sheet_column(sheet, cols, 'Period Method', length=len(names))[oi]
    table['Period_comment'] = get_sheet_column(sheet, cols, 'Period Comment', length=len(names))[oi]
    table['Period_ref'] = [treat_reference(r) for r in get_sheet_column(sheet, cols, 'Period Source', length=len(names))[oi]]

    # Disc state
    table['Disk_state'] = get_sheet_column(sheet, cols, 'State', length=len(names))[oi]
    

    # All references
    # In the dict below, put any references that are not of the form aaaaaa0000, and their publication year
    # If you put nan, it will be excluded from the table (use this for any priv comm, for instance)
    nonstandard_refs = {
        'RoelofsPhD':2007,
        'Kupfer (priv. comm.)':np.nan,
        'Aungwerojwit (priv. comm.)':np.nan,
        'Motsoaledi et al. (priv. comm.)':np.nan,
        'Green et al (in review)':np.nan,
        # '--':np.nan,
    } 
    phot_comment = get_sheet_column(sheet, cols, 'Photometry', length=len(names))[oi]
    spec_comment = get_sheet_column(sheet, cols, 'Spectroscopy', length=len(names))[oi]
    refs_targets = np.array([], dtype=str)
    all_all_refs = np.array([], dtype=str)
    for j, (pc, sc) in enumerate(zip(phot_comment, spec_comment)):
        extra_refs = re.findall(r'[A-Za-z]+\d\d\d\d[a-z]?', pc+' '+sc)
        extra_refs = np.append(extra_refs, [n for n in nonstandard_refs if n in pc+' '+sc])
        all_refs = np.append(extra_refs, treat_reference(table['Discovery_ref'][j]))
        if len(table['Period_ref'][j]) > 0:
            all_refs = np.append(all_refs, [treat_reference(s) for s in table['Period_ref'][j].split(',')])
        all_refs = np.unique(all_refs)
        ref_years = np.array([int(re.findall(r'\d\d\d\d', r)[0]) if len(re.findall(r'\d\d\d\d', r))>0 else \
                        np.nan if (r=='' or r=='--' or 'vsnet' in r) else 
                        nonstandard_refs[r] for r in all_refs])
        all_refs = all_refs[np.argsort(ref_years[~np.isnan(ref_years)])]
        refs_targets = np.append(refs_targets, ', '.join(all_refs[(all_refs!='--')]))
        all_all_refs = np.append(all_all_refs, all_refs)
    table['References'] = refs_targets

    # Find discovery year
    ref_years = np.array([
                        int(y) if y else \
                        int(re.findall(r'\d\d\d\d', r)[0]) if len(re.findall(r'\d\d\d\d', r))>0 else \
                        np.nan for r,y in zip(table['Discovery_ref'], \
                        get_sheet_column(sheet, cols, 'Discovery year (candidates only)', length=len(names))[oi])
    ])
    table['Discovery_year'] = ref_years


    # Check if any periods are missing references
    has_period = table['Period'] > 0
    for j, ref in enumerate(table['Period_ref'][has_period]):
        assert ref


    ### Output references 

    all_all_refs = np.unique(all_all_refs)
    # ref_years = np.array([int(re.findall(r'\d\d\d\d', r)[0]) if len(re.findall(r'\d\d\d\d', r))>0 else \
    #                     np.nan if r=='' else nonstandard_refs[r] for r in all_all_refs])
    # all_all_refs = all_all_refs[np.argsort(ref_years[~np.isnan(ref_years)])]

    with open('references.txt', 'w') as f:
        f.write("")
    with open('references.txt', 'a') as f:
        for ref in all_all_refs:
            f.write(ref+" & \\citet{"+ref+"}\\\\\n")



    ### Create temporary table for Gaia crossmatch

    temp_table = copy.copy(table)
    mag_strings = get_sheet_column(sheet, cols, 'Magnitude, band', length=len(names))[oi]
    mags = np.array([])
    for m in mag_strings:
        if m:
            mag = re.findall(r'\d+[.]?\d*', m)
            assert len(mag) == 1 or np.isnan(mag)
            mags = np.append(mags, float(mag[0]))
        else:
            mags = np.append(mags, np.nan)
    temp_table['Magnitude'] = mags

    temp_table.write('for_gaia_upload.fits', overwrite=True)


    ### Now go away and do the Gaia and Bailer-Jones crossmatches in Topcat

    ### Read in Gaia crossmatch and get data
    # Do the Gaia crossmatch with at least 12 arcseconds to include GP Com and V396 Hya
    # Then we will filter out all some bad crossmatches

    if os.path.isfile('gaia_crossmatch.fits'):
        gaia = Table.read('gaia_crossmatch.fits')
        gaia_ok = (gaia['angDist'] < 5) | (gaia['PM'] > 100)
        gaia = gaia[gaia_ok]

        inds_table_gaia = [np.argwhere(name.strip()==table['Name'])[0,0] for name in gaia['Name']]

        table['Gaia_ID'] = match_arrays(gaia['Source'], inds_table_gaia, len(table), str)
        table['Gaia_Gmag'] = match_arrays(gaia['Gmag'], inds_table_gaia, len(table), float)
        table['Gaia_Gmag_err'] = match_arrays(gaia['e_Gmag'], inds_table_gaia, len(table), float)
        table['Gaia_BPmag'] = match_arrays(gaia['BPmag'], inds_table_gaia, len(table), float)
        table['Gaia_BPmag_err'] = match_arrays(gaia['e_BPmag'], inds_table_gaia, len(table), float)
        table['Gaia_RPmag'] = match_arrays(gaia['RPmag'], inds_table_gaia, len(table), float)
        table['Gaia_RPmag_err'] = match_arrays(gaia['e_RPmag'], inds_table_gaia, len(table), float)
        table['Gaia_parallax'] = match_arrays(gaia['Plx'], inds_table_gaia, len(table), float)
        table['Gaia_parallax_err'] = match_arrays(gaia['e_Plx'], inds_table_gaia, len(table), float)


    ### Read in Bailer-Jones crossmatch and get data

    if os.path.isfile('bjdist_crossmatch.fits'):
        bjdist = Table.read('bjdist_crossmatch.fits')
        # bjdist_ok = (bjdist['angDist'] < 5) | (bjdist['PM'] > 100)
        bjdist_ok = np.array([n in gaia['Name'] for n in bjdist['Name']], dtype=bool)
        bjdist = bjdist[bjdist_ok]

        inds_table_bjdist = [np.argwhere(name.strip()==table['Name'])[0,0] for name in bjdist['Name']]

        table['Distance_mean'] = match_arrays(bjdist['rgeo'], inds_table_bjdist, len(table), float)
        table['Distance_16percentile'] = match_arrays(bjdist['b_rgeo_x'], inds_table_bjdist, len(table), float)
        table['Distance_84percentile'] = match_arrays(bjdist['B_rgeo_xa'], inds_table_bjdist, len(table), float)
    

    ### Few final formatting things, remove unwanted columns

    table.sort('Period')

    table.write('amcvn_catalogue_secret.fits', overwrite=True)

    table.remove_column('Discovery_year')

    confirmed_table = table[table['Confirmed']]
    confirmed_table.remove_column('Confirmed')
    confirmed_table.write('amcvn_catalogue_confirmed.fits', overwrite=True)

    candidate_table = table[~table['Confirmed']]
    candidate_table.remove_column('Confirmed')
    candidate_table.write('amcvn_catalogue_candidate.fits', overwrite=True)



    print(f"The catalogue contains {np.sum(table['Confirmed'])} AM CVns and {np.sum(~table['Confirmed'])} candidates")
    print(f"{np.sum(np.isfinite(table['Period'][table['Confirmed']]))} AM CVns and " + \
                    f"{np.sum(np.isfinite(table['Period'][~table['Confirmed']]))} candidates have periods")
    print(f"{np.sum(table['Gaia_ID'][table['Confirmed']]!='')} AM CVns and " + \
                    f"{np.sum(table['Gaia_ID'][~table['Confirmed']]!='')} candidates have a Gaia crossmatch")
    print(f"{np.sum(np.isfinite(table['Distance_mean'][table['Confirmed']]))} AM CVns and " + \
                    f"{np.sum(np.isfinite(table['Distance_mean'][~table['Confirmed']]))} candidates have a distance measurement")



    print(table)
    print('Done') 