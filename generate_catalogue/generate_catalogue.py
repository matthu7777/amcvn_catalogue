#!/usr/bin/env python

"""
"""
    
import numpy as np
# import matplotlib.pyplot as plt
from sys import argv, exit
# from astropy.io import fits
from astropy.table import Table
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import re
import os
import time
# import copy
# from itertools import zip_longest


def convertRA(hr,min,sec):
    hr = float(hr)
    min = float(min)
    sec = float(sec)
    return (hr + (min + sec/60.)/60.) * 360. / 24.
    
def convertDec(deg, min, sec):
    deg = float(deg)
    min = float(min)
    sec = float(sec)
    if np.size(deg) == 1:
        sign = +1 if deg==0 else np.sign(deg)
    else:
        sign = np.sign(deg)
        sign[deg==0] = +1
    return sign * (np.abs(deg) + (min + sec/60.)/60.)

def convertRADec(radecstr):
    radec = re.split('[:, ]', radecstr)
    ra = convertRA(radec[0],radec[1],radec[2])
    dec = convertDec(radec[3],radec[4],radec[5])
    return ra, dec

def get_sheet_column(sheets, cols, col_name, length=None):
    vals = np.array(sheets.col_values(np.argwhere(cols==col_name)[0,0]+1)[1:], dtype=str)
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


def findQMcAllisterB(eps,epserr):
    """ McAllister 2018's update to Knigge's eps-q relation (with errors) for phase B superhumps
    """
    def calculate(a,b,c,eps):
        return a + b*(eps-c)
    return functApp(calculate,0.118,0.003,4.45,0.28,0.025,0,eps,epserr,)
    

def functApp(funct, *args, verbose=False, plus=True, minus=True):
    """ Functional approach on an arbitrary function. Usage: functApp(funct, a, aerr, b, berr ...). 
        Use 'plus' and 'minus' to toggle whether you want to propogate by adding to the values or subtracting, or both and average.
    """
    var, err = [], []
    if len(args) % 2 != 0:
        # If odd number of args passed -- ie if they've missed off an error or something
        print ("Functional approach has been passed the wrong number of arguments (%d excluding the function pointer)".format(len(args)))
        exit()
    #split variables and errors into 2 arrays
    for i, arg in enumerate(args):
        if i % 2 == 0:
            var.append(float(arg))
        else:
            err.append(float(arg))
    var = np.array(var)
    err = np.array(err)
    
    #Find the 'expected' result
    result = funct(*var)
    
    # For each error, propogate it's effect through
    # Let user choose whether to add or take away error, or do both and average
    if plus and minus:
        diffs = []
        for j, (v, e) in enumerate(zip(var, err)):
            toPass1 = np.copy(var)
            toPass2 = np.copy(var)
            toPass1[j] = v + e   
            toPass2[j] = v - e   
            d1 = funct(*toPass1) - result
            d2 = funct(*toPass2) - result
            
            # Avoid any nans that might have come up if only on one side
            if isinstance(d1,float) or isinstance(d1,int):
                if np.isnan(d1):
                    d1 = d2
                if np.isnan(d2):
                    d2 = d1
            else:
                d1[np.isnan(d1)] = d2[np.isnan(d1)]
                d2[np.isnan(d2)] = d1[np.isnan(d2)]
            
            
            diffs.append((np.abs(d1)+np.abs(d2))/2.)
        diffs = np.array(diffs)
    elif plus:
        diffs = []
        for j, (v, e) in enumerate(zip(var, err)):
            toPass = np.copy(var)
            toPass[j] = v + e   
            d1 = funct(*toPass) - result
            diffs.append(d1)
        diffs = np.array(diffs)
    elif minus:
        diffs = []
        for j, (v, e) in enumerate(zip(var, err)):
            toPass = np.copy(var)
            toPass[j] = v - e   
            d1 = funct(*toPass) - result
            diffs.append(d1)
        diffs = np.array(diffs)
    
    if verbose:
        print ("Error contributions:", diffs)
    # combine error contributions in quadrature
    return result, quad(diffs)
    
def quad(*args):
    """ Combine a list of arguments in quadrature
    """
    args = np.array(args)
    return np.sqrt(np.sum(args**2))


if __name__ == '__main__':

    if "-h" in argv:
        print ("generate_catalogue.py usage")
        raise SystemExit


    # Add a delay between API requests to avoid the Google requests-per-minute limit
    use_cautious_timing = True


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
        j for j,p in enumerate(
            sheet.col_values(np.argwhere(cols=='Include in published table')[0,0]+1)[1:]
        ) if p=='y'
    ])

    # Initialise table
    table = Table()
    table['Name'] = names[oi]

    # Booleans
    # table['Include'] = np.array(get_sheet_column(sheet, cols, 'Include in published table', length=len(names))[oi] == 'y', dtype=bool)
    table['Confirmed'] = np.array(get_sheet_column(sheet, cols, 'Confirmed?', length=len(names))[oi] == 'y', dtype=bool)
    table['Has_spectrum'] = np.array(get_sheet_column(sheet, cols, 'Has spectrum', length=len(names))[oi] == 'y', dtype=bool)
    table['Has_optical_hydrogen'] = np.array(get_sheet_column(sheet, cols, 'Has optical hydrogen lines', length=len(names))[oi] == 'y', dtype=bool)
    table['Has_optical_donor'] = np.array(get_sheet_column(sheet, cols, 'Has visible donor', length=len(names))[oi] == 'y', dtype=bool)

    # Other names
    table['Other_names'] = get_sheet_column(sheet, cols, 'Other names', length=len(names))[oi]

    # Coords
    coords_hex = get_sheet_column(sheet, cols, 'RA/Dec')[oi]
    table['RA'] = np.array([convertRADec(radecstr)[0] for radecstr in coords_hex])
    table['Dec'] = np.array([convertRADec(radecstr)[1] for radecstr in coords_hex])
    for c in coords_hex:
        # assert len(c) == 23
        assert len(c.split())==2
        assert len(c.split(':'))==5
        assert c.split()[1].startswith('+') or c.split()[1].startswith('-')
    table['Coords_hex'] = coords_hex

    # Discovery paper
    table['Discovery_ref'] = [treat_reference(r) for r in get_sheet_column(sheet, cols, 'Discovery paper')[oi]]
    
    # Comment
    table['Notes'] = get_sheet_column(sheet, cols, 'Notes', length=len(names))[oi]

    # Period
    period_strs = get_sheet_column(sheet, cols, 'Period (min)', length=len(names))[oi]
    table['Period'] = [float(re.findall(r'\d+[.]?\d*', p)[0]) if len(p)>0 else np.nan for p in period_strs]
    table['Period_note'] = [p.replace(re.findall(r'\d+[.]?\d*', p)[0],'').strip() if len(p)>0 else '' for p in period_strs]
    table['Period_method'] = get_sheet_column(sheet, cols, 'Period Method', length=len(names))[oi]
    table['Period_comment'] = get_sheet_column(sheet, cols, 'Period Comment', length=len(names))[oi]
    table['Period_ref'] = [treat_reference(r) for r in get_sheet_column(sheet, cols, 'Period Source', length=len(names))[oi]]

    if use_cautious_timing:
        time.sleep(3)

    # Disc state
    table['Disk_state'] = get_sheet_column(sheet, cols, 'State', length=len(names))[oi]

    # Mass ratio
    table['Q'] = [float(v) if len(v)>0 else np.nan for v in get_sheet_column(sheet, cols, 'q', length=len(names))[oi]]
    table['Q_err'] = [float(v) if len(v)>0 else np.nan for v in get_sheet_column(sheet, cols, 'q err', length=len(names))[oi]]
    table['Q_method'] = get_sheet_column(sheet, cols, 'q Method', length=len(names))[oi]
    table['Q_ref'] = get_sheet_column(sheet, cols, 'q Source', length=len(names))[oi]
    table['Q_comment'] = get_sheet_column(sheet, cols, 'q Comment', length=len(names))[oi]

    if use_cautious_timing:
        time.sleep(3)

    # Superhump
    table['Superhump_excess'] = [float(v) if len(v)>0 else np.nan for v in get_sheet_column(sheet, cols, 'Superhump Excess', length=len(names))[oi]]
    table['Superhump_excess_err'] = [float(v) if len(v)>0 else np.nan for v in get_sheet_column(sheet, cols, 'Superhump err', length=len(names))[oi]]
    table['Superhump_ref'] = get_sheet_column(sheet, cols, 'Superhump Source', length=len(names))[oi]
    table['Superhump_comment'] = get_sheet_column(sheet, cols, 'Superhump Comment', length=len(names))[oi]

    if use_cautious_timing:
        time.sleep(3)

    # Masses
    table['M1'] = [float(v) if len(v)>0 else np.nan for v in get_sheet_column(sheet, cols, 'M1', length=len(names))[oi]]
    table['M1_err'] = [float(v) if len(v)>0 else np.nan for v in get_sheet_column(sheet, cols, 'M1 err', length=len(names))[oi]]
    table['M2'] = [float(v) if len(v)>0 else np.nan for v in get_sheet_column(sheet, cols, 'M2', length=len(names))[oi]]
    table['M2_err'] = [float(v) if len(v)>0 else np.nan for v in get_sheet_column(sheet, cols, 'M2 err', length=len(names))[oi]]
    table['Masses_ref'] = get_sheet_column(sheet, cols, 'M1 M2 Source', length=len(names))[oi]


    


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

    spec_refs = np.array([], dtype=str)
    refs_targets = np.array([], dtype=str)
    all_all_refs = np.array([], dtype=str)
    for j, (pc, sc) in enumerate(zip(phot_comment, spec_comment)):
        # read spectroscopy refs
        spec_refs_target = np.unique(re.findall(r'[A-Za-z]+\d\d\d\d[a-z]?', sc))
        spec_refs_target = np.append(spec_refs_target, [n for n in nonstandard_refs if n in pc+' '+sc])

        # read all refs from comment columns
        extra_refs = re.findall(r'[A-Za-z]+\d\d\d\d[a-z]?', pc+' '+sc)
        extra_refs = np.append(extra_refs, [n for n in nonstandard_refs if n in pc+' '+sc])
        
        # compile all references for this particular target
        all_refs = np.append(extra_refs, treat_reference(table['Discovery_ref'][j]))
        if len(table['Period_ref'][j]) > 0:
            all_refs = np.append(all_refs, [treat_reference(s) for s in table['Period_ref'][j].split(',')])
        if len(table['Q_ref'][j]) > 0:
            all_refs = np.append(all_refs, [treat_reference(s) for s in table['Q_ref'][j].split(',')])
        if len(table['Superhump_ref'][j]) > 0:
            all_refs = np.append(all_refs, [treat_reference(s) for s in table['Superhump_ref'][j].split(',')])
        if len(table['Masses_ref'][j]) > 0:
            all_refs = np.append(all_refs, [treat_reference(s) for s in table['Masses_ref'][j].split(',')])

        # tidy format
        all_refs = np.unique(all_refs)
        ref_years = np.array([int(re.findall(r'\d\d\d\d', r)[0]) if len(re.findall(r'\d\d\d\d', r))>0 else \
                        np.nan if (r=='' or r=='--' or 'vsnet' in r) else 
                        nonstandard_refs[r] for r in all_refs])
        all_refs = all_refs[np.argsort(ref_years[~np.isnan(ref_years)])]
        spec_refs_target = np.unique(spec_refs_target)
        ref_years_spec = np.array([int(re.findall(r'\d\d\d\d', r)[0]) if len(re.findall(r'\d\d\d\d', r))>0 else \
                        np.nan if (r=='' or r=='--' or 'vsnet' in r) else 
                        nonstandard_refs[r] for r in spec_refs_target])
        spec_refs_target = spec_refs_target[np.argsort(ref_years_spec[~np.isnan(ref_years_spec)])]

        # append to final columns
        if len(spec_refs_target) > 0:
            spec_refs = np.append(spec_refs, ', '.join(spec_refs_target[(spec_refs_target!='--')]))
        else:
            spec_refs = np.append(spec_refs, [''])
        refs_targets = np.append(refs_targets, ', '.join(all_refs[(all_refs!='--')]))
        all_all_refs = np.append(all_all_refs, all_refs)

    table['Spec_refs'] = spec_refs
    table['References'] = refs_targets


    ### Find discovery year
    ref_years = np.array([
                        int(y) if y else \
                        int(re.findall(r'\d\d\d\d', r)[0]) if len(re.findall(r'\d\d\d\d', r))>0 else \
                        np.nan for r,y in zip(table['Discovery_ref'], \
                        get_sheet_column(sheet, cols, 'Discovery year (unpublished systems only)', length=len(names))[oi])
    ])
    table['Discovery_year'] = ref_years



    ### Check if any values are missing references, errors, and a few other checks
    has_period = table['Period'] > 0
    has_q = table['Q'] > 0
    has_sh = table['Superhump_excess'] > 0
    has_m1 = table['M1'] > 0
    has_m2 = table['M2'] > 0
    # period
    for j, ref in enumerate(table['Period_ref'][has_period]):
        assert ref
    # spectrum 
    for j, ref in enumerate(table['Spec_refs'][table['Has_spectrum']]):
        assert ref
    # q
    for j, ref in enumerate(table['Q_ref'][has_q]):
        assert ref
        assert table['Q_err'][has_q][j] > 0
        # If there's a superhump-derived q, check that we have the original superhump excess as well
        if table['Q_method'][has_q][j] == 'Superhumps':
            assert has_sh[has_q][j]
    # superhump
    for j, ref in enumerate(table['Superhump_ref'][has_sh]):
        assert ref
        assert table['Superhump_excess_err'][has_sh][j] > 0
    # M1
    for j, ref in enumerate(table['Masses_ref'][has_m1]):
        assert ref
        assert table['M1_err'][has_m1][j] > 0
    # M2
    for j, ref in enumerate(table['Masses_ref'][has_m2]):
        assert ref
        assert table['M2_err'][has_m2][j] > 0
    # Any with both M1 and M2 should have Q as well
    for j, ref in enumerate(table['Masses_ref'][has_m1&has_m2]):
        assert table['Q'][has_m1&has_m2][j] > 0
    

    ### Recalculate superhump mass ratios
    for j, (sh, sh_err) in enumerate(zip(table['Superhump_excess'], table['Superhump_excess_err'])):
        if (not np.isfinite(table['Superhump_excess'][j])) or (table['Q_method'][j] and table['Q_method'][j] != 'Superhumps'):
            continue
        table['Q'][j], table['Q_err'][j] = findQMcAllisterB(sh, sh_err)
        table['Q_ref'][j] = 'Recalculated in this work'
        pass

        



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





    ######### Gaia crossmatch to go here


    # Old version:
            

    # ### Create temporary table for Gaia crossmatch

    # temp_table = copy.copy(table)
    # mag_strings = get_sheet_column(sheet, cols, 'Magnitude, band', length=len(names))[oi]
    # mags = np.array([])
    # for m in mag_strings:
    #     if m:
    #         mag = re.findall(r'\d+[.]?\d*', m)
    #         assert len(mag) == 1 or np.isnan(mag)
    #         mags = np.append(mags, float(mag[0]))
    #     else:
    #         mags = np.append(mags, np.nan)
    # temp_table['Magnitude'] = mags








    # temp_table.write('for_gaia_upload.fits', overwrite=True)


    # ### Now go away and do the Gaia and Bailer-Jones crossmatches in Topcat

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

        table['Distance'] = match_arrays(bjdist['rgeo'], inds_table_bjdist, len(table), float)
        table['Distance_16percentile'] = match_arrays(bjdist['b_rgeo_x'], inds_table_bjdist, len(table), float)
        table['Distance_84percentile'] = match_arrays(bjdist['B_rgeo_xa'], inds_table_bjdist, len(table), float)
    









    ### Few final formatting things, remove unwanted columns

    table['not confirmed'] = ~table['Confirmed']        # a quick fudge to make sorting easier
    table.sort(['Period', 'not confirmed', 'Name'])
    table.remove_column('not confirmed')

    # #### To check for non-ASCII characters that cause the file-write to crash (bug-fixing only, you should normally comment these lines)
    # for row in table:
    #     Table(row).write('temp.fits', overwrite=True)
    # ####

    table.write('amcvn_catalogue.fits', overwrite=True)



    print(f"The catalogue contains {np.sum(table['Confirmed'])} AM CVns and {np.sum(~table['Confirmed'])} candidates")
    print(f"{np.sum(np.isfinite(table['Period'][table['Confirmed']]))} AM CVns and " + \
                    f"{np.sum(np.isfinite(table['Period'][~table['Confirmed']]))} candidates have periods")
    # print(f"{np.sum(table['Gaia_ID'][table['Confirmed']]!='')} AM CVns and " + \
    #                 f"{np.sum(table['Gaia_ID'][~table['Confirmed']]!='')} candidates have a Gaia crossmatch")
    # print(f"{np.sum(np.isfinite(table['Distance_mean'][table['Confirmed']]))} AM CVns and " + \
    #                 f"{np.sum(np.isfinite(table['Distance_mean'][~table['Confirmed']]))} candidates have a distance measurement")



    print(table)
    print('Done') 