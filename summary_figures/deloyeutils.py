#!/usr/bin/env python


""" A messy collection of functions to use for plotting Deloye mass-radius diagrams
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sys import argv
import string
from deloyelines import *
from donor_mr_plot import *


#### A load of constants

#obvious ones
c = 299792458.
G = 6.67408e-11
h = 6.626e-34   #planck
planck = h
H0 = 2.25e-18   # Hubble (in s^-1)
hubble = H0
H0Mpc = 67.6 # Hubble constant in km/s/Mpc
#subatomic
e = 1.602e-19
eV = e
me = 9.109e-31
mp = 1.673e-27
#thermodynamics
kB = 1.3806e-23     #boltzmann
stefanBoltzmann = 5.67e-8
#standard conditions
g = 9.80665
atm = 101325


#time
secondsPerDay = 86400.
secondsPerYear = 86400. * 365.25
minutesPerDay = 86400./60
tH = 4.55e17    #Hubble time in seconds


#distance
ly = 9.461e15
pc = 3.086e16
au = 1.496e11

#energy
mJy = 1e-29

##properties of objects
#mass
mEarth = 5.972e24
mJup = 1.898e27
mNep = 1.024e26
mSun = 1.98855e30
mSol = mSun
msol = mSun
msun = mSun
#radius
rEarth = 6.371e6
rSun = 695.7e6
rSol = rSun
rsol = rSun
rsun = rSun
#others
lSun = 3.828e26
lSol = lSun




### Functions for reading from xfig and plotting Deloye lines

# 7944 corresponds to -1.5 on y axis
# 6256 to -1.3
# 8788 to -1.6
# 9164 to -0.5 on x axis
# 2564 to -2.4 

xf1 = 9164
xf2 = 2564
xr1 = -0.5
xr2 = -2.4
yf1 = 6256
yf2 = 8788
yr1 = -1.3
yr2 = -1.6

    
comments = {\
"HP Lib":"Additional constraints on the donor masses in these systems were derived based on the magnitude and parallax measurements",
"CR Boo":"Additional constraints on the donor masses in these systems were derived based on the magnitude and parallax measurements",
"V803 Cen":"Additional constraints on the donor masses in these systems were derived based on the magnitude and parallax measurements",
"YZ LMi":"The tabulated $q$ for YZ\\,LMi is derived from eclipse photometry; its superhump excess gives $q = 0.049 \pm 0.008$ by Equation~\\ref{eq:mrknigge}",
"AM CVn":"The tabulated $q$ for AM\\,CVn is derived from spectroscopy; its superhump excess gives $q = 101 \pm 0.005$ by Equation~\\ref{eq:mrknigge}",
}

def formatGraph(fignum=1,suptitle=None,xlabel=None,ylabel=None,returnFig=False,poster=False,*args,**kwargs):
    """ General formatting things like style sheet, axes etc. Title and axis labels are applied to the current figure and axes.
    """
    matplotlib.rc('font',**{'family':'serif', 'size':16})
    matplotlib.rc('text', usetex=True)
    
    fig = plt.figure(fignum, *args,**kwargs)
    if xlabel is not None:
        plt.gca().set_xlabel(xlabel)
    if ylabel is not None:
        plt.gca().set_ylabel(ylabel)
    if suptitle is not None:
        plt.suptitle(suptitle)
    
    
    
    if returnFig:
        return fig
    else:
        return plt.gca()



def fitLine(x1,y1,x2,y2):
    """ Fits two points with an y = mx + c line.
    """
    m = (y2 - y1) / (x2 - x1)
    c = y1 - m * x1
    return m, c

def convert(coords, figCoord1, realCoord1, figCoord2, realCoord2):
    """ Convert from xfig coordinates to real coordinates.
    """
    m, c = fitLine(figCoord1, realCoord1, figCoord2, realCoord2)
    
    return coords * m + c
    

def makeCoordinates(string):
    """ Takes a string of form "x y x y" and returns 2d numpy array [[x,y],[x,y]]
    """
    string = string.replace('\n','')
    return np.reshape(string.split(),(len(string.split())//2,2)).astype(np.float)


def plot(string,fmt='b-', **kwargs):
    """ Wrapper around all the things that need to happen to plot each deloye line.
    """
    a = makeCoordinates(string)
    x = convert(a[:,0], xf1, xr1, xf2, xr2)
    y = convert(a[:,1], yf1, yr1, yf2, yr2)
    plt.plot(x, y, fmt, **kwargs)
    return x,y


def plotDeloyeModels(fmt=None, label=None, twoDeloye=False, **kwargs):
    if fmt is None:
        x,y = plot(line1, 'y-', **kwargs)     #yellow
        x,y = plot(line2, 'y-', **kwargs)     #yellow
        x,y = plot(line3, 'r-', **kwargs)     #red
        x,y = plot(line4, 'r-', **kwargs)     #red
        #plot(line5,'rD')   #??
        x,y = plot(line6, 'g-', **kwargs)     #green
        x,y = plot(line7, 'g-', **kwargs)     #green
        x,y = plot(line8, 'b-', **kwargs)     #blue
        x,y = plot(line9, 'b-', **kwargs)     #blue
        #plot(line10,'rD')  #??
        x,y = plot(line11, 'y-', **kwargs)    #yellow
        x,y = plot(line12, 'y-', **kwargs)    #yellow
        x,y = plot(line13, 'r-', **kwargs)    #red
        x,y = plot(line14, 'r-', **kwargs)    #red
        #plot(line15,'rD')  #??
        x,y = plot(line16, 'g-', **kwargs)    #green
        x,y = plot(line17, 'g-', **kwargs)    #green
        x,y = plot(line18, 'b-', **kwargs)    #blue
        x,y = plot(line19, 'b-', **kwargs)    #blue
    else:
        x1,y1 = plot(line1, fmt, label=label, **kwargs)        #yellow
        x2,y2 = plot(line2, fmt, **kwargs)      #yellow
        np.savetxt("deloye_y1.dat", np.column_stack((np.append(x1,x2),np.append(y1,y2))))
        x1,y1 = plot(line3, fmt, **kwargs)      #red
        x2,y2 = plot(line4, fmt, **kwargs)      #red
        np.savetxt("deloye_r1.dat", np.column_stack((np.append(x1,x2),np.append(y1,y2))))
        #plot(line5,fmt)    #??
        x1,y1 = plot(line6, fmt, **kwargs)      #green
        x2,y2 = plot(line7, fmt, **kwargs)      #green
        np.savetxt("deloye_g1.dat", np.column_stack((np.append(x1,x2),np.append(y1,y2))))
        x1,y1 = plot(line8, fmt, **kwargs)      #blue
        x2,y2 = plot(line9, fmt, **kwargs)      #blue
        np.savetxt("deloye_b1.dat", np.column_stack((np.append(x1,x2),np.append(y1,y2))))
        if twoDeloye:
            #plot(line10,fmt)   #??
            x1,y1 = plot(line11, fmt, **kwargs)    #yellow
            x2,y2 = plot(line12, fmt, **kwargs)    #yellow
            np.savetxt("deloye_y2.dat", np.column_stack((np.append(x1,x2),np.append(y1,y2))))
            x1,y1 = plot(line13, fmt, **kwargs)    #red
            x2,y2 = plot(line14, fmt, **kwargs)    #red
            np.savetxt("deloye_r2.dat", np.column_stack((np.append(x1,x2),np.append(y1,y2))))
            ##plot(line15,fmt)  #??
            x1,y1 = plot(line16, fmt, **kwargs)    #green
            x2,y2 = plot(line17, fmt, **kwargs)    #green
            np.savetxt("deloye_g2.dat", np.column_stack((np.append(x1,x2),np.append(y1,y2))))
            x1,y1 = plot(line18, fmt, **kwargs)    #blue
            x2,y2 = plot(line19, fmt, **kwargs)    #blue
            np.savetxt("deloye_b2.dat", np.column_stack((np.append(x1,x2),np.append(y1,y2))))
    

#### Functions for Yungleson model lines from paper


def readCsv(fname):
    x, y = np.loadtxt(fname, unpack=True, delimiter=',')
    return x,y
    
def plotYung(fname, fmt, *args, **kwargs):
    x, y = readCsv(fname)
    return plt.plot(x,y,fmt, *args, **kwargs)


### Functions for Goliasch and Nelson CV track

def plotGol(fname, *args, **kwargs):
    lr, lm = np.loadtxt(fname, unpack=True)
    return plt.plot(lm, lr, *args, **kwargs)




### Functions for other lines

def radZS(m):
    """ Cold sphere radii from Zapolsky and Salpeter 1969, valid for 0.002 < M < 0.45, accurate to w/in 3%.
    """
    return 0.0106 - 0.0064 * np.log(m) + 0.0015 * m**2

def radTF(m):
    """ Helium secondary m-r relation, from Tutukov and Fedorova 1989.
    """
    return 0.043 * m**(-0.062)

def radYung(m, porb0, m0):
    """ Linear approx. of helium secondaries, from Yungelson 2008, for init M2 between 0.34 and 0.4 and init Porb of 10-35 minutes.
        porb0 and m0 are init Porb and M2.
    """
    return 10**(-1.478) * (porb0 / 20.)**(-0.05) * m**(-0.16) * (0.35 / m0)**(0.345)
    
def radSav(m):
    """ Semi-degenerate He secondaries, from Savonije et al 1986
    """
    return 0.029 * m ** (-0.19)


### For reading / writing to files

def loadFromCsv(loadname, title=None, usecols=None, delimiter=';'):
    """ Load a whole table from a csv file. Should be fairly general
    """
    with open(loadname) as f:
        lines = f.read().splitlines()
    
    table = {"Title":title, "Length":len(lines)-1} 
    
    
    headers = lines[0].split(';')
    
    for j, header in enumerate(headers):
        # Find values in table and make array
        if usecols is not None and not j in usecols:
            continue
        arr = np.array([],dtype=str)
        for line in lines[1:]:
            vals = line.split(';')
            arr = np.append(arr, vals[j])
        # Add array to dictionary
        table[header] = arr
    return table

def makeFloat(array, verbose=False):
    """ Takes an array of string numbers and returns them as floats. Problem numbers are returned as Nans
    """
    vals = np.array([])
    for val in array:
        try:
            newval = float(val)
        except ValueError:
            if verbose:
                print ("Did not recognise:", val)
            newval = np.nan
        vals = np.append(vals, newval)
    if len(vals) == 1:
        return vals[0]
    return vals
    


def reference(source, useShorthandReferences=False, shorts=None):
    """ Tidy up the references -- given references of form Wevers2016, will turn into form you want to appear as in paper (eg. \citep{Wevers2016})
    """
    
    def applyShorthand(ret):
        ## See if source is already in the dictionary of shorthands. If it is, return number. Else, add and then return number
        if not ret in shorts:
            shorts[ret] = str(len(shorts) + 1)
        return shorts[ret]
    
    ## Checks, apply formatting
    if "?" in source:
        return ''
        
    else: 
        ret = ()
        for s in source.split(','):
            s = s.strip()
            if "vsnet" in s or "in prep" in s or "in press" in s or "in review" in s or "priv. comm." in s \
                                or "private communication" in s or "this work" in s.lower():
                if useShorthandReferences:
                    ret += (applyShorthand(s),)
                else:
                    ret += (s,)
            else:
                if useShorthandReferences:
                    ret += (applyShorthand("\citet{%s}"%(s)),)
                else:
                    ret += ("\citet{%s}"%(s),)
        return ",".join(ret)
    

def shorthandSummary(shorts):
    """ Produce a key telling reader which reference belongs to which number
    """
    key = []
    nums = []
    # turn the 'shorts' dict into a list of strings in correct format
    for j, (source, num) in enumerate(shorts.iteritems()):
        key += ("[%s]~%s"%(num, source),)
        nums += [int(num)]
    # sort strings -- this will work as long as tehy are of form 1: source, 2: ...
    key = np.array(key)[np.argsort(nums)]
    # combine all strings into one big string
    return "; ".join(key) + "."
        

def addComment(comment,comdict):
    """ Add a comment to a dictionary of comments, which will then be added to the footer
        Returns a latex-style superscript letter
    """
    if comment == None:
        return ''
    else:
        if comment not in comdict:
            comdict[comment] = len(comdict)
        return "$^{%s}$"%(string.ascii_lowercase[comdict[comment]])
    
    
    
def commentSummary(comdict):
    """ Produce a key of comments which could go in a table footer
    """
    if len(comdict) == 0:
        return ''
    
    key = []
    nums = []
    for j, (comment, num) in enumerate(comdict.iteritems()):
        key += (["$^{%s}$ %s"%(string.ascii_lowercase[num], comment)])
        nums += [int(num)]
    key = np.array(key)[np.argsort(nums)]
    return " \\newline ".join(key)
    
    
### Functions for plotting known AM CVns 





class Star:
    def __init__(self, name, period, mmin=None, mmax=None, xoffset=None, yoffset=None, endpoint=None, diagonoffset=None):
        self.name = name
        self.period = period * 60
        self.mmin = mmin
        self.mmax = mmax
        self.endpoint = endpoint
        if xoffset is not None:
            self.xoffset = xoffset
        else:
            self.xoffset = 0.0
        if yoffset is not None:
            self.yoffset = yoffset
        else:
            self.yoffset = -0.02
        if diagonoffset is not None:
            self.diagonoffset = diagonoffset
        else:
            self.diagonoffset = -0.2
        
    

def cvRad(period, mass):
    """ Given a period and mass or array of masses, returns corresponding radii for a CV with those properties.
        M/R from Solheim 2010 eq8, combination of Paczynski 1971's Roche-lobe approximation and Kepler's third. Valid for q<0.8.
    """
    return np.power((0.01**3 * (mass / 0.1) * np.power((period / 101.), 2)), 1./3)


def cvMassRadius(period, n=200, minmasslog=-2.5, maxmasslog=-0.1):
    """ Returns masses and radii corresponding to a CV or AM CVn of a certain period. Wrapper around cvRad.
    """
    mass = np.logspace(minmasslog, maxmasslog, num=n)
    rad = cvRad(period, mass)
    return mass, rad

def eggletonRoche(q,a):
    """ Find the volume-equivalent Roche lobe radius for a CV donor, given mass ratio and orbital separation
    """
    return (a*0.49*q**(2./3)) / (0.6*q**(2./3)+np.log(1.+q**(1./3)))
    
def eggletonMR(m1,m2,period):
    """ Returns donor radii calculated using the Eggleton Roche lobe radius formula. Period in seconds. Wrapper around eggletonRoche
    """
    a = orbital_separation(m1,m2,period/86400.)
    q = m2/m1
    return eggletonRoche(q,a)

def findEps(p,sh):
    return (sh-p)/p


def orbital_separation(m1, m2, period):
    """
    Returns orbital semi-major axis in solar radii given masses in solar masses
    and an orbital period in days, i.e. Kepler's third law. Works with scalar
    or numpy array inputs.  ##Taken from Tom's subs code
    """

    return (G*mSun*(m1+m2)*(secondsPerDay*period/(2.*np.pi))**2)**(1./3.)/rSun



fmts = {"Eclipses":'ro',"Spectroscopy":'gs',"Superhumps":'bo',}
def plotStar(table,name,fmt=None,markersize=None,label=None, defaultM1=0.7, redoSuperhumps=False, verbose=False, notex=False, overrideM2=None, sname=sname_dat, *args, **kwargs):
    """ Plot a known amcvn on the graph
    """
    oi = (table["Name"]==name)
    
    # Default formatting
    if fmt==None:
        try:
            fmt = fmts[ table["q Method"][oi][0] ]
        except KeyError:
            fmt = 'k.'
            if verbose:
                print ("Unable to find format for", name)
                
    if 'marker' in kwargs:          # handle case where the user is trying to override the default marker shape
        marker = kwargs.pop('marker')
    else:
        marker = fmt[1]
        
    if markersize == None:
        if marker == '*':
            markersize=18
        else:
            markersize=12
    

    
    # Read table
    p = makeFloat(table["Period (min)"][oi])#[0].replace('(sh)','').replace('?',''))
    q = makeFloat(table["q"][oi])
    qerr = makeFloat(table["q err"][oi])
    sh = makeFloat(table["Superhump Excess"][oi])
    sherr = makeFloat(table["Superhump err"][oi])
    m1 = makeFloat(table["M1"][oi])
    m1err = makeFloat(table["M1 err"][oi])
    m2 = makeFloat(table["M2"][oi])
    m2err = makeFloat(table["M2 err"][oi])
    
    if overrideM2 is not None:
        m2 = overrideM2
        m2err = 0
    
    ## Recalculate q from superhumps?
    if redoSuperhumps and not (redoSuperhumps != "strong" and table["q Method"][oi][0] != "Superhumps") and (not np.isnan(sh)):
        
        q, qerr = findQKnigge(sh, sherr)             ##!!! Change this once McAllister's paper is published
        print(q,qerr)
    
    if textable and not notex:
        with open(sname, 'a',) as f:
            
            def checknan(val, err, prec=3):
                if np.isnan(val):
                    return "--"
                else:
                    if prec==3:
                        return "$%.3f \\pm %.3f$"%(val, err)
                    elif prec==2:
                        return "$%.2f \\pm %.2f$"%(val, err)
                    else:
                        raise KeyError("Value of prec not recognised")
                    
            if name in comments:
                comment = comments[name]
            else:
                comment = None
            
            
            f.write("%s%s & %s & %s & %s & %s & %s\\\\ \n"%(name, addComment(comment, comdict), checknan(sh,sherr),checknan(q,qerr),checknan(m1,m1err,2), table["q Method"][oi][0], reference(table["q Source"][oi][0], True, shorts)))
    
    
    ## Missing value handling
    if np.isnan(m2) and (not np.isnan(m1)) and (not np.isnan(q)):
        m2, m2err = functApp(np.multiply, m1, m1err, q, qerr)
    elif np.isnan(m2) and np.isnan(m1) and (not np.isnan(q)):
        m1 = defaultM1
        m1err = 0.1
        m2, m2err = functApp(np.multiply, m1, m1err, q, qerr)
        print(m2,m2err)
    elif (not np.isnan(m2)) and np.isnan(m1) and (not np.isnan(q)):
        m1 = m2 / q
    
    
    
    ## Need period, m1 and m2
    p = p * 60.
    
    print(f"Plotting {name} with m2 = {m2} +/- {m2err}")
    
    plt.plot(np.log10([m2-m2err,m2+m2err]), np.log10(eggletonMR(m1,np.array([m2-m2err,m2+m2err]),p)), fmt[0]+'-', linewidth=3, *args, **kwargs)  #Error bar line
    plt.plot(np.log10(m2), np.log10(eggletonMR(m1,m2,p)), fmt, markersize=markersize, label=label, marker=marker, *args, **kwargs)    # Marker
    
    
    
    ## write data to a file
    if sname is not None:
        name_fixed = name.replace(' ','_')
        r2,r2err = functApp(eggletonMR, m1,m1err, m2,m2err, p,0)
        with open(sname, 'a') as f:
            f.write(f"{name_fixed} {p:.0f} {q:.4f} {qerr:.4f} {m1:.4f} {m1err:.4f} {m2:.4f} {m2err:.4f} {r2:.4f} {r2err:.4f}\n")



def plotStars(table, *args, **kwargs):
    
    for name in table["Name"]:
        plotStar(table,name, *args, **kwargs)

def plotStarManually(p,q,qerr,fmt=None,markersize=None,label=None, defaultM1=0.7, redoSuperhumps=False, verbose=False, notex=False, overrideM2=None, sname=sname_dat, *args, **kwargs):
    """ Add an object manually rather than reading from the table
    """
    
    # Default formatting
    if fmt==None:
        try:
            fmt = fmts[ table["q Method"][oi][0] ]
        except KeyError:
            fmt = 'k.'
            if verbose:
                print ("Unable to find format for", name)
                
    if 'marker' in kwargs:          # handle case where the user is trying to override the default marker shape
        marker = kwargs.pop('marker')
    else:
        marker = fmt[1]
        
    if markersize == None:
        if marker == '*':
            markersize=18
        else:
            markersize=12
            
    m1 = defaultM1
    m1err = 0.1
    m2, m2err = functApp(np.multiply, m1, m1err, q, qerr)
    
    
    ## Need period, m1 and m2
    p = p * 60.
    
    print(f"Plotting manual target with m2 = {m2} +/- {m2err}")
    
    plt.plot(np.log10([m2-m2err,m2+m2err]), np.log10(eggletonMR(m1,np.array([m2-m2err,m2+m2err]),p)), fmt[0]+'-', linewidth=3, *args, **kwargs)  #Error bar line
    plt.plot(np.log10(m2), np.log10(eggletonMR(m1,m2,p)), fmt, markersize=markersize, label=label, marker=marker, *args, **kwargs)    # Marker
    
    



def quad(*args):
    """ Combine a list of arguments in quadrature
    """
    args = np.array(args)
    return np.sqrt(np.sum(args**2))


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






#### Finding Q equations


def findQ(eps):
    """ Inverse of Patterson 2005 eqn
    """
    return quadraticFormula(0.29,0.18,-eps)[0]

def findEpsKato(q):
    """ Empirical relation from Kato 20??
    """
    epsAsterisk = 0.00027 + 0.402*q - 0.467*q**2 + 0.297*q**3
    eps = epsAsterisk / (1 - epsAsterisk)
    return eps


def findQKnigge(eps,epserr):
    """ Knigge 2006's eps-q relation (with errors)
    """
    def calculate(a,b,c,eps):
        return a + b*(eps-c)
    return functApp(calculate,0.114,0.005,3.97,0.41,0.025,0,eps,epserr,)

def findQMcAllisterB(eps,epserr):
    """ McAllister 2018's update to Knigge's eps-q relation (with errors) for phase B superhumps
    """
    def calculate(a,b,c,eps):
        return a + b*(eps-c)
    return functApp(calculate,0.118,0.003,4.45,0.28,0.025,0,eps,epserr,)
    
def findQMcAllisterC(eps,epserr):
    """ McAllister 2018's update to Knigge's eps-q relation (with errors) for phase C superhumps
    """
    def calculate(a,b,c,eps):
        return a + b*(eps-c)
    return functApp(calculate,0.135,0.004,5.0,0.7,0.025,0,eps,epserr,)

