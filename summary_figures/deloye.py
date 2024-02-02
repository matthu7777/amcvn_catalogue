#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sys import argv
import string

import mgutils as mg
from mgutils import sci
from deloyelines import *

""" Plot various models of amcvn donors against the known systems. 
    Also some functions that are useful for converting from the xfig of deloye to numpy arrays. Could have some things useful for general xfig reading.
"""
save = "nosave" not in argv

figsize = (10,12)
xlim = (-2.2, -0.3)
ylim = (-1.65, -0.75)
plot_stars = True
plot_hecvs = False
legend = True
annotate = False
periodDiags = False
twoDeloye = False        # Plot both sets of Deloye tracks (ie different starting masses)
textable = False         # Write a latex table of q values to a file
poster = True           # Use mgposter rather than mgplot (makes text bigger, but if you want to increase the legend text size you need to change it yourself)
sname = "latex_table.txt"       # Save name for latex table

if len(argv) > 1 and argv[1] == "change":
    # allows to overwrite gaia mass from command line
    change = True
    gaiaMass = float(argv[2])
    gaiaErr = float(argv[3])
else:
    change = False
    
comments = {\
"HP Lib":"Additional constraints on the donor masses in these systems were derived based on the magnitude and parallax measurements",
"CR Boo":"Additional constraints on the donor masses in these systems were derived based on the magnitude and parallax measurements",
"V803 Cen":"Additional constraints on the donor masses in these systems were derived based on the magnitude and parallax measurements",
"YZ LMi":"The tabulated $q$ for YZ\\,LMi is derived from eclipse photometry; its superhump excess gives $q = 0.049 \pm 0.008$ by Equation~\\ref{eq:mrknigge}",
"AM CVn":"The tabulated $q$ for AM\\,CVn is derived from spectroscopy; its superhump excess gives $q = 101 \pm 0.005$ by Equation~\\ref{eq:mrknigge}",
}


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


def plotDeloyeModels(fmt=None, label=None, **kwargs):
    if fmt is None:
        plot(line1, 'y-', **kwargs)     #yellow
        plot(line2, 'y-', **kwargs)     #yellow
        plot(line3, 'r-', **kwargs)     #red
        plot(line4, 'r-', **kwargs)     #red
        #plot(line5,'rD')   #??
        plot(line6, 'g-', **kwargs)     #green
        plot(line7, 'g-', **kwargs)     #green
        plot(line8, 'b-', **kwargs)     #blue
        plot(line9, 'b-', **kwargs)     #blue
        #plot(line10,'rD')  #??
        plot(line11, 'y-', **kwargs)    #yellow
        plot(line12, 'y-', **kwargs)    #yellow
        plot(line13, 'r-', **kwargs)    #red
        plot(line14, 'r-', **kwargs)    #red
        #plot(line15,'rD')  #??
        plot(line16, 'g-', **kwargs)    #green
        plot(line17, 'g-', **kwargs)    #green
        plot(line18, 'b-', **kwargs)    #blue
        plot(line19, 'b-', **kwargs)    #blue
    else:
        plot(line1, fmt, label=label, **kwargs)        #yellow
        plot(line2, fmt, **kwargs)      #yellow
        plot(line3, fmt, **kwargs)      #red
        plot(line4, fmt, **kwargs)      #red
        #plot(line5,fmt)    #??
        plot(line6, fmt, **kwargs)      #green
        plot(line7, fmt, **kwargs)      #green
        plot(line8, fmt, **kwargs)      #blue
        plot(line9, fmt, **kwargs)      #blue
        if twoDeloye:
            #plot(line10,fmt)   #??
            plot(line11, fmt, **kwargs)    #yellow
            plot(line12, fmt, **kwargs)    #yellow
            plot(line13, fmt, **kwargs)    #red
            plot(line14, fmt, **kwargs)    #red
            ##plot(line15,fmt)  #??
            plot(line16, fmt, **kwargs)    #green
            plot(line17, fmt, **kwargs)    #green
            plot(line18, fmt, **kwargs)    #blue
            plot(line19, fmt, **kwargs)    #blue
    

#### Functions for Yungleson model lines from paper


def readCsv(fname):
    x, y = np.loadtxt(fname, unpack=True, delimiter=',')
    return x,y
    
def plotYung(fname, fmt, *args, **kwargs):
    x, y = readCsv(fname)
    return plt.plot(x,y,fmt, *args, **kwargs)


### Functions for Goliasch and Nelson CV track

def plotGol(fname, *args, **kwargs):
    lm, lr = np.loadtxt(fname, unpack=True)
    return plt.plot(lr, lm, *args, **kwargs)




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
    a = sci.orbital_separation(m1,m2,period/86400.)
    q = m2/m1
    return eggletonRoche(q,a)

def findEps(p,sh):
    return (sh-p)/p


fmts = {"Eclipses":'ro',"Spectroscopy":'gs',"Superhumps":'bo',}
def plotStar(table,name,fmt=None,markersize=None,label=None, defaultM1=0.7, redoSuperhumps=False, verbose=False, notex=False, overrideM2=None, *args, **kwargs):
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
                
    if markersize == None:
        if '*' in fmt:
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
        
        #q, qerr = mg.findQMcAllisterB(sh, sherr)
        q, qerr = sci.findQKnigge(sh, sherr)             ##!!! Change this once McAllister's paper is published
        print (q)
    
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
        m2, m2err = mg.functApp(np.multiply, m1, m1err, q, qerr)
    elif np.isnan(m2) and np.isnan(m1) and (not np.isnan(q)):
        m1 = defaultM1
        m1err = 0.1
        m2, m2err = mg.functApp(np.multiply, m1, m1err, q, qerr)
    elif (not np.isnan(m2)) and np.isnan(m1) and (not np.isnan(q)):
        m1 = m2 / q
    
    ## Need period, m1 and m2
    p = p * 60.
    plt.plot(np.log10([m2-m2err,m2+m2err]), np.log10(eggletonMR(m1,np.array([m2-m2err,m2+m2err]),p)), fmt[0]+'-', linewidth=3, *args, **kwargs)  #Error bar line
    plt.plot(np.log10(m2), np.log10(eggletonMR(m1,m2,p)), fmt, markersize=markersize, label=label, *args, **kwargs)    # Marker
    
    print(name, m1,m1err, m2, m2err)



def plotStars(table, *args, **kwargs):
    
    for name in table["Name"]:
        plotStar(table,name, *args, **kwargs)
    

shorts = {}
comdict = {}
tableWidth = "12.5cm"

if __name__ == "__main__":
    ## Format graphs
    if not legend:
        mg.formatGraph(xlabel="log (M$_2$ / M$_\odot$)", ylabel="log (R$_2$ / R$_\odot$)", figsize=figsize, poster=poster)
        ax = plt.gca()
        # Defaults
        ax.set_ylim(*ylim)
        ax.set_xlim(*xlim)
    else:
        # Need larger figure
        mg.formatGraph(xlabel="log (M$_2$ / M$_\odot$)", ylabel="log (R$_2$ / R$_\odot$)", figsize=figsize, poster=poster)
        ax = plt.gca()
        ax.set_ylim(*ylim)
        ax.set_xlim(*xlim)
        
        
        
    if textable:
        ### Header
        with open(sname, 'w') as f:
            f.write("\\begin{table*} \n \\centering \n \\caption{A summary of the AM\\,CVn mass ratios used for Figure~\\ref{fig:mrdiag-j1351}. Where $q$ was derived by the superhump method, we recalculate it using Equation~\\ref{eq:mrknigge} for the sake of consistency.} \n")
            f.write("\\label{tab:masses-j1351} \n \\begin{tabular}{lccccc} \n ")
            f.write("\\hline \n Designation & $\epsilon$ & $q$ & $M_2$ ($M_\\odot$) & Method & Reference \\\\ \n  &&&&& \\\\ \n \\hline \n")
    
    
    ## Plot period diagonal lines if desired
    if periodDiags:
        colour = 'r'
        for p in [10,20,30,40,50,60]:
            m, r = cvMassRadius(p*60.)
            plt.plot(np.log10(m),np.log10(r),color=colour,)
        plt.annotate("10", (-1.23,-1.58),color=colour,)
        plt.annotate("20", (-1.62,-1.51),color=colour,)
        plt.annotate("30", (-1.85,-1.47),color=colour,)
        plt.annotate("40", (-2.01,-1.44),color=colour,)
        plt.annotate("50", (-2.16,-1.42),color=colour)
    
    
    ### plot Deloye models
    plotDeloyeModels(fmt='k-', label="White dwarf models (Deloye et al. 2007)")
    
    
    ### plot Yungelson models
    plotYung("yung_blue.csv",'k--', label="Helium-star models (Yungelson, 2008)")
    plotYung("yung_red.csv",'k--')
    plotYung("yung_cyan.csv",'k--')
    plotYung("yung_green.csv",'k--')
    plotYung("yung_black.csv",'k--')
    
    ### plot Goliasch model
    plotGol("goliasch.txt", 'k-.', label="Hydrogen-CV models (Goliasch \& Nelson, 2015)")
    plotGol("Goliasch0_85.dat", 'k:', label="Hydrogen-CV models (Goliasch \& Nelson, priv. comm.)")
    
    
    
    ## Plot other model lines
    mrange = np.logspace(-3,0,200)
    plt.plot(np.log10(mrange), np.log10(radZS(mrange)), 'k-', linewidth=2)  # Cold sphere radii
    #plt.plot(np.log10(mrange), np.log10(radTF(mrange)), 'k:')  # Nelemans Helium secondaries radii
    #plt.fill_between(np.log10(mrange), np.log10(radYung(mrange, 49, 0.4)), np.log10(radYung(mrange, 49, 0.34)), alpha=0.5, color='k', label="Helium models (Yungelson, 2008)")
    #plt.plot(np.log10(mrange), np.log10(radYung(mrange, 49, 0.4)), 'k--')
    #plt.plot(np.log10(mrange), np.log10(radYung(mrange, 49, 0.34)), 'k--')
    #plt.plot(np.log10(0.0234),np.log10(0.0594),'ko')   #Closest Goliasch model to what i originally sent him
    #plt.plot(np.log10(mrange), np.log10(radSav(mrange)), 'k:') # Savonije He radii
    
    
    
    
    
    
    
    
    
    ### Plot all star masses and mass ranges
    table = loadFromCsv("amcvns.csv")
#    plotStars(table, redoSuperhumps=True)

    if plot_stars and not plot_hecvs:
        
        plotStar(table, "SDSSJ1351-0643", color='grey')
        plotStar(table, "AM CVn", redoSuperhumps=True, fmt='gs')
        plotStar(table, "HP Lib", redoSuperhumps=True, color='grey')
        plotStar(table, "CXOGBS J1751-2940", redoSuperhumps=True, color='grey')
        
        plotStar(table, "CR Boo", redoSuperhumps=True, label="Superhump-derived $q$", color='grey')
        plotStar(table, "KL Dra", redoSuperhumps=True, color='grey')
        plotStar(table, "V803 Cen", redoSuperhumps=True, color='grey')
        plotStar(table, "CP Eri", redoSuperhumps=True, color='grey')
        plotStar(table, "NSV 1440", redoSuperhumps=False, color='grey')
        plotStar(table, "SDSSJ1240-0159", redoSuperhumps=True, fmt='gs', label="Spectroscopy-derived $q$")    #spec
        plotStar(table, "SDSSJ0129+3842", redoSuperhumps=True, color='grey')
        plotStar(table, "GP Com", redoSuperhumps=True, fmt='gs')            #spec
        plotStar(table, "SDSSJ0902+3819", redoSuperhumps=True, color='grey')
        plotStar(table, "V396 Hya", redoSuperhumps=True, fmt='gs')          #spec
        plotStar(table, "Gaia14aae", redoSuperhumps=True, fmt='b*', label="Gaia14aae (Green et al. 2018)")         #eclipse
        
        ####yz lmi has to go last or it gets covered by one of the superhump systems
        plotStar(table, "YZ LMi", redoSuperhumps=True, fmt='b*', label="YZ LMi (Copperwheat et al. 2011)", notex=True)      #eclipse
        #plotStar(table, "YZ LMi", redoSuperhumps=True, fmt='b*', label="SDSS J0926+3624 (Copperwheat et al. 2011)", notex=True)      #eclipse
        
        
    if plot_stars and plot_hecvs:
        
        plotStar(table, "SDSSJ1351-0643", color='grey', fmt='ko')
        plotStar(table, "AM CVn", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "HP Lib", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "CXOGBS J1751-2940", redoSuperhumps=True, color='grey', fmt='ko')
        
        plotStar(table, "CR Boo", redoSuperhumps=True, label="AM CVn", color='grey', fmt='ko')
        plotStar(table, "KL Dra", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "V803 Cen", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "CP Eri", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "NSV 1440", redoSuperhumps=False, color='grey', fmt='ko')
        plotStar(table, "SDSSJ1240-0159", redoSuperhumps=True, color='grey',fmt='ko')
        plotStar(table, "SDSSJ0129+3842", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "GP Com", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "SDSSJ0902+3819", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "V396 Hya", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "V406 Hya", redoSuperhumps=True, color='grey', fmt='ko')
        plotStar(table, "Gaia14aae", redoSuperhumps=True, color='grey', fmt='ko')
        
        ####yz lmi has to go last or it gets covered by one of the superhump systems
        plotStar(table, "YZ LMi", redoSuperhumps=True, color='grey', fmt='ko')      #eclipse
        
        #### He CVs
        plotStar(table, "CRTS J112253.3-111037", redoSuperhumps=True, fmt='rD', label="He CVs")
        plotStar(table, "CRTS J111126.9+571239", redoSuperhumps=True, fmt='rD')
        plotStar(table, "V485 Cen", redoSuperhumps=True, fmt='rD')
        plotStar(table, "EI Psc", redoSuperhumps=True, fmt='rD')
        
    
    if textable:
        with open(sname, 'a') as f:
            ### Footer
            f.write("\\hline \n \\multicolumn{6}{p{%s}}{References: %s}\\\\ \n"%(tableWidth, shorthandSummary(shorts)))
            
            if len(comdict) != 0:
                f.write("\\hline \n \\multicolumn{6}{p{%s}}{%s}\\\\ \n"%(tableWidth, commentSummary(comdict)))
            
            f.write("\\hline \n \\end{tabular} \n \\end{table*}")
    
    
    if legend:
        ###Legend
        handles, labels = ax.get_legend_handles_labels()
        handles, labels = np.array(handles), np.array(labels)
        print ("Legend labels:")
        for j, l in enumerate(labels):
            print (j, l)
        try:
            order = [4,5,0,1,2,3]
            ax.legend(handles[order], labels[order], loc=3, fontsize=14, numpoints=1, framealpha=1)
        except IndexError:
            order = np.arange(len(labels))
            ax.legend(handles[order], labels[order], loc=3, fontsize=14, numpoints=1, framealpha=1)
            
        print ("Current order:", order)
    
    if annotate:
        #plt.annotate("J1351",(-1.1,-1.42), color='r')
        #plt.annotate("AM CVn",(-0.9,-1.33), color='g')
        #plt.annotate("HP Lib",(-1.14,-1.37), color='grey')
        #plt.annotate("CXOGBS J1751",(-1.5,-1.38), color='grey')
        #plt.annotate("SDSS J0926+3624",(-1.45,-1.38), color='b')
        plt.annotate("Gaia14aae",(-1.55,-1.24), color='r')
    
    
    plt.tight_layout()
    if save:
        plt.savefig("M-R_plot.pdf")
        plt.savefig("M-R_plot.png")
        
    plt.show()
    
