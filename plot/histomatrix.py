#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as plticker
import sys, os, json
import argparse, subprocess
from scipy.stats import linregress

font = {'family' : 'computer modern',
       # 'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)


#The size is given in inches for this function, but the program itelf takes cm
def plot_maps(prob_map,name, mini=None,maxi=None,tickfreq=0,title="",odpi=600,size=4, noshow=False, highlights="", timed=False,printrow=-1):
    sizex=size
    sizey=size
    if printrow>=0:
        firstline=""
        secondline=""
        print("Printing row: %3d"%printrow)
        for i,v in enumerate(prob_map[printrow]):
            firstline+="%6d "%(i+1)
            secondline+="%6.3f "%(v)
        print(firstline)
        print(secondline)
    if timed:
        prop=len(prob_map)/len(prob_map[0])
        sizey=prop*size

    fig, ax = plt.subplots(1, 1, figsize=(sizex,sizey))
    #We ensure that the range is symmetrical around 0
    print(mini, maxi) ##############################
   # highs=parse_highlights(highlights,len(prob_map))
   # prob_map=add_highlights(prob_map,highs)
    im = ax.imshow(prob_map,norm=None,cmap='bwr',vmin=mini,vmax=maxi) # interpolation=None
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Color Bar", rotation=-90, va="bottom")
    title = "AD Distance histograms integration"
    ax.set_title(title)
    if tickfreq>0:
        if tickfreq<50:
            ax.tick_params(axis='x', labelsize=6, labelrotation = 45)
            ax.tick_params(axis='y', labelsize=6)
        loc = plticker.MultipleLocator(base=tickfreq) # From stackoverflow 
        ax.xaxis.set_major_locator(loc)
        if not timed:
            ax.yaxis.set_major_locator(loc)
    #SOD
    form1=list(range(1,152,10))
    form2=list(range(8,152,10))
    form=form1+form2
    form = plticker.FixedFormatter(form)
    loc2=list(range(0,306,10))
    print(form) #######
    loc = plticker.FixedLocator(loc2) # From stackoverflow 
    ax.xaxis.set_major_formatter(form)
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_formatter(form)
    ax.yaxis.set_major_locator(loc)
    #ENDSOD
    fig.tight_layout()
    plt.savefig(name, dpi=odpi)
    if not noshow:
        plt.show()


def parsesel(sele):
    fstr=sele.split(",")
    filt=[]
    for i in fstr:
        if "-" in i:
            ex=i.split("-")
            t=list(range(int(ex[0]),int(ex[1])+1))
            filt=filt+t
        else:
            filt.append(int(i))
    return filt






p = argparse.ArgumentParser()
p.add_argument('finp', metavar='gomd_output.dat', type=str, help='The data file in goMD format')
#p.add_argument("cols",type=int, help="Columns in the file, not counting the first one")
p.add_argument("--dirname",type=str, help="Name of the directory. If empty, the name is taken from the input file",default="")


p.add_argument("--lfilter",type=float, help="If >0, set to 0 anything with an absolute value smaller than this",default=-1)
p.add_argument("--ufilter",type=float, help="If >0,  set to 0 anything with an absolute value larger than this",default=-1)
p.add_argument("--size",type=int, help="Size for the figures, in cm. The figure is always a square",default=10)
p.add_argument("--dpi",type=int, help="Resolution (in dots per inch) of the plot images",default=600)



p.add_argument("--tickfreq",type=int, help="if >0, sets the tick frequency in the plot to its value",default=10)


p.add_argument("--noshow",type=bool, help="Don't show the plots, just save the images",default=False)
p.add_argument("--toone",type=bool, help="All values != 0 are taken to 1 or -1",default=False)

p.add_argument("--extraplot",type=bool, help="Plot p-values and standard errors, if available",default=False)
p.add_argument("--corr",type=bool, help="Obtains correlations instead of slopes",default=False)

p.add_argument("--highlights", type=str, help="A string, enclosed by \" \" and separates with spaces, of the residues to highlight in the final plot, or empty if nothing is to be highlighted",default="")
p.add_argument("--printrow", type=int, help="If >= 0 print the numerical values for the n-1th row ",default=-1)


p.add_argument("--min",type=float, help="Minimum in the data range, for plotting",default=0)
p.add_argument("--max",type=float, help="Maximum in the data range, for plotting",default=1)
p.add_argument("--sym",type=bool, help="Fill a symmetric matrix",default=False)

a = p.parse_args()

cm2inches=0.3937

fin=open(a.finp,"r")
resmatrix=np.array(json.load(fin))
        
for k,v in enumerate(resmatrix):
    for j,w in enumerate(v):
        if (a.lfilter>0 and np.abs(w)<a.lfilter) or (a.ufilter>0 and np.abs(w)>a.ufilter):
            resmatrix[k][j]=0.0
        if a.toone and  resmatrix[k][j]!=0:
            if resmatrix[k][j]<0:
                resmatrix[k][j]=-1
            else:
                resmatrix[k][j]=1


            

if a.sym:
    print("symmetric matrix")
    for i,v in enumerate(resmatrix):
        for j,w in enumerate(v):
            if resmatrix[i][j] == 0.0:
                print(i,j,resmatrix[i][j],resmatrix[j][i],"filled")
                resmatrix[i][j]=resmatrix[j][i]

plot_maps(resmatrix,"histogram",mini=a.min,maxi=a.max,tickfreq=a.tickfreq,title="histogram",odpi=a.dpi,size=a.size*cm2inches,noshow=a.noshow,highlights=a.highlights,printrow=a.printrow)

