#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import sys, os, json
import argparse, subprocess
from scipy.stats import linregress


def parse_highlights(highstr,l):
    if highstr=="":
        return []
    hs=highstr.split(" ")
    highs=[]
    for i,v in enumerate(hs):
        if "+" in v:
            s1,s2=v.split("+")
            highs.append(int(s1)+int(s2))
            continue
        highs.append(int(v))
    r=list(range(l))
    h=1.0
    for i,v in enumerate(r):
        if i+1 in highs:
            r[i]=h
            h*=-1.0
        else:
            r[i]=0
    return r


def add_highlights(prob_map, highs):
    ret=[]
    if len(highs)==0:
        return prob_map
    for i,v in enumerate(prob_map):
     #   print("ori", len(v),highs[i])
        ret.append(list(v)+[highs[i]])
    #    print("final",len(ret[-1])) ############
    ret.append(highs+[1.0])
   # print(len(ret[-1])) ############
    return np.array(ret)


#The size is given in inches for this function, but the program itelf takes cm
def plot_maps(prob_map,name, mini=None,maxi=None,tickfreq=0,title="",odpi=600,size=4, noshow=False, highlights=""):
    fig, ax = plt.subplots(1, 1, figsize=(size,size))
    #We ensure that the range is symmetrical around 0
    print(mini, maxi) ##############################
    highs=parse_highlights(highlights,len(prob_map))
    prob_map=add_highlights(prob_map,highs)
    im = ax.imshow(prob_map,norm=None,cmap='bwr',vmin=mini,vmax=maxi) # interpolation=None
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Color Bar", rotation=-90, va="bottom")
    
    ax.set_title(title)
    if tickfreq>0:
        if tickfreq<50:
            ax.tick_params(axis='x', labelsize=4, labelrotation = 45)
            ax.tick_params(axis='y', labelsize=5)
        loc = plticker.MultipleLocator(base=tickfreq) # From stackoverflow 
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)

    

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
p.add_argument("cols",type=int, help="Columns in the file, not counting the first one")
p.add_argument("--dirname",type=str, help="Name of the directory. If empty, the name is taken from the input file",default="")

p.add_argument("--filter",type=float, help="filter (set to zero) any value with absolute value smaller than the given",default=0.0)

p.add_argument("--exec",type=str, help="Full path to the dep executable",default="$GOMD/dep/dep")
p.add_argument("--dpi",type=int, help="Resolution (in dots per inch) of the plot images",default=600)


p.add_argument("--filterd",type=float, help="Only print values with absolute value larger than this",default=0.1)
p.add_argument("--di",type=int, help="Delay for the first figure",default=0)
p.add_argument("--df",type=int, help="Delay for the last figure",default=0)
p.add_argument("--deltad",type=int, help="Delay step between figures",default=0)
p.add_argument("--size",type=int, help="Size for the figures, in cm. The figure is always a square",default=10)
p.add_argument("--chunklen",type=int, help="If > 0, the trajectory is taken as formed by independent 'chunks' of chunklen frames each",default=0)

p.add_argument("--blur",type=int, help="Use not the i-delay element, but the average between i-delay-(blur/2) and i-delay+(blur/2). If an odd number is given, the previous even number will be used dinstead",default=10)

p.add_argument("--tickfreq",type=int, help="if >0, sets the tick frequency in the plot to its value",default=0)

p.add_argument("--noshow",type=bool, help="Don't show the plots, just save the images",default=False)
p.add_argument("--novideo",type=bool, help="Do not call ffmpeg to make a video out of the png files",default=False)

p.add_argument("--skipcalc",type=bool, help="Don't obtain the dependence matrix, read it from a previously written file",default=False)
p.add_argument("--extraplot",type=bool, help="Plot p-values and standard errors, if available",default=False)
p.add_argument("--corr",type=bool, help="Obtains correlations instead of slopes",default=False)

p.add_argument("--det1", type=str, help="Print the correlations involving the given residues, and those given in --det2, which must be given separated by spaces, surounded by quotes o double quotes",default="")
p.add_argument("--det2", type=str, help="Print the correlations involving the given residues, and those given in --det2, which must be given separated by spaces, surounded by quotes o double quotes",default="")
p.add_argument("--highlights", type=str, help="A string, enclosed by \" \" and separates with spaces, of the residues to highlight in the final plot, or empty if nothing is to be highlighted",default="")


p.add_argument("--min",type=float, help="Minimum in the data range, for plotting",default=1)
p.add_argument("--max",type=float, help="Maximum in the data range, for plotting",default=-1)
p.add_argument("--pvalfilter",type=float, help="points with p-values lareger than this will be set to 1",default=0.01)


a = p.parse_args()

corr=""
if a.corr:
    corr="corr"


oridir=os.getcwd()
newdir=a.finp.replace(".dat","")+corr
if a.dirname:
    newdir=a.dirname
try:
    os.mkdir(newdir)
    if not a.skipcalc:
        os.system("cp %s %s/GOMDDATA.dat"%(a.finp,newdir))
except:
    pass
os.chdir(newdir)



cont=0
for i in range(a.di,a.df,a.deltad):
    name="%04d.json"%(cont) #we use a counter and not the delay to make it easier to 
    #transform all the PNGs in a video. The program will print each delay anyway
    print("Delay: %d  Filename: %s"%(i,name))
    cont+=1
    corr=""
    if a.corr:
        corr="--corr"
    if not a.skipcalc:
        print("%s --delay %d --delayblur %s %s -out %s -c %d --chunklen %d GOMDDATA.dat"%(a.exec,i,a.blur,corr,name, a.cols, a.chunklen))
        os.system("%s --delay %d --delayblur %s %s -out %s -c %d --chunklen %d GOMDDATA.dat"%(a.exec,i,a.blur,corr,name, a.cols, a.chunklen))
    fin=open(name,"r")
    resmatrix=np.array(json.load(fin))
    fin.close()
    filtername=""
    if a.filter>0:
        for k,v in enumerate(resmatrix):
            for j,w in enumerate(v):
                if np.abs(w)<a.filter:
                    resmatrix[k][j]=0.0
                elif a.filter>0.0: #so, if we _are_ filtering, the ones that survived get promoted towards 1, to make them more
                    #noticeable
                    x=resmatrix[k][j]
                    resmatrix[k][j]*=(1+np.exp(-2*(x-1)))/2 #This is just a way of getting the unfiltered values towards 1
                    #It's just an inverted logistic-like function.
    if a.det1!="":
        d1=parsesel(a.det1)
        d2=parsesel(a.det2)
        for i,v in enumerate(resmatrix):
            for j,w in enumerate(v):
                if (i+1 in d1 and j+1 in d2) or (i+1 in d2 and j+1 in d1):
                    if np.abs(w)>=a.filterd:
                        print("i:%03d j:%03d %5.3f"%(i+1,j+1,w))
    if a.filter>0:
        filtername="filt%01.3f"%(a.filter)
    plotname=name.replace(".json","%s.png"%(filtername))
    title="Delay: %d Blur %d"%(i,a.blur)
    cm2inches=0.3937
    plot_maps(resmatrix,plotname,mini=a.min,maxi=a.max,tickfreq=a.tickfreq,title=title,odpi=a.dpi,size=a.size*cm2inches,noshow=a.noshow,highlights=a.highlights)

if not a.novideo:
    if a.filter:
        os.system("ffmpeg -framerate 2 -pattern_type glob  -i '*.png' -s 1920x1080 -c:v libx264 out_filt.mp4")
    else:
        os.system("ffmpeg -framerate 2 -i '%04d.png' -s 1920x1080 -c:v libx264 out.mp4")
os.chdir(oridir)
       




