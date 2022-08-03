#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from scipy.stats import linregress
import json

def plot_maps(prob_map,name, mini=None,maxi=None):
    fig, ax = plt.subplots(1, 1, figsize=(20,20))
    #We ensure that the range is symmetrical around 0
    if mini==None and maxi==None:
        maxi=np.amax(prob_map)
        mini=np.amin(prob_map)
        if np.abs(maxi)>=np.abs(mini):
            mini=-1*maxi
        else:
            maxi=-1*mini
    im = ax.imshow(prob_map,cmap='Spectral',interpolation=None,vmin=mini,vmax=maxi)
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Color Bar", rotation=-90, va="bottom")

    ax.set_title("Per-atom RMSD correlation")

    fig.tight_layout()
    plt.savefig(name, dpi=600)
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


#Winsize has to be smaller than len(vals[i])
#the whole reason for this function is the "blur" feature
#When blur is set to zero, makewindow(vals,delay,0) is 
#the same as simply writting vals[:len(vals)-delay]
def makewindow(vals,delay,blur):
    ret =[]
    for i,v in enumerate(vals):
        ret.append([])
        for j,w in enumerate(v):
            if j>=delay:
                ret[-1].append(np.mean(v[j-delay:j-delay+blur+1]))

    return np.array(ret)





p = argparse.ArgumentParser()
p.add_argument("-f", type=str, help="Input file, gomd format. Either this or the -jn option must be given",default="")
p.add_argument("-j", type=str, help="Read the already-computed correlation matrix from the JSON file given.", default="")
p.add_argument("--nomax1",type=bool, help="if slope >1 for the regression x vs y, the slope for regression y vs x. This option prevent that behavior. ",default=False)
p.add_argument("-n",type=bool, help="Normalize the RMSDs of each atom along the trajectory by dividing each value for the largest in the whole trajectory ",default=False)

p.add_argument("-c",type=int, help="Columns in the gomd files, not counting the time (i.e first) column",default=2)
p.add_argument("--filter",type=float, help="filter (set to zero) any value with absolute value smaller than the given",default=0.0)

p.add_argument("--det1", type=str, help="Print the correlations involving the given residues, and those given in --det2, which must be given separated by spaces, surounded by quotes o double quotes",default="")
p.add_argument("--det2", type=str, help="Print the correlations involving the given residues, and those given in --det2, which must be given separated by spaces, surounded by quotes o double quotes",default="")

p.add_argument("--filterd",type=float, help="If the --details flag is set, only print values with absolute value larger than this",default=0.1)

p.add_argument("--delay",type=int, help="if >0, obtain a delayed-correlation index with a delay of the given number of frames",default=0)
p.add_argument("--delayblur",type=int, help="if --delay is used, for the 'delayed' data, use not the i-delay element, but the average between i-delay and i-delay+delayblur",default=0)


p.add_argument("--noplot",type=bool, help="Don't plot the data",default=False)
p.add_argument("--extraplot",type=bool, help="Plot p-values and standard errors, if available",default=False)
p.add_argument("--checksymmetry",type=bool, help="Checks and prints warnings if the slope matrix is not symmetric",default=False)



p.add_argument("--min",type=float, help="Minimum in the data range, for plotting",default=None)
p.add_argument("--max",type=float, help="Maximum in the data range, for plotting",default=None)
p.add_argument("--pvalfilter",type=float, help="points with p-values lareger than this will be set to 1",default=0.01)



a = p.parse_args()

#ideas
#correlacion con "promedio ventana" para captar correlaciones con tiempos diferidos
#
# Indices root mean inter-sector correlation. 

vals=[]
resmatrix=[]
rvalmatrix=[]
pvalmatrix=[]

if a.f!="":
    for i in range(a.c):
        vals.append([])
        resmatrix.append([])
        pvalmatrix.append([])
        rvalmatrix.append([])
        for j in range(a.c):
            resmatrix[-1].append(-1000) #just any obviously wrong number
            pvalmatrix[-1].append(-1000) #just any obviously wrong number
            rvalmatrix[-1].append(-1000) #just any obviously wrong number
    fin=open(a.f,"r")
    for i in fin:
        l=i.split()
        for j,v in enumerate(vals):
            vals[j].append(float(l[j+1])) #We skip the first "frame" column, so we count columns from 1
    vals=np.array(vals)
    if a.n:
        for i,v in enumerate(vals):
            vals[i]=v/v.max()
    disp=[]
    if a.delay>0:
        disp=makewindow(vals,a.delay,a.delayblur)
    
##################
##DEBUG####3
#
#    x=vals[200]
#    y=vals[130]
#
#    print(x, np.std(x))
#    print(y, np.std(y))
#    
#    d = linregress(x, y)
 #   print(d.slope,d.rvalue)
 #   d2 = linregress(y,x)
 #   print(d2.slope,d2.rvalue)
 #   axes = plt.gca()
 #   axes.set_ylim(np.array([0,1]))
 #   axes.set_xlim(np.array([0,1]))
 #   plt.plot(x,y,"bo", alpha=0.1) 
 #   plt.plot(x,d.intercept+d.slope*x,"r")
 #   plt.savefig("xy.png")
 #   plt.show()
 #   axes = plt.gca()
 #   axes.set_ylim(np.array([0,1]))
 #   axes.set_xlim(np.array([0,1]))
 #
 #   plt.plot(y,x,"ro", alpha=0.1) 
 #   plt.plot(y,d2.intercept+d2.slope*y,"b")
 #   plt.savefig("yx.png")
 #
 #    plt.show()
 #    sys.exit(0)
##END DEBUG
    jsonending=".json"
    if a.delay<=0:
    #have the values prepared, now to build the correlation matrix
        for i,v in enumerate(vals):
            for j,w in enumerate(vals):
                stdv=np.std(v)
                stdw=np.std(w)
                if stdv>=stdw:
                    reg = linregress(v, w)
                else:
                    reg = linregress(w,v)
                if np.abs(reg.slope)>1.0 and not a.nomax1:
                    reg = linregress(w,v)
                pvalmatrix[i][j]=reg.pvalue
                rvalmatrix[i][j]=reg.rvalue
                resmatrix[i][j]=reg.slope
    else:
        for i,v in enumerate(vals):
            for j,w in enumerate(vals):
                y=w[a.delay:] #x is the delayed one, so I starts from a larger index to compensate
                x=disp[i]     #say, at the first position with delay 5 we will be comparing x[0] and y[4]
                stdx=np.std(x)
                stdy=np.std(y)
                if stdx>=stdy:
                    reg = linregress(x, y)
                else:
                    reg = linregress(y,x)
                if np.abs(reg.slope)>1.0 and not a.nomax1:
                    print("shashu") ##############################
                    reg = linregress(w,v)
                pvalmatrix[i][j]=reg.pvalue
                rvalmatrix[i][j]=reg.rvalue
                resmatrix[i][j]=reg.slope
                jsonending="_del%0db%0d.json"%(a.delay,a.delayblur)


    #yeah yeah, I know
    jn=a.f.split(".")[0]+jsonending
    jout=open(jn,"w")
    json.dump(resmatrix,jout)
    jout.close()
    jnpval=jn.replace(".json","_pval.json")
    jout=open(jnpval,"w")
    json.dump(pvalmatrix,jout)
    jout.close()
    jnrval=jn.replace(".json","_rval.json")
    jout=open(jnrval,"w")
    json.dump(rvalmatrix,jout)
    jout.close()



elif a.j!="":
    fin=open(a.j,"r")
    resmatrix=np.array(json.load(fin))
    fin.close()
    pj=a.j.replace(".json","_pval.json")
    fin=open(pj,"r")
    pvalmatrix=np.array(json.load(fin))
    fin.close()
    rj=a.j.replace(".json","_rval.json")
    fin=open(rj,"r")
    rvalmatrix=np.array(json.load(fin))
    fin.close()



for i,v in enumerate(resmatrix):
    for j,w in enumerate(v):
        if np.abs(w)<a.filter:
            resmatrix[i][j]=0.0
        if len(pvalmatrix)>0 and pvalmatrix[i][j]>=a.pvalfilter:
            resmatrix[i][j]=0.0
            rvalmatrix[i][j]=0.0



if a.checksymmetry:
    for i,v in enumerate(resmatrix):
        for j,w in enumerate(v):
            if w != resmatrix[j][i]:
                print("Assymetry found in slope matrix! %d,%d: %5.3f,  %d,%d: %5.3f"%(i,j,resmatrix[i][j],j,i,resmatrix[j][i]))





#The details

if a.det1!="":
    d1=parsesel(a.det1)
    d2=parsesel(a.det2)
    for i,v in enumerate(resmatrix):
        for j,w in enumerate(v):
            if (i+1 in d1 and j+1 in d2) or (i+1 in d2 and j+1 in d1):
                if np.abs(w)>=a.filterd:
                    print("i:%03d j:%03d %5.3f"%(i+1,j+1,w))


if a.noplot:
    sys.exit(0)




#The plot
ending=".png"
if a.delay>0:
    ending="_del%0db%0d.png"%(a.delay,a.delayblur)
base=""
if a.f!="":
    base=a.f
elif a.j!="":
    base=a.j
if a.filter>0:
    base=base+"_fil"
name=base.split(".")[0]+ending
pname=name.replace(".png","_pval.png")
rname=name.replace(".png","_rval.png")

if a.extraplot:
    if len(pvalmatrix)>0:
        plot_maps(pvalmatrix,pname,mini=0,maxi=a.pvalfilter)
    if len(rvalmatrix)>0:
        plot_maps(rvalmatrix,rname,mini=0,maxi=1) #I could set a max for this


plot_maps(resmatrix,name,mini=a.min,maxi=a.max)


        




