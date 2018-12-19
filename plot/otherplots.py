#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

p = argparse.ArgumentParser()

p.add_argument("fname", type=str, help="Input file, gomd format")
p.add_argument("tbf", type=float, help="delta time between frames")
p.add_argument("columns",type=int, help="Columns in the gomd files, not counting the time (i.e first) column")
tagsstring=p.add_argument("--tags", type=str, help="Tags for the Y axis and for each plotted thing,",default="")
p.add_argument("--onlyplot",type=int,help="Print only the Nth column versus the first", default=-1)
p.add_argument("--runav",type=int,help="Use a Running average with N-values window", default=-1)
p.add_argument("--histogram", type=bool,help="Plot an histogram of the values", default=False)
p.add_argument("--cdf", type=bool,help="Plot the CDF of the values, only applicable if the --histogram flag is also given"  , default=False)
p.add_argument("--tu", type=str,help="Time units" , default="ns")
p.add_argument("--forceyrange", type=str,help="Force the y range to go from two given values" , default="")


a = p.parse_args()

forcey=[]
if a.forceyrange!="":
    for i in a.forceyrange.split():
        forcey.append(float(i))


#The following is not very pretty, sorry about that.
prop="Property"
#if you don't want to tag the different numbers, just give a "" in this option.

tagslist=[]
if a.tags=="":
    for i in range(a.columns):
        tagslist.append("Prop. "+str(i+1))
elif len(a.tags.split())==1:
    prop=a.tags
    for i in range(a.columns):
        tagslist.append("")

else:
    prop=a.tags.split()[0]
    tagslist=a.tags.split()[1:]

runav=False
window=-1
if a.runav >0:
    runav=True
    window=a.runav
hist=False
cdf=False
density=True
if a.histogram:
    hist=True
    if a.cdf:
        density=False
#End the horrible user interface



fin=open(a.fname,"r")

x=[]
ys=[]


if a.onlyplot==-1:
    for i in range(a.columns):
        ys.append([])
else:
    ys.append([])

for line in fin:
    fields=line.split()
    x.append(float(fields[0])*a.tbf)
    for i,v in enumerate(fields[1:]):
        if a.onlyplot==-1:
            ys[i].append(float(v))
        elif (i+1)==a.onlyplot:
            ys[0].append(float(v))


glyphs=["b-","r-","g-","k-","c-","m-","k--","b^-","ro-","g.-","c:"]        


if a.columns>len(glyphs):
    print "only up to ", len(glyphs), "properties can be plot simultaneously"

z = b = np.arange(0, 3, .02)
c = np.exp(z)
d = c[::-1]

# Create plots with pre-defined labels.
fig, ax = plt.subplots()

if forcey!=[]:
    axes = plt.gca()
    axes.set_ylim([forcey[0],forcey[1]])

plt.xlabel('Simulation time ('+a.tu+')')
plt.ylabel(prop)
if hist:
    plt.xlabel("Property")
    plt.ylabel("Population")
    
x2=x
for i,y in enumerate(ys):
    if a.runav>0:
        y=np.convolve(y, np.ones((window,))/window, mode='valid')
        ac=len(x)-len(y)
        x2=x[ac:]
    if a.histogram:
        print(len(x),len(y))
        if tags[i]!="":
            plt.hist(y,bins="auto",histtype="step",label=tagslist[i],cumulative=cdf,normed=density)
        else:
            plt.hist(y,bins="auto",histtype="step",cumulative=cdf,normed=density)
        continue
    if tagslist!="":
        ax.plot(x2,y,glyphs[i],label=tagslist[i])
    else:
        ax.plot(x2,y,glyphs[i])
ax.legend(loc='right')


plt.show()


