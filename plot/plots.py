#!/usr/bin/env python3

#This is on the verge of what is the acceptable size for a "simple" script.
#More functionality, and I'll have to refactor the code into functions and whatnot.

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

p = argparse.ArgumentParser()

p.add_argument("fname", type=str, help="Input file, gomd format")
p.add_argument("--f2", type=str, help="Second input file, gomd format. This option will only work if both input files have one column, and both files must have the same amount of data points.", default="")
p.add_argument("--unitc", type=float, help="Unit conversion factor for _all_ the Y-axis data", default=1)

p.add_argument("-c",type=int, help="Columns in the gomd files, not counting the time (i.e first) column",default=1)
p.add_argument("--tbf", type=float, help="delta time between frames, if applicable", default=1)
p.add_argument("--xlabel", type=str, help="Label for the X axis", default="Time")
tagsstring=p.add_argument("--tags", type=str, help="Tags for the Y axis and for each plotted thing, separated by commas",default="")
p.add_argument("--onlyplot",type=int,help="Print only the Nth column versus the first", default=-1)
p.add_argument("--runav",type=int,help="Use a Running average with N-values window", default=-1)
p.add_argument("--histogram", type=bool,help="Plot an histogram of the values", default=False)
p.add_argument("--cdf", type=bool,help="Plot the CDF of the values, only applicable if the --histogram flag is also given"  , default=False)
p.add_argument("--diff", type=bool, help="Plot the difference between 2 sets of values (1st - 2nd)",default=False)

p.add_argument("--tu", type=str,help="Time units" , default="ns")
p.add_argument("--rmsf", type=bool, help="plot a rmsf run",default=False)
p.add_argument("--highlight", type=str, help="X values to highlight, give as a list of numbers separated by spaces, all surrounded by a set of quotation marks",default="")

p.add_argument("--forcexrange", type=str,help="Give two numbers separated by a space. Forces the boundary of the x axis to be twose two numbers",default="")

p.add_argument("--forceyrange", type=str,help="Give two numbers separated by a space. Forces the boundary of the y axis to be twose two numbers",default="")
a = p.parse_args()

yrange=[]
xrange=[]
force=a.forceyrange
if force!="":
    f=force.split()
    yrange.append(float(f[0]))
    yrange.append(float(f[1]))

force=a.forcexrange
if force!="":
    f=force.split()
    xrange.append(float(f[0]))
    xrange.append(float(f[1]))



#process highlight, we get the x coordinates that should be highlighted.
highstr=a.highlight.split()
xhigh=[]
for i in highstr:
    if a.highlight=="":
        break
    xhigh.append(float(i))

#The following is not very pretty, sorry about that.
prop="Property"
#if you don't want to tag the different numbers, just give a "" in this option.

tagslist=[]
if a.tags=="":
    for i in range(a.c):
        tagslist.append(str(i+1))
elif len(a.tags.split(","))==1:
    prop=a.tags
    for i in range(a.c):
        tagslist.append("")

else:
    prop=a.tags.split(",")[0]
    tagslist=a.tags.split(",")[1:]

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


extension=a.fname.split(".")[-1]
fin=open(a.fname,"r")

fins=[fin]


x=[[]]
ys=[]






if a.onlyplot==-1:
    for i in range(a.c):
        ys.append([])

else:
    ys.append([])


# Here we handle the case of a second data file
if a.f2!="":
    print("You have given 2 input files. Note that this is only supported if both files have one data column each, and the same number of data points. Both files also need to have the same extension.")
    fin2=open(a.f2,"r")
    ys.append([])
    x=[[],[]]
    if a.tags=="":
       tagslist.append(str(i+1))
    f1=a.fname.replace("."+extension,"")
    f2=a.f2.replace("."+extension,"")
    a.fname=f1+"_"+f2+"."+extension #a little Frankenstein monster.
    fins.append(fin2)


for j,fin in enumerate(fins):
    for line in fin:
        if line.startswith("@") or line.startswith("&") or line.startswith("#"):
            continue
        fields=line.split()
        x[j].append(float(fields[0])*a.tbf)
        for i,v in enumerate(fields[1:]):
            if a.onlyplot==-1:
                ys[j+i].append(float(v)*a.unitc)
            elif (i+1)==a.onlyplot:
                ys[0].append(float(v)*a.unitc)


x=x[0] #We assume that both files have the same amount of datapoints, 
#from the x coordinates we got before, we
#now get the y coordinates that are to be highlighted
#then we plot this pair with some exotic gliph
yhigh=[]
for i,val in enumerate(x):
    if a.highlight=="":
        break
    if val in xhigh:
        yhigh.append(ys[0][i]) #only the first column is highlightened

            

glyphs=["b-","r-","g-","m-","k-","c-","k--","b^-","ro-","g.-","c:"]    
histcolors=[]
for i,v in enumerate(glyphs):
	histcolors.append(v[0])


if a.c>len(glyphs):
    print("only up to ", len(glyphs), "properties can be plot simultaneously")

z = b = np.arange(0, 3, .02)
c = np.exp(z)
d = c[::-1]

# changing the rc parameters and plotting a line plot
plt.rcParams['figure.figsize'] = [10, 5]
# Create plots with pre-defined labels.
fig, ax = plt.subplots()


if  a.xlabel=="Time":
    plt.xlabel('Time '+a.tu)
else:
    plt.xlabel(a.xlabel)

if yrange:
    axes = plt.gca()
    axes.set_ylim(yrange)
if xrange:
    axes = plt.gca()
    axes.set_xlim(xrange)


if a.rmsf:
    plt.xlabel("Residue number")
    aafreq=20
    plt.xticks(np.arange(1,len(x),aafreq)) # We want more ticks in this case, as we need to pinpoint each residue
   # plt.rcParams['figure.figsize'] = [8, 4] # we need more space also




#plot differences between 2 sets.
#You have to give exactly 2 sets (1 file with 2 columns, or 2 files with 1 column)
#both sets need to have the same number of data points.
#print(len(ys[0])) #############3
if a.diff:
    newy=[]
    for i,v in enumerate(ys[0]):
        newy.append(v -ys[1][i])
      #  newy[1].append(0) ##We'll plot zeros as reference.
    ys=np.array([newy])
    plt.hlines(0, x[0], x[-1], colors="r", linestyles='dashed')
    a.fname=a.fname.replace("."+extension,"-diff."+extension) #This will be later use for the png file name.
#print(len(ys[0]))


plt.ylabel(prop)
if hist:
    plt.xlabel(prop)
    plt.ylabel("Normalized frequency")
    
x2=x


ax.plot(xhigh,yhigh,"g*")

for i,y in enumerate(ys):
    if a.runav>0:
        y=np.convolve(y, np.ones((a.runav,))/a.runav, mode='valid')
        ac=len(x)-len(y)
        x2=x[ac:]
    if a.histogram:
        print(len(x),len(y))
        if tagslist[i]!="":
            plt.hist(y,bins="auto",histtype="step",label=tagslist[i],cumulative=cdf, density=True,color=histcolors[i])
        else:
            plt.hist(y,bins="auto",histtype="step",cumulative=cdf,density=True,color=histcolors[i])
        continue
    if  tagslist[i]!="":
        ax.plot(x2,y,glyphs[i],label=tagslist[i])
    else:
        ax.plot(x2,y,glyphs[i])
ax.legend(loc='upper right')
plt.savefig(a.fname.replace("."+extension,".png"),dpi=600)

plt.show()


