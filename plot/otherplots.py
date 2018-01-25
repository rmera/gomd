#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys


#The following is not very pretty, sorry about that.
fname=sys.argv[1]
tbf=float(sys.argv[2]) #time between frames
columns=int(sys.argv[3]) #not counting the first one (frame number)
toplot=-1
if "--onlyplot" in sys.argv:
    toplot=int(sys.argv[-1])
runav=False
window=-1
if "--runav" in sys.argv:
    runav=True
    window=int(sys.argv[-1])
hist=False
cdf=False
density=True
if "--histogram" in sys.argv:
    hist=True
    if "--cdf" in sys.argv:
        cdf=True
        density=False
#End the horrible user interface



fin=open(fname,"r")

x=[]
ys=[]


if toplot==-1:
    for i in range(columns):
        ys.append([])
else:
    ys.append([])

for line in fin:
    fields=line.split()
    x.append(float(fields[0])*tbf)
    for i,v in enumerate(fields[1:]):
        if toplot==-1:
            ys[i].append(float(v))
        elif (i+1)==toplot:
            ys[0].append(float(v))


glyphs=["k-","b-","r-","g-","c-","m-","k--","b^-","ro-","g.-","c:"]        


if columns>len(glyphs):
    print "only up to ", len(glyphs), "properties can be plot simultaneously"

a = b = np.arange(0, 3, .02)
c = np.exp(a)
d = c[::-1]

# Create plots with pre-defined labels.
fig, ax = plt.subplots()

plt.xlabel('time/frame number')
plt.ylabel('Property')
if hist:
    plt.xlabel("Property")
    plt.ylabel("Population")
    
x2=x
for i,y in enumerate(ys):
    if runav:
        y=np.convolve(y, np.ones((window,))/window, mode='valid')
        ac=len(x)-len(y)
        x2=x[ac:]
    if hist:
        print(len(x),len(y))
        plt.hist(y,bins="auto",histtype="step",label="Prop. "+str(i+1),cumulative=cdf,normed=density)
        continue
    ax.plot(x2,y,glyphs[i],label="Prop. "+str(i+1))

ax.legend(loc='right')


plt.show()



