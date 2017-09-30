#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

fname=sys.argv[1]
tbf=float(sys.argv[2]) #time between frames
columns=int(sys.argv[3]) #not counting the first one (frame number)

fin=open(fname,"r")

x=[]
ys=[]


for i in range(columns):
    ys.append([])

for line in fin:
    fields=line.split()
    x.append(float(fields[0])*tbf)
    for i,v in enumerate(fields[1:]):
        ys[i].append(v)


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

for i,y in enumerate(ys):
    ax.plot(x,y,glyphs[i],label="Prop. "+str(i+1))

ax.legend(loc='right')


plt.show()
